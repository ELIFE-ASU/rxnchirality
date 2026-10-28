"""Map-free CIP ligand ranking on the hierarchical digraph (Rules 1a, 1b, 2).

Only atomic numbers, mass numbers and the bond graph are used; atom-map numbers
never influence a rank. Rules 3-5 (stereo-dependent) are not implemented: ties
they would break are reported as ties.

Digraph construction follows the 2013 IUPAC recommendations:

* ring closures and multiple bonds create *duplicate* nodes that carry the
  atomic number of the duplicated atom and have only phantom (Z = 0) children;
* aromatic (mancude) rings are not kekulized; instead every aromatic atom that
  takes part in a ring double bond gets one duplicate whose atomic number is the
  mean of its possible double-bond partners (so the result does not depend on
  which Kekulé structure RDKit happens to pick);
* spheres are compared breadth-first, branches explored in order of rank.

Rules are applied hierarchically: Rule 1a over the whole digraph, then 1b, then 2.
"""

from __future__ import annotations

from collections.abc import Sequence
from functools import cmp_to_key

from rdkit import Chem

from .ligands import Ligand

DEFAULT_NODE_BUDGET = 200_000

_PT = Chem.GetPeriodicTable()


class CIPBudgetExceeded(RuntimeError):
    pass


class _Node:
    __slots__ = ("atom", "children", "dist", "dup", "mass", "parent", "path", "sorted", "z")

    def __init__(self, atom, z, mass, dist, path, parent, dup):
        self.atom = atom  # atom index, or None for H/LP/phantom pseudo-nodes
        self.z = z
        self.mass = mass
        self.dist = dist  # root distance used by Rule 1b
        self.path = path  # atom indices from root to this node (real nodes only)
        self.parent = parent  # parent atom index
        self.dup = dup
        self.children = None
        self.sorted = {}


def _bond_order(bond: Chem.Bond) -> int:
    t = bond.GetBondType()
    if t == Chem.BondType.DOUBLE:
        return 2
    if t == Chem.BondType.TRIPLE:
        return 3
    return 1  # single, aromatic (handled by the mancude duplicate), dative, ...


class _Digraph:
    _KEYS = (
        lambda n: (n.z,),
        lambda n: (n.z, -n.dist),
        lambda n: (n.z, -n.dist, n.mass),
    )
    _PHANTOM_KEYS = ((0,), (0, 0), (0, 0, 0))

    def __init__(self, mol: Chem.Mol, root: int, budget: int):
        self.mol = mol
        self.root = root
        self.budget = budget
        self.n_nodes = 0
        self._mancude = self._mancude_duplicates(mol)

    @staticmethod
    def _mass(atom: Chem.Atom) -> int:
        iso = atom.GetIsotope()
        return iso if iso else _PT.GetMostCommonIsotope(atom.GetAtomicNum())

    @staticmethod
    def _mancude_duplicates(mol: Chem.Mol) -> dict[int, tuple[float, float]]:
        """(Z, mass) of the averaged duplicate for atoms with an aromatic double bond."""
        arom = [a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic()]
        if not arom:
            return {}
        kek = Chem.Mol(mol)
        try:
            Chem.Kekulize(kek, clearAromaticFlags=False)
        except Chem.KekulizeException:
            return {}
        has_double: set[int] = set()
        for b in kek.GetBonds():
            if (
                b.GetBondType() == Chem.BondType.DOUBLE
                and mol.GetBondWithIdx(b.GetIdx()).GetIsAromatic()
            ):
                has_double.update((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))
        out: dict[int, tuple[float, float]] = {}
        for idx in arom:
            if idx not in has_double:
                continue
            atom = mol.GetAtomWithIdx(idx)
            partners = [
                b.GetOtherAtom(atom)
                for b in atom.GetBonds()
                if b.GetIsAromatic() and b.GetOtherAtomIdx(idx) in has_double
            ]
            if partners:
                out[idx] = (
                    sum(p.GetAtomicNum() for p in partners) / len(partners),
                    sum(_Digraph._mass(p) for p in partners) / len(partners),
                )
        return out

    def _new(self, *args) -> _Node:
        self.n_nodes += 1
        if self.n_nodes > self.budget:
            raise CIPBudgetExceeded(f"CIP digraph exceeded {self.budget} nodes")
        return _Node(*args)

    def _real(self, idx: int, path: tuple[int, ...], parent: int | None) -> _Node:
        a = self.mol.GetAtomWithIdx(idx)
        new_path = path + (idx,)
        return self._new(idx, a.GetAtomicNum(), self._mass(a), len(path), new_path, parent, False)

    def _duplicate(self, idx: int, dist: int) -> _Node:
        a = self.mol.GetAtomWithIdx(idx)
        return self._new(None, a.GetAtomicNum(), self._mass(a), dist, (), None, True)

    def ligand_node(self, lig: Ligand) -> _Node:
        if lig.kind == "LP":
            return self._new(None, 0, 0, 1, (), None, True)
        if lig.kind == "H":
            return self._new(None, 1, lig.isotope or 1, 1, (), None, True)
        return self._real(lig.atom_idx, (self.root,), self.root)

    def children(self, node: _Node) -> list[_Node]:
        if node.children is not None:
            return node.children
        kids: list[_Node] = []
        if not node.dup and node.atom is not None:
            atom = self.mol.GetAtomWithIdx(node.atom)
            depth = len(node.path) - 1
            for bond in atom.GetBonds():
                nbr = bond.GetOtherAtomIdx(node.atom)
                order = _bond_order(bond)
                if nbr == node.parent:
                    kids.extend(self._duplicate(nbr, depth - 1) for _ in range(order - 1))
                    continue
                if nbr in node.path:  # ring closure
                    kids.append(self._duplicate(nbr, node.path.index(nbr)))
                    kids.extend(self._duplicate(nbr, depth + 1) for _ in range(order - 1))
                    continue
                kids.append(self._real(nbr, node.path, node.atom))
                kids.extend(self._duplicate(nbr, depth + 1) for _ in range(order - 1))
            for _ in range(atom.GetTotalNumHs()):
                kids.append(self._new(None, 1, 1, depth + 1, (), None, True))
            if node.atom in self._mancude:
                z, mass = self._mancude[node.atom]
                kids.append(self._new(None, z, mass, depth + 1, (), None, True))
        node.children = kids
        return kids

    def sorted_children(self, node: _Node, rule: int) -> list[_Node]:
        cached = node.sorted.get(rule)
        if cached is None:
            cached = sorted(
                self.children(node),
                key=cmp_to_key(lambda a, b: self.compare(b, a, rule)),
            )
            node.sorted[rule] = cached
        return cached

    def compare(self, a: _Node, b: _Node, rule: int) -> int:
        """Compare two branches under one rule; >0 if ``a`` ranks higher."""
        key = self._KEYS[rule]
        ka, kb = key(a), key(b)
        if ka != kb:
            return 1 if ka > kb else -1
        phantom = self._PHANTOM_KEYS[rule]
        level_a, level_b = [a], [b]
        while level_a or level_b:
            next_a: list[_Node] = []
            next_b: list[_Node] = []
            for x, y in zip(level_a, level_b):
                ca = self.sorted_children(x, rule)
                cb = self.sorted_children(y, rule)
                for i in range(max(len(ca), len(cb))):
                    kx = key(ca[i]) if i < len(ca) else phantom
                    ky = key(cb[i]) if i < len(cb) else phantom
                    if kx != ky:
                        return 1 if kx > ky else -1
                next_a.extend(ca)
                next_b.extend(cb)
            level_a, level_b = next_a, next_b
        return 0

    def compare_all_rules(self, a: _Node, b: _Node) -> int:
        for rule in range(len(self._KEYS)):
            c = self.compare(a, b, rule)
            if c:
                return c
        return 0


def rank_ligands(
    mol: Chem.Mol, center_idx: int, ligands: Sequence[Ligand], budget: int = DEFAULT_NODE_BUDGET
) -> list[int]:
    """CIP ranks of ``ligands`` around ``center_idx`` (higher = higher priority).

    Tied ligands share a rank. Raises :class:`CIPBudgetExceeded` if the digraph
    grows beyond ``budget`` nodes before the ligands are resolved.
    """
    g = _Digraph(mol, center_idx, budget)
    nodes = [g.ligand_node(lig) for lig in ligands]
    order = sorted(
        range(len(nodes)), key=cmp_to_key(lambda i, j: g.compare_all_rules(nodes[i], nodes[j]))
    )
    ranks = [0] * len(nodes)
    rank = 0
    for pos, i in enumerate(order):
        if pos and g.compare_all_rules(nodes[i], nodes[order[pos - 1]]) != 0:
            rank += 1
        ranks[i] = rank
    return ranks

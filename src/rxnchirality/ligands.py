"""Ligand tokens and oriented tetrahedral centers read from RDKit molecules.

RDKit stores a tetrahedral tag relative to the order of the atom's bonds
(``atom.GetBonds()``), with an implicit hydrogen or lone pair counted as the
*last* neighbor. Reading tags this way makes the result independent of how the
SMILES was written; ``tests/test_ligands.py`` pins the convention down.
"""

from __future__ import annotations

from dataclasses import dataclass

from rdkit import Chem

from .types import Token

LONE_PAIR = "LP"

#: Elements that can carry a stereogenic lone pair when 3-coordinate and untagged
#: (phosphines, sulfoxides, sulfonium, ...). Amines are excluded because they invert
#: rapidly; a tagged 3-coordinate N is still accepted.
LONE_PAIR_ELEMENTS = frozenset({15, 16, 33, 34, 52})

_TAG_SIGN = {
    Chem.ChiralType.CHI_TETRAHEDRAL_CCW: 1,
    Chem.ChiralType.CHI_TETRAHEDRAL_CW: -1,
}


def hydrogen_token(isotope: int = 0) -> str:
    """Hydrogens are indistinguishable, so they share one token (per isotope)."""
    return "H" if not isotope else f"[{isotope}H]"


@dataclass(frozen=True)
class Ligand:
    token: Token
    kind: str  # "atom", "H" or "LP"
    atom_idx: int | None = None  # neighbor index for "atom" (and explicit H atoms)
    isotope: int = 0  # for implicit "H"


@dataclass(frozen=True)
class CenterSide:
    """One side's view of a center: its ligands in RDKit order and, if tagged, the sign."""

    atom_idx: int
    ligands: tuple[Ligand, ...]
    sign: int | None  # +1 CCW, -1 CW, None if untagged

    @property
    def tokens(self) -> tuple[Token, ...]:
        return tuple(lig.token for lig in self.ligands)

    @property
    def tagged(self) -> bool:
        return self.sign is not None


def _neighbor_ligand(nbr: Chem.Atom, side_prefix: str) -> Ligand:
    if nbr.GetAtomicNum() == 1:
        return Ligand(hydrogen_token(nbr.GetIsotope()), "atom", nbr.GetIdx())
    m = nbr.GetAtomMapNum()
    token: Token = m if m else f"{side_prefix}:{nbr.GetIdx()}"
    return Ligand(token, "atom", nbr.GetIdx())


def tag_sign(atom: Chem.Atom) -> int | None:
    return _TAG_SIGN.get(atom.GetChiralTag())


def center_side(atom: Chem.Atom, side_prefix: str) -> CenterSide | None:
    """Build the ligand list of a potential tetrahedral center, or None.

    ``side_prefix`` ("r" or "p") makes tokens of unmapped heavy neighbors unique
    to their side. Returns None when the atom cannot be a tetrahedral center
    (wrong coordination, or two indistinguishable hydrogens).
    """
    sign = tag_sign(atom)
    ligands = [_neighbor_ligand(b.GetOtherAtom(atom), side_prefix) for b in atom.GetBonds()]
    n_h = atom.GetTotalNumHs()
    total = len(ligands) + n_h
    if total == 4 and n_h <= 1:
        if n_h:
            ligands.append(Ligand(hydrogen_token(), "H"))
    elif total == 3 and n_h == 0:
        if sign is None:
            if atom.GetAtomicNum() not in LONE_PAIR_ELEMENTS or atom.GetIsAromatic():
                return None
            if (
                atom.GetFormalCharge() == 0
                and atom.GetAtomicNum() == 16
                and not any(b.GetBondType() == Chem.BondType.DOUBLE for b in atom.GetBonds())
            ):
                return None  # neutral 3-coordinate S without S=X is not a plain center
        ligands.append(Ligand(LONE_PAIR, "LP"))
    else:
        return None
    if len({lig.token for lig in ligands}) != 4:
        return None
    return CenterSide(atom.GetIdx(), tuple(ligands), sign)

"""Mapping-aware comparison of stereocenters across a reaction."""

from __future__ import annotations

from dataclasses import dataclass

from rdkit import Chem

from . import types as T
from .cip import DEFAULT_NODE_BUDGET, CIPBudgetExceeded, rank_ligands
from .io import parse_reaction_smiles
from .ligands import CenterSide, center_side, tag_sign
from .parity import descriptor, same_handedness
from .types import CenterResult, ReactionResult, SkippedAtom, Status, sort_tokens


@dataclass
class _SideInfo:
    center: CenterSide | None  # None: atom cannot be a tetrahedral center here
    ranks: dict | None  # token -> CIP rank, None if not computed
    stereogenic: bool  # passes the stereogenicity test in use
    cip: str | None  # standard CIP label (tagged, all ranks distinct)
    note: str | None = None  # e.g. CIP budget exceeded

    @property
    def tagged(self) -> bool:
        return self.center is not None and self.center.tagged


def _side_info(mol: Chem.Mol, idx: int, prefix: str, require_stereogenic: bool, budget: int):
    center = center_side(mol.GetAtomWithIdx(idx), prefix)
    if center is None:
        return _SideInfo(None, None, False, None)
    note = None
    try:
        rank_list = rank_ligands(mol, idx, center.ligands, budget)
        ranks = dict(zip(center.tokens, rank_list))
        distinct = len(set(rank_list)) == 4
    except CIPBudgetExceeded:
        ranks, distinct, note = None, False, T.REASON_CIP_BUDGET
    stereogenic = distinct if require_stereogenic and note is None else True
    cip = descriptor(center.tokens, center.sign, ranks) if center.tagged and distinct else None
    return _SideInfo(center, ranks, stereogenic, cip, note)


def _compare(map_num: int, element: str, r: _SideInfo, p: _SideInfo) -> CenterResult:
    rt, pt = r.center.tokens, p.center.tokens
    conserved = set(rt) & set(pt)
    leaving = [t for t in rt if t not in conserved]
    entering = [t for t in pt if t not in conserved]
    common = {
        "map_num": map_num,
        "element": element,
        "reactant_cip": r.cip,
        "product_cip": p.cip,
        "conserved": sort_tokens(conserved),
        "leaving": sort_tokens(leaving),
        "entering": sort_tokens(entering),
    }
    if len(leaving) > 1:
        return CenterResult(status=Status.AMBIGUOUS, reason=T.REASON_MULTIPLE_EXCHANGED, **common)
    # Slot inheritance: the entering ligand takes the place of the leaving one.
    sub = {entering[0]: leaving[0]} if leaving else {}
    p_tokens = tuple(sub.get(t, t) for t in pt)
    same = same_handedness(rt, r.center.sign, p_tokens, p.center.sign)
    rcip = None
    reason = None
    if r.ranks is not None and len(set(r.ranks.values())) == 4:
        rcip = descriptor(p_tokens, p.center.sign, r.ranks)
    else:
        reason = r.note or T.REASON_CIP_TIE
    status = Status.RETAINED if same else Status.INVERTED
    return CenterResult(status=status, reason=reason, product_rcip=rcip, **common)


def analyze_reaction(
    rxn_smiles: str,
    *,
    require_stereogenic: bool = True,
    cip_node_budget: int = DEFAULT_NODE_BUDGET,
) -> ReactionResult:
    """Classify every mapped tetrahedral stereocenter of an atom-mapped reaction.

    A center counts as *inverted* when its ligand arrangement, read through the
    atom mapping, has the opposite handedness in the product. If exactly one
    ligand is exchanged, the entering ligand takes the leaving ligand's slot;
    with two or more exchanged the center is *ambiguous*.

    With ``require_stereogenic`` (default) a side counts as a stereocenter only
    if map-free CIP Rules 1-2 give four distinct ligand ranks, so that atom maps
    cannot turn e.g. an isopropyl CH into a "stereocenter".
    """
    rxn = parse_reaction_smiles(rxn_smiles)
    result = ReactionResult(rxn_smiles=rxn_smiles)

    def info(mol, idx, prefix):
        return _side_info(mol, idx, prefix, require_stereogenic, cip_node_budget)

    # Tagged atoms that cannot be compared: unmapped, or mapped on one side only.
    for mol, maps, other, side in (
        (rxn.reactants, rxn.reactant_maps, rxn.product_maps, "reactant"),
        (rxn.products, rxn.product_maps, rxn.reactant_maps, "product"),
    ):
        for atom in mol.GetAtoms():
            if tag_sign(atom) is None:
                continue
            m = atom.GetAtomMapNum()
            if m and m in other:
                continue
            reason = T.REASON_MAP_ONE_SIDE if m else T.REASON_UNMAPPED_CENTER
            result.skipped.append(
                SkippedAtom(side, atom.GetIdx(), m or None, atom.GetSymbol(), reason)
            )

    for m in sorted(rxn.reactant_maps.keys() & rxn.product_maps.keys()):
        ri, pi = rxn.reactant_maps[m], rxn.product_maps[m]
        r_atom = rxn.reactants.GetAtomWithIdx(ri)
        p_atom = rxn.products.GetAtomWithIdx(pi)
        r_tag, p_tag = tag_sign(r_atom) is not None, tag_sign(p_atom) is not None
        if not (r_tag or p_tag):
            continue
        element = r_atom.GetSymbol()
        r = info(rxn.reactants, ri, "r")
        p = info(rxn.products, pi, "p")
        for tagged, side_info, side, idx in ((r_tag, r, "reactant", ri), (p_tag, p, "product", pi)):
            if tagged and side_info.center is None:
                result.skipped.append(SkippedAtom(side, idx, m, element, T.REASON_INVALID_CENTER))

        r_ok = r.tagged and r.stereogenic
        p_ok = p.tagged and p.stereogenic
        if r_ok and p_ok:
            result.centers.append(_compare(m, element, r, p))
        elif r_ok:
            if p.stereogenic:
                status, reason = Status.UNSPECIFIED, T.REASON_UNTAGGED_PRODUCT
            else:
                status, reason = Status.DESTROYED, p.note or T.REASON_NOT_STEREOGENIC_PRODUCT
            result.centers.append(
                CenterResult(m, element, status, reason, reactant_cip=r.cip, product_cip=p.cip)
            )
        elif p_ok:
            if r.stereogenic:
                status, reason = Status.UNSPECIFIED, T.REASON_UNTAGGED_REACTANT
            else:
                status, reason = Status.CREATED, r.note or T.REASON_NOT_STEREOGENIC_REACTANT
            result.centers.append(
                CenterResult(m, element, status, reason, reactant_cip=r.cip, product_cip=p.cip)
            )
        else:
            for tagged, side_info, side, idx in (
                (r_tag, r, "reactant", ri),
                (p_tag, p, "product", pi),
            ):
                if tagged and side_info.center is not None:
                    reason = side_info.note or T.REASON_NOT_STEREOGENIC
                    result.skipped.append(SkippedAtom(side, idx, m, element, reason))
    return result

"""Shared test utilities: gold data loading, map-free RDKit CIP oracle, SMILES randomization."""

from __future__ import annotations

import random
from collections import defaultdict
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

DATA = Path(__file__).parent / "data"


def _none(v: str):
    return None if v == "-" else v


def load_gold() -> dict[str, dict]:
    """name -> {"rxn": smiles, "centers": {map_num: (status, r_cip, p_cip, p_rcip)}}"""
    gold: dict[str, dict] = defaultdict(lambda: {"rxn": None, "centers": {}})
    for line in (DATA / "gold_reactions.tsv").read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        name, m, status, r_cip, p_cip, p_rcip, rxn = line.split("\t")
        gold[name]["rxn"] = rxn
        if m != "-":
            gold[name]["centers"][int(m)] = (status, _none(r_cip), _none(p_cip), _none(p_rcip))
    return dict(gold)


def oracle_cip(side_smiles: str, map_num: int) -> str | None:
    """RDKit CIP label of a mapped atom, computed on a map-free re-parsed copy.

    Re-parsing matters: RDKit keeps tags on centers that only map numbers make
    stereogenic and would label them even after the maps are cleared.
    """
    mol = Chem.MolFromSmiles(side_smiles)
    target = next(a.GetIdx() for a in mol.GetAtoms() if a.GetAtomMapNum() == map_num)
    for a in mol.GetAtoms():
        a.SetAtomMapNum(0)
    smi = Chem.MolToSmiles(mol)
    order = list(mol.GetPropsAsDict(True, True)["_smilesAtomOutputOrder"])
    fresh = Chem.MolFromSmiles(smi)
    rdCIPLabeler.AssignCIPLabels(fresh)
    atom = fresh.GetAtomWithIdx(order.index(target))
    return atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else None


def randomize_side(side: str, rng: random.Random, relabel: dict[int, int]) -> str:
    """Rewrite one side: relabel maps, shuffle atom order, random SMILES, shuffle molecules."""
    pieces = []
    for frag in side.split("."):
        mol = Chem.MolFromSmiles(frag)
        for a in mol.GetAtoms():
            if a.GetAtomMapNum():
                a.SetAtomMapNum(relabel[a.GetAtomMapNum()])
        perm = list(range(mol.GetNumAtoms()))
        rng.shuffle(perm)
        mol = Chem.RenumberAtoms(mol, perm)
        pieces.append(Chem.MolToSmiles(mol, canonical=False, doRandom=True))
    rng.shuffle(pieces)
    return ".".join(pieces)


def randomize_reaction(rxn: str, rng: random.Random) -> tuple[str, dict[int, int]]:
    """Random equivalent of ``rxn``; returns it with the map relabeling used (old -> new)."""
    reactants, agents, products = rxn.split(">")
    maps = sorted(
        {a.GetAtomMapNum() for s in (reactants, products) for a in Chem.MolFromSmiles(s).GetAtoms()}
        - {0}
    )
    new = rng.sample(range(1, 10 * len(maps) + 10), len(maps))
    relabel = dict(zip(maps, new))
    return (
        f"{randomize_side(reactants, rng, relabel)}>{agents}>{randomize_side(products, rng, relabel)}",
        relabel,
    )

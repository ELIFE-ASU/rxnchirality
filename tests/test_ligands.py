"""Pin down how RDKit tags are read: bond order, implicit H / lone pair last."""

import random

import pytest
from rdkit import Chem

from rxnchirality.ligands import LONE_PAIR, center_side
from rxnchirality.parity import same_handedness


def _side(smiles: str, map_num: int):
    mol = Chem.MolFromSmiles(smiles)
    atom = next(a for a in mol.GetAtoms() if a.GetAtomMapNum() == map_num)
    return center_side(atom, "r")


def _same(s1: str, s2: str, map_num: int = 1) -> bool:
    a, b = _side(s1, map_num), _side(s2, map_num)
    return same_handedness(a.tokens, a.sign, b.tokens, b.sign)


def test_implicit_h_position_in_smiles():
    # The H sits right after the preceding atom; a leading center has it first.
    assert _same("[F:2][C@H:1]([Cl:3])[Br:4]", "[C@@H:1]([F:2])([Cl:3])[Br:4]")
    assert not _same("[F:2][C@H:1]([Cl:3])[Br:4]", "[C@H:1]([F:2])([Cl:3])[Br:4]")
    assert _same("[F:2][C@H:1]([Cl:3])[Br:4]", "[Cl:3][C@@H:1]([F:2])[Br:4]")


def test_explicit_h_matches_implicit_h():
    mol = Chem.MolFromSmiles("[F:2][C@H:1]([Cl:3])[Br:4]")
    a = center_side(mol.GetAtomWithIdx(1), "r")
    molh = Chem.AddHs(mol)
    b = center_side(molh.GetAtomWithIdx(1), "r")
    assert "H" in a.tokens and "H" in b.tokens
    assert same_handedness(a.tokens, a.sign, b.tokens, b.sign)


@pytest.mark.parametrize(
    "smiles,map_num",
    [
        ("[F:2][C@H:1]([Cl:3])[Br:4]", 1),
        ("[CH3:1][S@:2](=[O:3])[c:4]1[cH:5][cH:6][c:7]([CH3:8])[cH:9][cH:10]1", 2),
        ("[CH3:1][P@:2]([CH2:3][CH3:4])[c:5]1[cH:6][cH:7][cH:8][cH:9][cH:10]1", 2),
        ("[CH3:1][C@@:2]1([OH:3])[CH2:4][CH2:5][CH2:6][C@H:7]1[NH2:8]", 2),
        ("[CH3:1][C@@:2]1([OH:3])[CH2:4][CH2:5][CH2:6][C@H:7]1[NH2:8]", 7),
    ],
)
def test_handedness_survives_rewriting(smiles, map_num):
    rng = random.Random(7)
    mol = Chem.MolFromSmiles(smiles)
    for _ in range(40):
        perm = list(range(mol.GetNumAtoms()))
        rng.shuffle(perm)
        other = Chem.MolToSmiles(Chem.RenumberAtoms(mol, perm), canonical=False, doRandom=True)
        assert _same(smiles, other, map_num), other


def test_lone_pair_token():
    side = _side("[CH3:1][S@:2](=[O:3])[CH2:4][CH3:5]", 2)
    assert side.tokens[-1] == LONE_PAIR


def test_unmapped_neighbors_get_side_tokens():
    side = _side("C[C@H:1](F)Cl", 1)
    assert sum(isinstance(t, str) and t.startswith("r:") for t in side.tokens) == 3


@pytest.mark.parametrize(
    "smiles",
    [
        "[CH3:2][CH2:1][CH3:3]",  # two hydrogens
        "[CH3:2][CH:1]=[CH2:3]",  # sp2 carbon
        "[CH3:2][S:1][CH3:3]",  # sulfide: only two neighbors
        "[cH:2]1[cH:3][cH:1][cH:4][cH:5][cH:6]1",  # aromatic CH
    ],
)
def test_not_a_center(smiles):
    assert _side(smiles, 1) is None

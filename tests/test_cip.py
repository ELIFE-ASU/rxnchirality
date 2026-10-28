"""CIP engine (Rules 1a, 1b, 2) against RDKit's rdCIPLabeler on map-free molecules."""

import random

import pytest
from rdkit import Chem
from rdkit.Chem import rdCIPLabeler
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

from rxnchirality.cip import CIPBudgetExceeded, rank_ligands
from rxnchirality.ligands import center_side
from rxnchirality.parity import descriptor

from helpers import DATA


def _skeletons():
    for line in (DATA / "cip_molecules.smi").read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            yield line.strip()


def _our_labels(mol):
    out = {}
    for atom in mol.GetAtoms():
        side = center_side(atom, "x")
        if side is None or not side.tagged:
            continue
        ranks = dict(zip(side.tokens, rank_ligands(mol, atom.GetIdx(), side.ligands)))
        out[atom.GetIdx()] = descriptor(side.tokens, side.sign, ranks)
    return out


@pytest.mark.parametrize("skeleton", list(_skeletons()))
def test_agrees_with_rdkit(skeleton):
    opts = StereoEnumerationOptions(maxIsomers=8, unique=True, onlyUnassigned=True)
    for iso in EnumerateStereoisomers(Chem.MolFromSmiles(skeleton), options=opts):
        mol = Chem.MolFromSmiles(Chem.MolToSmiles(iso))
        rdCIPLabeler.AssignCIPLabels(mol)
        for idx, ours in _our_labels(mol).items():
            atom = mol.GetAtomWithIdx(idx)
            theirs = atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else None
            if ours is None:
                continue  # Rule 1-2 tie: only Rules 3-5 (not implemented) can resolve it
            assert ours == theirs, (Chem.MolToSmiles(mol), idx)


def test_known_rule_3_5_cases_are_ties():
    # Pseudoasymmetric center (r/s) between two enantiomorphic CHFCH3 branches.
    mol = Chem.MolFromSmiles("C[C@H](F)[C@H](F)[C@@H](C)F")
    assert _our_labels(mol)[3] is None


def test_labels_independent_of_atom_order():
    rng = random.Random(3)
    for skeleton in ["CC(N)C(=O)O", "OC(c1ccccc1)c1ccccn1", "CS(=O)c1ccc(C)cc1", "[2H]C(C)O"]:
        mol = next(iter(EnumerateStereoisomers(Chem.MolFromSmiles(skeleton))))
        mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
        ref = sorted(_our_labels(mol).values(), key=str)
        for _ in range(20):
            perm = list(range(mol.GetNumAtoms()))
            rng.shuffle(perm)
            smi = Chem.MolToSmiles(Chem.RenumberAtoms(mol, perm), canonical=False, doRandom=True)
            assert sorted(_our_labels(Chem.MolFromSmiles(smi)).values(), key=str) == ref


def test_map_numbers_do_not_create_stereocenters():
    # RDKit's own labeler calls this R because the two methyls carry different maps.
    mol = Chem.MolFromSmiles("[CH3:1][C@H:2]([CH3:3])[Cl:4]")
    side = center_side(mol.GetAtomWithIdx(1), "r")
    ranks = rank_ligands(mol, 1, side.ligands)
    assert len(set(ranks)) == 3


def test_isotopes_rule_2():
    mol = Chem.MolFromSmiles("[2H][C@@H](C)O")
    assert _our_labels(mol)[1] in ("R", "S")


def test_budget():
    cage = "C12C3C4C5C1C6C7C2C8C3C9C4C%10C5C6C%11C7C8C9C%10%11"
    mol = Chem.MolFromSmiles(f"O[C@H]({cage}){cage}")
    side = center_side(mol.GetAtomWithIdx(1), "r")
    with pytest.raises(CIPBudgetExceeded):
        rank_ligands(mol, 1, side.ligands, budget=1000)

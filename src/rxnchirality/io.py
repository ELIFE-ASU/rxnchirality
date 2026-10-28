"""Parsing and validation of atom-mapped reaction SMILES."""

from __future__ import annotations

from dataclasses import dataclass

from rdkit import Chem


class ReactionParseError(ValueError):
    pass


@dataclass
class ParsedReaction:
    """Reactant and product sides, each parsed as one (possibly multi-fragment) molecule.

    Agents (the middle field of ``r>a>p``) are ignored. ``reactant_maps`` and
    ``product_maps`` map atom-map number to atom index.
    """

    reactants: Chem.Mol
    products: Chem.Mol
    reactant_maps: dict[int, int]
    product_maps: dict[int, int]


def _parse_side(smiles: str, label: str) -> Chem.Mol:
    if not smiles:
        return Chem.Mol()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ReactionParseError(f"could not parse {label} SMILES: {smiles!r}")
    return mol


def _map_index(mol: Chem.Mol, label: str) -> dict[int, int]:
    maps: dict[int, int] = {}
    for atom in mol.GetAtoms():
        m = atom.GetAtomMapNum()
        if not m:
            continue
        if m in maps:
            raise ReactionParseError(f"atom-map number {m} occurs more than once in {label}")
        maps[m] = atom.GetIdx()
    return maps


def parse_reaction_smiles(rxn_smiles: str) -> ParsedReaction:
    """Parse ``reactants>agents>products`` (or ``reactants>>products``).

    Chiral tags are kept as written; no stereo cleanup is applied beyond RDKit's
    default sanitization. A CXSMILES extension (``|...|``) is ignored.
    """
    text = rxn_smiles.strip().split()[0] if rxn_smiles.strip() else ""
    parts = text.split(">")
    if len(parts) != 3:
        raise ReactionParseError(f"expected 'reactants>agents>products', got {rxn_smiles!r}")
    reactants = _parse_side(parts[0], "reactants")
    products = _parse_side(parts[2], "products")
    rmaps = _map_index(reactants, "reactants")
    pmaps = _map_index(products, "products")
    for m in rmaps.keys() & pmaps.keys():
        zr = reactants.GetAtomWithIdx(rmaps[m]).GetAtomicNum()
        zp = products.GetAtomWithIdx(pmaps[m]).GetAtomicNum()
        if zr != zp:
            raise ReactionParseError(f"atom-map number {m} changes element ({zr} -> {zp})")
    return ParsedReaction(reactants, products, rmaps, pmaps)

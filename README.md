# rxnchirality

Count how many stereocenters **keep** or **change handedness** in an atom-mapped
reaction, without artifacts from SMILES canonicalization or from CIP priority changes.

```console
$ pip install .
$ rxnchirality "[CH3:1][C@H:2]([CH2:3][OH:4])[Br:5].[C-:6]#[N:7]>>[CH3:1][C@@H:2]([CH2:3][OH:4])[C:6]#[N:7].[Br-:5]"
id  map_num  element  status    reactant_cip  product_rcip  product_cip  leaving  entering  reason
1   2        C        inverted  R             S             R            5        6
```

This SN2 reaction inverts the center, but its plain CIP label stays **R** on both sides
because cyanide ranks below CH₂OH while bromide ranked above it. `rxnchirality` reports
`inverted`, and the reactant-frame product label (`product_rcip`) flips to **S**.

## Why comparing labels is not enough

| Shortcut | What goes wrong |
|---|---|
| Compare RDKit `@`/`@@` tags | Tags are relative to neighbor order, which changes whenever a SMILES is written differently. |
| Compare CIP R/S | Labels follow priorities, not geometry. Hydrogenating the vinyl group of (R)-3-methylhex-1-ene gives (S)-3-methylhexane with full retention. |
| Compute CIP on mapped molecules | RDKit treats map numbers as distinguishing: `[CH3:1][C@H:2]([CH3:3])Cl` gets the label "R" and keeps its tag even after the maps are cleared. |

## How it works

**1. Mapped parity (decides retention or inversion).** Each ligand of a center is
identified by its atom-map number, plus `H` for hydrogen (hydrogens are
indistinguishable; isotopes get their own token) and `LP` for the lone pair of a
3-coordinate center. Ligands are read in RDKit's internal order, with an implicit H or
lone pair last. The ligands are then split into conserved, leaving and entering:

* nothing exchanged: the two orientations are compared directly;
* one ligand exchanged: the entering ligand takes the leaving ligand's slot
  (*slot inheritance*). Walden inversion gives `inverted` and front-side substitution
  gives `retained`. An LP ↔ O exchange (phosphine or sulfide oxidation) is handled the
  same way;
* two or more exchanged: `ambiguous`, and left out of the counts.

The center is *inverted* when the permutation parity between the two ligand orders
disagrees with the tags. No ranking is involved, so the answer cannot depend on atom
order, molecule order or the actual map numbers.

**2. Map-free CIP (`cip.py`).** The engine builds a hierarchical digraph using CIP Rules
1a, 1b and 2. It uses duplicate atoms and averaged mancude duplicates for aromatic
rings, so the result doesn't depend on which Kekulé form is chosen. It uses only atomic
numbers, masses and connectivity, never map numbers. It does two jobs:

* *Stereogenicity guard.* A side counts as a stereocenter only if its four ligands rank
  differently. Mapped "pseudo-centers" such as a mapped isopropyl CH are therefore
  skipped, and a desymmetrization is reported as `created`. Turn the guard off with
  `require_stereogenic=False` / `--no-stereo-guard`.
* *Labels.* `reactant_cip` and `product_cip` are standard R/S. `product_rcip` is the
  product read in the **reactant frame**: conserved ligands keep their reactant ranks
  and an entering ligand inherits the leaving ligand's rank. So `product_rcip !=
  reactant_cip` exactly when the center is inverted.

Rules 3–5 (stereo-dependent priorities, pseudoasymmetric r/s) are not implemented yet.
Ties that only those rules could break are reported as ties (`cip_tie_rules_1_2`), and
the parity verdict is unaffected.

## Statuses

| status | meaning |
|---|---|
| `retained` / `inverted` | tagged stereocenter on both sides; compared through the mapping |
| `created` | tagged stereocenter in the product; not stereogenic in the reactant |
| `destroyed` | tagged stereocenter in the reactant; not stereogenic in the product |
| `unspecified` | stereogenic on both sides but untagged on one |
| `ambiguous` | two or more ligands exchanged; no unique slot correspondence |

Tagged atoms that cannot be compared (unmapped, mapped on one side only, tagged but not
stereogenic) are listed in `ReactionResult.skipped`, each with a reason code.

## Python API

```python
from rxnchirality import analyze_reaction

res = analyze_reaction(
    "[CH3:1][C@H:2]([OH:3])[C@@H:4]([CH3:5])[Br:6].[N-:7]=[N+:8]=[N-:9]"
    ">>[CH3:1][C@H:2]([OH:3])[C@H:4]([CH3:5])[N:7]=[N+:8]=[N-:9].[Br-:6]"
)
res.n_inverted, res.n_retained          # (1, 1)
for c in res.centers:
    print(c.map_num, c.status, c.reactant_cip, c.product_rcip, c.leaving, c.entering)
# 2 retained S S () ()
# 4 inverted R S (6,) (7,)
res.to_dict()                           # JSON-ready
```

## Command line

```console
rxnchirality "<rxn smiles>" ["<rxn smiles>" ...]
rxnchirality -i reactions.smi            # one "rxn_smiles [id]" per line, '-' = stdin
rxnchirality -i reactions.smi --summary  # one row of counts per reaction
rxnchirality -i reactions.smi --format json
```

## Input

The input is atom-mapped reaction SMILES (`reactants>agents>products`), fully or
partly mapped. Agents are ignored, and so is a CXSMILES `|...|` suffix. Map numbers must
be unique on each side and must keep their element across the arrow.

## Development

```console
pip install -e ".[dev]"
pytest
ruff check src tests
```

The tests include:

* a hand-derived gold set of mechanisms (`tests/data/gold_reactions.tsv`);
* an independent check of CIP labels against RDKit's `rdCIPLabeler` on map-free
  molecules;
* invariance tests that rewrite every gold reaction 30 times (random SMILES, shuffled
  molecules and atom order, relabeled maps) and require identical results.

## License

MIT. See [LICENSE](LICENSE).

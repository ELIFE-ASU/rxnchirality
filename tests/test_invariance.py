"""Results must not depend on how the reaction is written or how its maps are numbered."""

import random

import pytest

from rxnchirality import analyze_reaction

from helpers import load_gold, randomize_reaction

GOLD = load_gold()
N_VARIANTS = 30


def _canonical(res, inverse: dict[int, int]):
    """Result with map numbers translated back and side-local tokens anonymized."""

    def tok(t):
        if isinstance(t, int):
            return inverse[t]
        return t.split(":")[0] + ":*" if ":" in t else t

    centers = sorted(
        (
            inverse[c.map_num],
            c.status.value,
            c.reason,
            c.reactant_cip,
            c.product_cip,
            c.product_rcip,
            tuple(sorted(map(str, map(tok, c.conserved)))),
            tuple(sorted(map(str, map(tok, c.leaving)))),
            tuple(sorted(map(str, map(tok, c.entering)))),
        )
        for c in res.centers
    )
    skipped = sorted((s.side, inverse.get(s.map_num), s.element, s.reason) for s in res.skipped)
    return centers, skipped


@pytest.mark.parametrize("name", sorted(GOLD))
@pytest.mark.parametrize("guard", [True, False])
def test_invariant_under_rewriting(name, guard):
    rxn = GOLD[name]["rxn"]
    identity = {m: m for m in range(1, 1000)}
    ref = _canonical(analyze_reaction(rxn, require_stereogenic=guard), identity)
    rng = random.Random(f"{name}-{guard}")
    for _ in range(N_VARIANTS):
        variant, relabel = randomize_reaction(rxn, rng)
        inverse = {new: old for old, new in relabel.items()}
        got = _canonical(analyze_reaction(variant, require_stereogenic=guard), inverse)
        assert got == ref, variant

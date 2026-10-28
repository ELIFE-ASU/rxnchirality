"""Permutation parity and orientation comparison.

An oriented center is a 4-tuple of ligand tokens plus a sign: +1 when, looking
from the first ligand toward the center, the other three run counterclockwise
(RDKit ``CHI_TETRAHEDRAL_CCW``, SMILES ``@``), and -1 for clockwise. Swapping any
two ligands in the tuple flips the sign, so two oriented centers over the same
token set describe the same handedness iff ``sign_a * sign_b * parity(a -> b) == +1``.
"""

from __future__ import annotations

from collections.abc import Hashable, Sequence


def permutation_parity(src: Sequence[Hashable], dst: Sequence[Hashable]) -> int:
    """Return +1 if ``dst`` is an even permutation of ``src``, -1 if odd.

    Both sequences must hold the same distinct items.
    """
    if len(src) != len(dst) or set(src) != set(dst) or len(set(src)) != len(src):
        raise ValueError(f"not a permutation of distinct items: {src!r} vs {dst!r}")
    pos = {item: i for i, item in enumerate(dst)}
    perm = [pos[item] for item in src]
    sign = 1
    seen = [False] * len(perm)
    for i in range(len(perm)):
        if seen[i]:
            continue
        j, length = i, 0
        while not seen[j]:
            seen[j] = True
            j = perm[j]
            length += 1
        if length % 2 == 0:
            sign = -sign
    return sign


def reorient(tokens: Sequence[Hashable], sign: int, new_order: Sequence[Hashable]) -> int:
    """Sign of the same center when its ligands are listed in ``new_order``."""
    return sign * permutation_parity(tokens, new_order)


def same_handedness(
    tokens_a: Sequence[Hashable], sign_a: int, tokens_b: Sequence[Hashable], sign_b: int
) -> bool:
    return reorient(tokens_a, sign_a, tokens_b) == sign_b


def descriptor(tokens: Sequence[Hashable], sign: int, rank: dict) -> str | None:
    """R/S from an oriented center and ligand ranks (higher rank = higher priority).

    Returns None if the ranks of the four ligands are not all distinct.
    """
    ranks = [rank[t] for t in tokens]
    if len(set(ranks)) != len(ranks):
        return None
    ordered = sorted(tokens, key=lambda t: rank[t])
    # Lowest priority first, then highest -> lowest. With the lowest ligand in
    # front, CCW seen from it means CW seen from the opposite side: R.
    view = [ordered[0], ordered[3], ordered[2], ordered[1]]
    return "R" if reorient(tokens, sign, view) == 1 else "S"

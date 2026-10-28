"""Result types, status values and reason codes."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from enum import Enum

#: A ligand token: an atom-map number (int) or a special token (str) such as
#: ``"H"``, ``"LP"`` or ``"r:7"`` (unmapped heavy atom, reactant atom index 7).
Token = int | str


class Status(str, Enum):
    RETAINED = "retained"
    INVERTED = "inverted"
    CREATED = "created"
    DESTROYED = "destroyed"
    AMBIGUOUS = "ambiguous"
    UNSPECIFIED = "unspecified"

    def __str__(self) -> str:
        return self.value


# Reason codes attached to centers and skipped atoms.
REASON_MULTIPLE_EXCHANGED = "multiple_ligands_exchanged"
REASON_UNTAGGED_REACTANT = "reactant_center_untagged"
REASON_UNTAGGED_PRODUCT = "product_center_untagged"
REASON_NOT_STEREOGENIC_REACTANT = "not_stereogenic_in_reactant"
REASON_NOT_STEREOGENIC_PRODUCT = "not_stereogenic_in_product"
REASON_NOT_STEREOGENIC = "tagged_but_not_stereogenic"
REASON_UNMAPPED_CENTER = "center_atom_unmapped"
REASON_MAP_ONE_SIDE = "center_mapped_on_one_side_only"
REASON_INVALID_CENTER = "tag_on_invalid_tetrahedral_center"
REASON_CIP_TIE = "cip_tie_rules_1_2"
REASON_CIP_BUDGET = "cip_digraph_budget_exceeded"


def _token_key(t: Token) -> tuple[int, int, str]:
    return (1, 0, t) if isinstance(t, str) else (0, t, "")


def sort_tokens(tokens) -> tuple[Token, ...]:
    return tuple(sorted(tokens, key=_token_key))


@dataclass(frozen=True)
class CenterResult:
    """Outcome for one mapped stereocenter.

    ``reactant_cip`` and ``product_cip`` are standard CIP labels (Rules 1-2) of
    each side on its own. ``product_rcip`` is the product label read in the
    reactant frame: conserved ligands keep their reactant ranks and an entering
    ligand inherits the rank of the ligand it replaced, so it differs from
    ``reactant_cip`` exactly when the center is inverted.
    """

    map_num: int
    element: str
    status: Status
    reason: str | None = None
    reactant_cip: str | None = None
    product_cip: str | None = None
    product_rcip: str | None = None
    conserved: tuple[Token, ...] = ()
    leaving: tuple[Token, ...] = ()
    entering: tuple[Token, ...] = ()

    def to_dict(self) -> dict:
        d = asdict(self)
        d["status"] = self.status.value
        d["conserved"] = list(self.conserved)
        d["leaving"] = list(self.leaving)
        d["entering"] = list(self.entering)
        return d


@dataclass(frozen=True)
class SkippedAtom:
    """A chiral-tagged or mapped atom that was not reported as a center."""

    side: str  # "reactant", "product" or "both"
    atom_idx: int | None
    map_num: int | None
    element: str
    reason: str

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class ReactionResult:
    rxn_smiles: str
    centers: list[CenterResult] = field(default_factory=list)
    skipped: list[SkippedAtom] = field(default_factory=list)

    def count(self, status: Status) -> int:
        return sum(1 for c in self.centers if c.status is status)

    @property
    def n_retained(self) -> int:
        return self.count(Status.RETAINED)

    @property
    def n_inverted(self) -> int:
        return self.count(Status.INVERTED)

    @property
    def n_created(self) -> int:
        return self.count(Status.CREATED)

    @property
    def n_destroyed(self) -> int:
        return self.count(Status.DESTROYED)

    @property
    def n_ambiguous(self) -> int:
        return self.count(Status.AMBIGUOUS)

    @property
    def n_unspecified(self) -> int:
        return self.count(Status.UNSPECIFIED)

    def center(self, map_num: int) -> CenterResult | None:
        for c in self.centers:
            if c.map_num == map_num:
                return c
        return None

    def summary(self) -> dict[str, int]:
        return {s.value: self.count(s) for s in Status}

    def to_dict(self) -> dict:
        return {
            "rxn_smiles": self.rxn_smiles,
            "summary": self.summary(),
            "centers": [c.to_dict() for c in self.centers],
            "skipped": [s.to_dict() for s in self.skipped],
        }

"""Command-line interface: ``rxnchirality "<rxn smiles>"`` or ``rxnchirality -i reactions.smi``."""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Iterator

from .compare import analyze_reaction
from .io import ReactionParseError
from .types import Status

CENTER_COLUMNS = [
    "id",
    "map_num",
    "element",
    "status",
    "reactant_cip",
    "product_rcip",
    "product_cip",
    "leaving",
    "entering",
    "reason",
]
SUMMARY_COLUMNS = ["id"] + [s.value for s in Status] + ["error"]


def _read_inputs(args) -> Iterator[tuple[str, str]]:
    """Yield (id, rxn_smiles). Input lines are ``rxn_smiles[<tab or space>id]``."""
    for i, rxn in enumerate(args.reaction or [], 1):
        yield str(i), rxn
    if args.input == "-":
        yield from _parse_lines(sys.stdin)
    elif args.input:
        with open(args.input) as fh:
            yield from _parse_lines(fh)


def _parse_lines(lines) -> Iterator[tuple[str, str]]:
    for n, line in enumerate(lines, 1):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(None, 1)
        yield (parts[1] if len(parts) > 1 else str(n)), parts[0]


def _fmt(value) -> str:
    if value is None:
        return ""
    if isinstance(value, (list, tuple)):
        return ",".join(str(v) for v in value)
    return str(value)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="rxnchirality",
        description="Count stereocenters that keep or change handedness in atom-mapped reactions.",
    )
    p.add_argument("reaction", nargs="*", help="atom-mapped reaction SMILES")
    p.add_argument(
        "-i",
        "--input",
        help="file with one reaction SMILES (and optional id) per line; '-' for stdin",
    )
    p.add_argument("--format", choices=["tsv", "json"], default="tsv")
    p.add_argument("--summary", action="store_true", help="one row of counts per reaction")
    p.add_argument(
        "--no-stereo-guard",
        action="store_true",
        help="accept any tagged center, even if map-free CIP Rules 1-2 find tied ligands",
    )
    args = p.parse_args(argv)
    if not args.reaction and not args.input:
        p.error("give a reaction SMILES or --input")

    out = sys.stdout
    records = []
    if args.format == "tsv":
        out.write("\t".join(SUMMARY_COLUMNS if args.summary else CENTER_COLUMNS) + "\n")
    n_errors = 0
    for rid, rxn in _read_inputs(args):
        try:
            res = analyze_reaction(rxn, require_stereogenic=not args.no_stereo_guard)
        except ReactionParseError as e:
            n_errors += 1
            if args.format == "json":
                records.append({"id": rid, "rxn_smiles": rxn, "error": str(e)})
            elif args.summary:
                out.write("\t".join([rid] + [""] * len(Status) + [str(e)]) + "\n")
            else:
                print(f"{rid}: {e}", file=sys.stderr)
            continue
        if args.format == "json":
            records.append({"id": rid, **res.to_dict()})
        elif args.summary:
            counts = res.summary()
            out.write("\t".join([rid] + [str(counts[s.value]) for s in Status] + [""]) + "\n")
        else:
            for c in res.centers:
                row = {"id": rid, **c.to_dict()}
                out.write("\t".join(_fmt(row[k]) for k in CENTER_COLUMNS) + "\n")
    if args.format == "json":
        json.dump(records, out, indent=2)
        out.write("\n")
    return 1 if n_errors else 0


if __name__ == "__main__":
    raise SystemExit(main())

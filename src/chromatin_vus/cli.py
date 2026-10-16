import argparse
import sys
from . import __version__
from .input import InputError
from .pipeline import run, demo


def main(argv=None):
    parser = argparse.ArgumentParser(description="Traceable research annotation and assay suggestions for chromatin-regulator missense SNVs.")
    parser.add_argument("--version", action="version", version=__version__)
    commands = parser.add_subparsers(dest="command", required=True)
    annotate = commands.add_parser("annotate", help="Annotate a human GRCh38 variant list using public services")
    annotate.add_argument("input")
    annotate.add_argument("--out", required=True)
    annotate.add_argument("--cache", default=".cache/chromatin-vus")
    annotate.add_argument("--panel", help="Custom panel CSV, replacing the bundled panel")
    annotate.add_argument("--offline", action="store_true", help="Replay only cached queries")
    annotate.add_argument("--refresh", action="store_true", help="Ignore cached responses")
    annotate.add_argument("--revel-threshold", type=float, default=0.5)
    annotate.add_argument("--af-threshold", type=float, default=0.001)
    annotate.add_argument("--context", choices=["neurodevelopment", "general"], default="neurodevelopment")
    annotate.add_argument("--limit", type=int, default=200)
    example = commands.add_parser("demo", help="Generate a clearly marked synthetic report without network access")
    example.add_argument("--out", required=True)
    args = parser.parse_args(argv)
    try:
        if args.command == "demo":
            manifest = demo(args.out)
        else:
            manifest = run(args.input, args.out, cache=args.cache, panel_path=args.panel, offline=args.offline, refresh=args.refresh,
                           revel_threshold=args.revel_threshold, af_threshold=args.af_threshold, context=args.context, limit=args.limit,
                           progress=lambda i, n, r: print(f"[{i}/{n}] {r['input_id']}: {r['status']}", flush=True))
        print(f"Report: {args.out}/report.html")
        print(f"Annotated {manifest['n_annotated']}/{manifest['n_input_alleles']} alleles; see provenance.json and row warnings for incomplete sources.")
        return 0 if manifest["n_annotated"] == manifest["n_input_alleles"] else 2
    except (InputError, OSError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

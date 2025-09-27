#!/usr/bin/env python3
import argparse
from . import __version__

def main():
    parser = argparse.ArgumentParser(
        description="metaProtease: A pipeline for the identification of protease genes from metagenome data"
    )
    parser.add_argument("--version", action="version", version=f"metaProtease {__version__}")

    subparsers = parser.add_subparsers(dest="command")

    # Placeholder reference subcommand
    ref = subparsers.add_parser("reference", help="Run reference-based workflow (WIP)")
    ref.set_defaults(func=lambda args: print("Reference-based workflow not yet implemented."))

    # Placeholder exploratory subcommand
    exp = subparsers.add_parser("exploratory", help="Run exploratory workflow (WIP)")
    exp.set_defaults(func=lambda args: print("Exploratory workflow not yet implemented."))

    args = parser.parse_args()
    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()


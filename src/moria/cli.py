from pathlib import Path
import argparse

from . import (
    data_prep_early,
    run_xgf_conversion,
    data_prep,
    matchup_files,
    loc_trans,
    extract_psf_1,
    extract_psf_2,
)


def main():
    parser = argparse.ArgumentParser(
        prog="moria",
        description="MORIA HST fly star reduction pipeline",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    prep = subparsers.add_parser("prep", help="Copy binaries and prepare early data")
    prep.add_argument(
        "project_dir",
        type=Path,
        help="Root directory with 00.DATA, 01.XYM, etc.",
    )

    xgf = subparsers.add_parser("xgf", help="Convert _flc files to _WJ2")
    xgf.add_argument("project_dir", type=Path)

    dp = subparsers.add_parser("data-prep", help="Prepare IN.* files")
    dp.add_argument("project_dir", type=Path)

    match = subparsers.add_parser("matchup", help="Create MATCHUP files")
    match.add_argument("project_dir", type=Path)

    loc = subparsers.add_parser("loc-trans", help="Run local transformation stage")
    loc.add_argument("project_dir", type=Path)

    psf1 = subparsers.add_parser("extract-psf-1", help="First PSF extraction stage")
    psf1.add_argument("project_dir", type=Path)

    psf2 = subparsers.add_parser("extract-psf-2", help="Second PSF extraction stage")
    psf2.add_argument("project_dir", type=Path)

    args = parser.parse_args()

    project = args.project_dir

    if args.command == "prep":
        # Use the installed package as the source of fortran binaries
        data_prep_early(source=".", destination=project)
    elif args.command == "xgf":
        run_xgf_conversion(project)
    elif args.command == "data-prep":
        data_prep(project)
    elif args.command == "matchup":
        matchup_files(project)
    elif args.command == "loc-trans":
        loc_trans(project)
    elif args.command == "extract-psf-1":
        extract_psf_1(project)
    elif args.command == "extract-psf-2":
        # placeholder; requires good_psf input in Python usage
        raise SystemExit(
            "extract-psf-2 requires good_psf selection; "
            "call moria.extract_psf_2(...) from Python."
        )


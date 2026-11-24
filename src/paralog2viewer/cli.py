# SPDX-License-Identifier: GPL-3.0-or-later
"""Command-line interface: run, build-data, check."""

import argparse
import json
import logging
import os
import shutil
import sys
from collections.abc import Sequence
from pathlib import Path

from paralog2viewer import (
    __version__,
    config,
    dataprep,
    rccx,
    runner,
)
from paralog2viewer.datalayout import DataLayout

LOGGER = logging.getLogger("paralog2viewer")


def _formats(text: str) -> tuple[str, ...]:
    try:
        return config.parse_formats(text)
    except ValueError as err:
        raise argparse.ArgumentTypeError(str(err)) from err


REQUIRED_TOOLS = (
    "minimap2",
    "samtools",
    "bcftools",
    "tabix",
    "bgzip",
    "run_clair3.sh",
    "sawfish",
    "vep",
    "blastn",
    "perl",
    "python3",
    "/opt/annovar/table_annovar.pl",
    "/opt/intervar/Intervar.py",
    "/opt/clair3_model/pileup.pt",
    "/opt/clair3_model/full_alignment.pt",
)


def _threads(text: str) -> int:
    try:
        value = int(text)
    except ValueError:
        value = 0
    if value < 1:
        raise argparse.ArgumentTypeError("--threads must be an integer >= 1")
    return value


def _memory(text: str) -> float:
    try:
        value = float(text)
    except ValueError:
        value = 0.0
    if value <= 0:
        raise argparse.ArgumentTypeError("--memory must be > 0")
    return value


def build_parser() -> argparse.ArgumentParser:
    """Builds the argument parser."""
    parser = argparse.ArgumentParser(prog="paralog2viewer")
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )
    sub = parser.add_subparsers(dest="command", required=True)

    run = sub.add_parser("run", help="process one Paraphase result")
    run.add_argument(
        "--bam",
        type=Path,
        required=True,
        help="Paraphase BAM (adjacent .bai and .json required)",
    )
    run.add_argument(
        "--reference",
        type=Path,
        help="GRCh38 FASTA with .fai and .dict (default: the copy in the "
        "data directory made by build-data)",
    )
    run.add_argument(
        "--data",
        type=Path,
        help=f"data directory from build-data "
        f"(fallback: ${config.DATA_ENV_VAR})",
    )
    run.add_argument(
        "--outdir",
        type=Path,
        default=Path("results"),
        help="output directory (default: %(default)s)",
    )
    run.add_argument("--sample", help="sample ID (default: from BAM name)")
    run.add_argument(
        "--tmp",
        type=Path,
        help="scratch directory; own subdirectory is created",
    )
    run.add_argument(
        "--format",
        dest="formats",
        type=_formats,
        default=config.VALID_FORMATS,
        help="comma list of vcf,mviewer,summary (default: all three)",
    )
    run.add_argument(
        "--figures",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="draw figures (default: on)",
    )
    run.add_argument(
        "--threads",
        type=_threads,
        help="maximum CPU threads (default: every CPU this process may use; "
        "pin with taskset); above that count it is clamped with a warning",
    )
    run.add_argument(
        "--memory",
        type=_memory,
        help="maximum memory in GB (default: 80%% of the available memory)",
    )
    run.add_argument(
        "--cache",
        type=Path,
        help="annotation cache (default: <data>/annotation_cache.sqlite)",
    )
    run.add_argument(
        "--pseudogene-tsv",
        type=Path,
        help="pseudogene allele table (default: "
        "known_gene_pseudogene_variants.tsv shipped with the package)",
    )
    run.add_argument(
        "--vcf-annotation",
        choices=config.VCF_ANNOTATION_MODES,
        default="default",
        help="VCF annotation content (default: %(default)s)",
    )
    run.add_argument(
        "--dry-run",
        action="store_true",
        help="check inputs, print the planned steps and exit",
    )
    run.add_argument(
        "--work-dir",
        type=Path,
        help="reuse this work directory to resume a run (never deleted)",
    )
    run.add_argument(
        "--resume",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="skip steps with a valid done-marker (default: on)",
    )
    run.add_argument(
        "--verify-checksums",
        action="store_true",
        help="hash every data file at the start (slow; default: sizes only)",
    )

    build = sub.add_parser("build-data", help="build the annotation data")
    build.add_argument(
        "--out", type=Path, required=True, help="output data directory"
    )
    build.add_argument(
        "--reference",
        type=Path,
        required=True,
        help="GRCh38 FASTA with .fai and .dict next to it (required)",
    )
    build.add_argument(
        "--vep-cache",
        type=Path,
        required=True,
        help="VEP cache root containing homo_sapiens_merged/ and plugins/ "
        "(plugins/ holds the VEP plugin code, alphamissense/ and "
        "maxentscan/) (required)",
    )
    build.add_argument(
        "--annovar-db",
        type=Path,
        required=True,
        help="ANNOVAR humandb directory with the hg38 annotation tables "
        "(required)",
    )
    build.add_argument(
        "--intervar-db",
        type=Path,
        required=True,
        help="InterVar intervardb directory (required)",
    )
    build.add_argument(
        "--mviewer-profile",
        type=Path,
        required=True,
        help="MViewer profile data directory containing ClinVar, Taiwan "
        "Biobank, HGMD, OMIM genemap2 and the small tables (required)",
    )
    build.add_argument(
        "--regions",
        type=Path,
        nargs="+",
        required=True,
        help="one BED file, or one or more Paraphase JSON files",
    )
    build.add_argument(
        "--threads",
        type=_threads,
        default=4,
        help="threads for copying, subsetting and BLAST (default: %(default)s)",
    )
    build.add_argument(
        "--padding",
        type=int,
        default=50000,
        help="bp added around each locus (default: %(default)s)",
    )

    check = sub.add_parser(
        "check", help="check required tools and, with --data, the data"
    )
    check.add_argument(
        "--data",
        type=Path,
        help=f"data directory to verify (fallback: ${config.DATA_ENV_VAR})",
    )
    check.add_argument(
        "--quick",
        action="store_true",
        help="compare file sizes only, skip checksums",
    )
    return parser


def check_tools(names: Sequence[str] = REQUIRED_TOOLS) -> dict[str, str | None]:
    """Looks up each tool; None marks a missing one.

    A bare name is searched on PATH; an absolute path must be a file.
    """
    return {
        name: (
            (name if os.path.isfile(name) else None)
            if name.startswith("/")
            else shutil.which(name)
        )
        for name in names
    }


def _run_check(args: argparse.Namespace) -> int:
    """Prints the tool table and the data check; returns the exit code."""
    missing = 0
    for name, path in check_tools().items():
        label = name if len(name) < 14 else f"{name} "
        print(f"{label:<14}{path or 'MISSING'}")
        missing += path is None
    data = config.resolve_data_dir(args.data)
    problems: list[str] = []
    if data is not None:
        problems = runner.data_problems(data, checksums=not args.quick)
        print(f"{'data':<14}{data}: {'ok' if not problems else 'PROBLEMS'}")
        for problem in problems:
            print(f"  {problem}")
    return 1 if missing or problems else 0


def _run_build_data(args: argparse.Namespace) -> int:
    """Builds the data directory, then the RCCX discriminating sites."""
    sources = dataprep.DataSources(
        reference=args.reference,
        vep_cache=args.vep_cache,
        annovar_humandb=args.annovar_db,
        intervar_db=args.intervar_db,
        profile=args.mviewer_profile,
    )
    dataprep.build_data(
        args.out, sources, args.regions, args.padding, args.threads
    )
    layout = DataLayout(args.out)
    layout.rccx_sites.parent.mkdir(parents=True, exist_ok=True)
    count = rccx.build_rccx_sites(
        layout.reference_fasta,
        rccx.default_genes_tsv(),
        layout.rccx_sites,
        threads=args.threads,
    )
    _register_in_manifest(layout, layout.rccx_sites)
    print(f"data directory ready: {args.out} ({count} RCCX sites)")
    return 0


def _register_in_manifest(layout: DataLayout, path: Path) -> None:
    """Adds a file made after ``build_data`` to the manifest."""
    manifest = layout.read_manifest()
    manifest["files"][path.relative_to(layout.root).as_posix()] = {
        "size": path.stat().st_size,
        "sha256": dataprep.sha256_file(path),
    }
    with open(layout.manifest, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)


def main(argv: Sequence[str] | None = None) -> int:
    """Runs the command line.

    Returns:
        Exit code: 0 on success, 1 if a check failed, a step failed or the
        run could not start; argument errors exit 2 through argparse.
    """
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    if args.command == "check":
        return _run_check(args)
    if args.command == "build-data":
        return _run_build_data(args)
    data = config.resolve_data_dir(args.data)
    if data is None:
        parser.error(f"--data is required (or set {config.DATA_ENV_VAR})")
    reference = args.reference or DataLayout(data).reference_fasta
    threads, clamped = config.clamp_threads(args.threads)
    if clamped:
        LOGGER.warning(
            "--threads %d exceeds the %d CPUs this process may use; using %d",
            args.threads,
            threads,
            threads,
        )
    try:
        settings = config.RunSettings(
            bam=args.bam,
            reference=reference,
            data=data,
            outdir=args.outdir,
            sample=args.sample,
            tmp=args.tmp,
            work_dir=args.work_dir,
            formats=args.formats,
            figures=args.figures,
            threads=threads,
            memory_gb=args.memory,
            cache=args.cache,
            pseudogene_tsv=args.pseudogene_tsv,
            vcf_annotation=args.vcf_annotation,
            dry_run=args.dry_run,
            resume=args.resume,
            verify_checksums=args.verify_checksums,
        )
    except ValueError as err:
        parser.error(str(err))
    try:
        result = runner.run_pipeline(settings)
    except runner.PipelineError as err:
        print(f"error: {err}", file=sys.stderr)
        return 1
    if result.dry_run:
        return 0
    if not result.ok:
        print("FAILED steps:", file=sys.stderr)
        for key, message in result.failures.items():
            print(f"  {key}: {message}", file=sys.stderr)
        return 1
    for path in result.outputs:
        print(path)
    return 0


if __name__ == "__main__":
    sys.exit(main())

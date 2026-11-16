# SPDX-License-Identifier: GPL-3.0-or-later
"""Run settings and parsing helpers."""

import dataclasses
import os
from pathlib import Path

DATA_ENV_VAR = "PARALOG2VIEWER_DATA"
VALID_FORMATS = ("vcf", "mviewer")
VCF_ANNOTATION_MODES = ("default", "full")
# gs4 rule: stay below 50% of 96 logical CPUs and of the memory.
DEFAULT_THREADS = 40
DEFAULT_MEMORY_GB = 800.0


def parse_formats(text: str) -> tuple[str, ...]:
    """Parses a comma-separated output format list.

    Args:
      text: For example "vcf,mviewer". Case and spaces are ignored.

    Returns:
      Unique formats in the order given.

    Raises:
      ValueError: If the list is empty or has an unknown format.
    """
    items = [x.strip().lower() for x in text.split(",") if x.strip()]
    if not items:
        raise ValueError("--format needs at least one of: vcf, mviewer")
    bad = [x for x in items if x not in VALID_FORMATS]
    if bad:
        raise ValueError(
            f"unknown --format value(s): {', '.join(bad)}; "
            f"choose from: {', '.join(VALID_FORMATS)}"
        )
    return tuple(dict.fromkeys(items))


def resolve_data_dir(
    cli_value: str | Path | None, env: dict[str, str] | None = None
) -> Path | None:
    """Returns the data directory; the command line wins over the env var."""
    if cli_value:
        return Path(cli_value)
    env = os.environ if env is None else env
    value = env.get(DATA_ENV_VAR)
    return Path(value) if value else None


def packaged_resource(name: str) -> Path:
    """Returns a data file shipped inside the package (``resources/``)."""
    return Path(__file__).parent / "resources" / name


@dataclasses.dataclass(frozen=True)
class RunSettings:
    """Settings of one `run` invocation.

    Attributes:
        bam: Paraphase BAM (``.bai`` and ``.json`` beside it).
        reference: GRCh38 FASTA with ``.fai``.
        data: Data directory made by ``build-data``.
        outdir: Parent of ``<outdir>/<sample>/``.
        sample: Sample ID; None derives it from the BAM name.
        tmp: Scratch parent; the run makes its own subdirectory in it.
        work_dir: Existing or new work directory to reuse for resuming; it
            is never deleted by the run (overrides the ``tmp`` subdirectory).
        gpu: Let Clair3 use a GPU.
        formats: Output formats (``vcf``, ``mviewer``).
        figures: Draw figures.
        threads: CPU budget in logical CPUs.
        cpus: Logical CPUs to pin the run to; None leaves the mask alone.
        max_memory_gb: Memory budget for the task estimates.
        cache: Annotation cache file; None uses
            ``<data>/annotation_cache.sqlite``.
        pseudogene_tsv: Pseudogene allele table; None uses the packaged table.
        vcf_annotation: ``default`` or ``full`` (more ANNOVAR fields).
        dry_run: Print the plan and exit.
        resume: Skip steps whose done-marker is valid.
        hap_threads: Threads given to one haplotype step.
        verify_checksums: Hash every data file at the start of the run.
    """

    bam: Path
    reference: Path
    data: Path
    outdir: Path
    sample: str | None = None
    tmp: Path | None = None
    work_dir: Path | None = None
    gpu: bool = False
    formats: tuple[str, ...] = VALID_FORMATS
    figures: bool = True
    threads: int = DEFAULT_THREADS
    cpus: tuple[int, ...] | None = None
    max_memory_gb: float = DEFAULT_MEMORY_GB
    cache: Path | None = None
    pseudogene_tsv: Path | None = None
    vcf_annotation: str = "default"
    dry_run: bool = False
    resume: bool = True
    hap_threads: int = 4
    verify_checksums: bool = False

    def __post_init__(self) -> None:
        """Validates numeric limits and choices.

        Raises:
            ValueError: On an invalid value.
        """
        if self.threads < 1 or self.hap_threads < 1:
            raise ValueError("threads must be >= 1")
        if self.max_memory_gb <= 0:
            raise ValueError("--max-memory-gb must be > 0")
        if self.vcf_annotation not in VCF_ANNOTATION_MODES:
            raise ValueError(
                f"--vcf-annotation must be one of {VCF_ANNOTATION_MODES}"
            )
        if self.cpus is not None and self.threads > len(self.cpus):
            raise ValueError(
                f"--threads {self.threads} exceeds the {len(self.cpus)} "
                "CPUs of --cpus"
            )

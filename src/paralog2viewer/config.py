# SPDX-License-Identifier: GPL-3.0-or-later
"""Run settings and parsing helpers."""

import dataclasses
import os
from pathlib import Path

DATA_ENV_VAR = "PARALOG2VIEWER_DATA"
VALID_FORMATS = ("vcf", "mviewer")
VCF_ANNOTATION_MODES = ("default", "full")
# Share of the available memory used when --memory is not given.
MEMORY_FRACTION = 0.8
# cgroup limits above this are "no limit" (v1 reports a huge number).
_NO_LIMIT_BYTES = 1 << 60
CGROUP_V2_MEMORY = Path("/sys/fs/cgroup/memory.max")
CGROUP_V1_MEMORY = Path("/sys/fs/cgroup/memory/memory.limit_in_bytes")
PROC_MEMINFO = Path("/proc/meminfo")


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


def default_threads() -> int:
    """Returns the number of CPUs this process may use (its affinity)."""
    return len(os.sched_getaffinity(0))


def clamp_threads(requested: int | None) -> tuple[int, bool]:
    """Resolves the thread budget against the inherited CPU affinity.

    Args:
        requested: Requested maximum threads; None means all usable CPUs.

    Returns:
        The budget and whether the request was clamped down to the affinity.
    """
    usable = default_threads()
    if requested is None:
        return usable, False
    return min(requested, usable), requested > usable


def _read_bytes_limit(path: Path) -> int | None:
    """Reads a cgroup byte limit; None when absent, ``max`` or huge."""
    try:
        text = path.read_text(encoding="ascii").strip()
    except OSError:
        return None
    if not text.isdigit() or int(text) >= _NO_LIMIT_BYTES:
        return None
    return int(text)


def _mem_available_bytes(meminfo: Path) -> int | None:
    """Reads ``MemAvailable`` of a meminfo file in bytes."""
    try:
        with open(meminfo, encoding="ascii") as handle:
            for line in handle:
                if line.startswith("MemAvailable:"):
                    return int(line.split()[1]) * 1024
    except (OSError, ValueError, IndexError):
        return None
    return None


def default_memory_gb(
    cgroup_v2: Path = CGROUP_V2_MEMORY,
    cgroup_v1: Path = CGROUP_V1_MEMORY,
    meminfo: Path = PROC_MEMINFO,
) -> float:
    """Returns 80% of the memory this process can use, in GB (1e9 bytes).

    The smaller of the cgroup limit (v2, else v1; ``max`` and huge values
    are ignored) and ``MemAvailable`` is used.

    Raises:
        ValueError: If no memory limit can be determined.
    """
    limits = [
        v
        for v in (
            _read_bytes_limit(cgroup_v2),
            _read_bytes_limit(cgroup_v1),
            _mem_available_bytes(meminfo),
        )
        if v is not None
    ]
    if not limits:
        raise ValueError("cannot determine memory; pass --memory")
    return MEMORY_FRACTION * min(limits) / 1e9


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
        formats: Output formats (``vcf``, ``mviewer``).
        figures: Draw figures.
        threads: Maximum CPU threads; None uses every CPU in the affinity.
        memory_gb: Maximum memory in GB; None uses 80% of the available.
        cache: Annotation cache file; None uses
            ``<data>/annotation_cache.sqlite``.
        pseudogene_tsv: Pseudogene allele table; None uses the packaged table.
        vcf_annotation: ``default`` or ``full`` (more ANNOVAR fields).
        dry_run: Print the plan and exit.
        resume: Skip steps whose done-marker is valid.
        verify_checksums: Hash every data file at the start of the run.
    """

    bam: Path
    reference: Path
    data: Path
    outdir: Path
    sample: str | None = None
    tmp: Path | None = None
    work_dir: Path | None = None
    formats: tuple[str, ...] = VALID_FORMATS
    figures: bool = True
    threads: int | None = None
    memory_gb: float | None = None
    cache: Path | None = None
    pseudogene_tsv: Path | None = None
    vcf_annotation: str = "default"
    dry_run: bool = False
    resume: bool = True
    verify_checksums: bool = False

    def __post_init__(self) -> None:
        """Resolves the budgets and validates choices.

        ``threads`` is clamped to the inherited CPU affinity and ``memory_gb``
        defaults to 80% of the available memory.

        Raises:
            ValueError: On an invalid value.
        """
        if self.threads is not None and self.threads < 1:
            raise ValueError("threads must be >= 1")
        if self.memory_gb is not None and self.memory_gb <= 0:
            raise ValueError("--memory must be > 0")
        if self.vcf_annotation not in VCF_ANNOTATION_MODES:
            raise ValueError(
                f"--vcf-annotation must be one of {VCF_ANNOTATION_MODES}"
            )
        threads, _ = clamp_threads(self.threads)
        memory = (
            default_memory_gb() if self.memory_gb is None else self.memory_gb
        )
        object.__setattr__(self, "threads", threads)
        object.__setattr__(self, "memory_gb", memory)

# SPDX-License-Identifier: GPL-3.0-or-later
"""End-to-end ``run`` workflow.

Steps:

* [0] check inputs: BAM, index, JSON, reference, BAM dictionary against the
  FASTA, sample sex, data manifest.
* [1] ``split.split_haplotypes`` (global step).
* [2-4] per haplotype, in the worker pool, longest source BAM first:
  ``realign.prepare_haplotype`` then ``call.call_haplotype``. A failing
  haplotype is recorded; the others continue.
* [5] ``merge.merge_haplotype_calls``; the calling BEDs are checked against
  the data regions.
* [6] ``annotate.annotate_alleles`` with the SQLite cache.
* [7] outputs: annotated VCF, MViewer table, read-pileup figures and, when
  the Paraphase JSON has an ``rccx`` region, the RCCX module figure.

Work directory: ``<tmp>/paralog2viewer.XXXXXX`` made by ``mkdtemp`` with
``--tmp`` (removed after a successful run, guarded), else ``./work/<sample>``.
Only final files go to ``<outdir>/<sample>/``. No symbolic links are made.
"""

from __future__ import annotations

import dataclasses
import json
import logging
import os
import re
import shutil
import tempfile
import time
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import pysam

from paralog2viewer import (
    __version__,
    annotate,
    call,
    config,
    dataprep,
    merge,
    mviewer_table,
    plot,
    rccx,
    realign,
    scheduler,
    split,
    summary,
    vcf_out,
)
from paralog2viewer.datalayout import DataLayout

LOGGER = logging.getLogger(__name__)

SAMPLE_PATTERN = re.compile(r"[A-Za-z0-9_-]+")
WORK_PREFIX = "paralog2viewer."
CACHE_NAME = "annotation_cache.sqlite"
PSEUDOGENE_FILE = "known_gene_pseudogene_variants.tsv"
RCCX_KEY = "rccx"
# One thread per haplotype task: Clair3 and Sawfish each use about one core
# on these small BAMs, so the budget is spent on many haplotypes at once.
HAP_THREADS = 1
# Memory estimates per task in GB (call: observed peak 3-5 GB, with margin;
# the others are rough guesses).
MEMORY_GB = {
    "split": 4.0,
    "prepare": 4.0,
    "call": 6.0,
    "merge": 4.0,
    "annotate": 64.0,
    "output": 8.0,
}
_JOB_PATHS = ("source_bam",)
_PREPARED_PATHS = ("bam", "bed", "expected_cn", "audit")


class PipelineError(RuntimeError):
    """The run cannot start or a global step failed."""


@dataclasses.dataclass(frozen=True)
class Inputs:
    """Validated inputs of a run.

    Attributes:
        sample: Sample ID used in output names.
        sex: ``male`` or ``female`` from the Paraphase JSON.
        json_path: Paraphase JSON beside the BAM.
        genes: Gene keys of the Paraphase JSON.
        rccx_haplotypes: True when the RCCX region lists haplotypes (a
            region without any, such as a failed Paraphase call, has
            only null fields).
    """

    sample: str
    sex: str
    json_path: Path
    genes: tuple[str, ...]
    rccx_haplotypes: bool = False

    @property
    def has_rccx(self) -> bool:
        """True when the JSON has an RCCX region with haplotypes."""
        return RCCX_KEY in self.genes and self.rccx_haplotypes


@dataclasses.dataclass
class RunResult:
    """Outcome of :func:`run_pipeline`.

    Attributes:
        sample: Sample ID.
        out_dir: ``<outdir>/<sample>``.
        work_dir: Work directory (None for a dry run).
        outputs: Published files.
        failures: Error text per failed step or haplotype.
        dry_run: True if nothing was run.
        work_removed: True if the ``--tmp`` subdirectory was deleted.
    """

    sample: str
    out_dir: Path
    work_dir: Path | None = None
    outputs: list[Path] = dataclasses.field(default_factory=list)
    failures: dict[str, str] = dataclasses.field(default_factory=dict)
    dry_run: bool = False
    work_removed: bool = False

    @property
    def ok(self) -> bool:
        """True when no step failed."""
        return not self.failures


# --- [0] input checks ------------------------------------------------------


def sample_id_from(bam: Path, sample: str | None) -> str:
    """Returns the sample ID: ``--sample`` or the BAM name without suffix.

    Raises:
        PipelineError: On a name outside ``[A-Za-z0-9_-]``, or a BAM name
            that does not end in ``.paraphase.bam`` when no ID is given.
    """
    if sample is None:
        if not bam.name.endswith(".paraphase.bam"):
            raise PipelineError("input must be a Paraphase .paraphase.bam file")
        sample = bam.name[: -len(".paraphase.bam")]
    if not SAMPLE_PATTERN.fullmatch(sample):
        raise PipelineError(f"unsafe sample name {sample!r}")
    return sample


def read_sex(data: Mapping[str, Any]) -> str:
    """Returns the one sample sex of a Paraphase JSON.

    Returns:
        ``male``, ``female`` or ``unknown`` (no region reports a sex).

    Raises:
        PipelineError: If the regions disagree or report another value.
    """
    sexes = {
        v.get("sample_sex")
        for v in data.values()
        if isinstance(v, dict) and v.get("sample_sex") is not None
    }
    if not sexes:
        LOGGER.warning("Paraphase JSON reports no sample sex")
        return "unknown"
    if len(sexes) != 1 or next(iter(sexes)) not in ("male", "female"):
        raise PipelineError(
            f"Paraphase JSON has inconsistent sample sex: "
            f"{sorted(map(str, sexes))}"
        )
    return str(next(iter(sexes)))


def check_bam_dictionary(bam: Path, reference: Path) -> None:
    """Checks that every BAM contig exists in the FASTA with equal length.

    Raises:
        PipelineError: On a missing @SQ line or a mismatch.
    """
    with (
        pysam.AlignmentFile(str(bam), "rb", check_sq=False) as handle,
        pysam.FastaFile(str(reference)) as fasta,
    ):
        bad = split.dictionary_mismatches(handle.header.to_dict(), fasta)
    if bad is not None:
        raise PipelineError(
            "BAM sequence dictionary does not match the GRCh38 FASTA"
            + (f" (e.g. {bad[0]})" if bad else " (no @SQ lines)")
        )


def data_problems(
    data: Path,
    calling_bed: Sequence[Path] | None = None,
    *,
    checksums: bool = False,
) -> list[str]:
    """Runs ``dataprep.verify_data`` and ignores the annotation cache file.

    The default cache lives in the data directory but is not part of the
    manifest, so SQLite files there are not reported as unlisted.
    """
    problems = dataprep.verify_data(
        data, calling_bed=calling_bed, checksums=checksums
    )
    return [
        p for p in problems if not p.startswith(f"unlisted file: {CACHE_NAME}")
    ]


def pseudogene_path(settings: config.RunSettings) -> Path:
    """Returns the pseudogene allele table of this run."""
    return settings.pseudogene_tsv or config.packaged_resource(PSEUDOGENE_FILE)


def check_inputs(settings: config.RunSettings) -> Inputs:
    """Step [0]: validates all inputs before any work starts.

    Raises:
        PipelineError: Listing the first problem found.
    """
    bam = settings.bam
    sample = sample_id_from(bam, settings.sample)
    json_path = bam.with_suffix(".json")
    required = {
        "BAM": bam,
        "BAM index": Path(f"{bam}.bai"),
        "Paraphase JSON": json_path,
        "reference": settings.reference,
        "reference index": Path(f"{settings.reference}.fai"),
    }
    for label, path in required.items():
        if not path.is_file():
            hint = (
                " (pass --reference, or run build-data to create the "
                "data directory copy)"
                if label.startswith("reference")
                else ""
            )
            raise PipelineError(f"{label} not found: {path}{hint}")
    for name in call.CLAIR3_MODEL_FILES:
        if not (Path(call.CLAIR3_MODEL) / name).is_file():
            raise PipelineError(
                f"Clair3 model file not found: {call.CLAIR3_MODEL}/{name} "
                "(download the model and bind-mount it there; see README)"
            )
    if settings.tmp is not None and not settings.tmp.is_dir():
        raise PipelineError(f"--tmp is not a directory: {settings.tmp}")
    if settings.figures and not pseudogene_path(settings).is_file():
        raise PipelineError(
            f"pseudogene table not found: {pseudogene_path(settings)}"
        )
    try:
        data = json.loads(json_path.read_text(encoding="utf-8"))
    except ValueError as err:
        raise PipelineError(f"unreadable Paraphase JSON: {err}") from err
    if not isinstance(data, dict):
        raise PipelineError("Paraphase JSON is not an object")
    stale = sorted(
        k
        for k, v in data.items()
        if isinstance(v, dict)
        and v.get("final_haplotypes")
        and not v.get("phase_region")
    )
    if stale:
        raise PipelineError(
            "Paraphase JSON has haplotypes without phase_region (written "
            f"by an older Paraphase version?): {', '.join(stale[:5])}"
        )
    sex = read_sex(data)
    check_bam_dictionary(bam, settings.reference)
    problems = data_problems(settings.data, checksums=settings.verify_checksums)
    if problems:
        raise PipelineError("data directory: " + "; ".join(problems[:5]))
    region = data.get(RCCX_KEY)
    has_haps = isinstance(region, dict) and bool(
        region.get("haplotype_details")
    )
    return Inputs(sample, sex, json_path, tuple(data), has_haps)


def plan_lines(settings: config.RunSettings, info: Inputs) -> list[str]:
    """Returns the planned steps for ``--dry-run``."""
    sample = info.sample
    out = settings.outdir / sample
    work = (
        settings.tmp / f"{WORK_PREFIX}XXXXXX"
        if settings.tmp
        else Path("work") / sample
    )
    lines = [
        f"sample {sample} (sex {info.sex}); {len(info.genes)} genes in JSON",
        f"budget: {settings.threads} threads, "
        f"{settings.memory_gb:.0f} GB; resume "
        f"{'on' if settings.resume else 'off'}",
        f"work directory: {work}",
        "[0] check inputs: done",
        "[1] split_haplotypes (global)",
        f"[2-4] per haplotype, longest first: prepare_haplotype "
        f"({HAP_THREADS} thread) -> call_haplotype (Clair3 + Sawfish)",
        "[5] merge_haplotype_calls; calling BEDs checked against data regions",
        "[6] annotate_alleles (cache "
        f"{settings.cache or settings.data / CACHE_NAME})",
    ]
    if "vcf" in settings.formats:
        lines.append(
            f"[7] write_annotated_vcf ({settings.vcf_annotation}) -> "
            f"{out / f'{sample}.paralog.vcf.gz'}"
        )
    if "mviewer" in settings.formats:
        lines.append(
            f"[7] write_mviewer_table -> {out / f'{sample}.paralog.tsv.gz'}"
        )
    if settings.figures:
        lines.append(f"[7] plot_regions -> {out / 'figures'}/")
        if info.has_rccx:
            lines.append(
                f"[7] classify_rccx + plot_rccx_module -> "
                f"{out / 'figures' / 'RCCX_module.jpeg'}"
            )
    if "summary" in settings.formats:
        lines.append(f"[8] write_summary -> {out / summary.SUMMARY_NAME}")
    return lines


# --- work directory --------------------------------------------------------


def make_work_dir(settings: config.RunSettings, sample: str) -> Path:
    """Creates the work directory (own mkdtemp subdirectory with --tmp)."""
    if settings.work_dir is not None:
        settings.work_dir.mkdir(parents=True, exist_ok=True)
        return settings.work_dir
    if settings.tmp is not None:
        return Path(tempfile.mkdtemp(prefix=WORK_PREFIX, dir=str(settings.tmp)))
    work = Path("work") / sample
    work.mkdir(parents=True, exist_ok=True)
    return work


def remove_work_dir(work: Path, tmp: Path, created: Path) -> bool:
    """Deletes only the subdirectory this run created inside ``tmp``.

    Args:
        work: Directory to remove.
        tmp: The ``--tmp`` parent. Nothing else in it is touched.
        created: The path returned by :func:`make_work_dir`.

    Returns:
        True if removed; False if a guard refused.
    """
    if (
        work != created
        or work.is_symlink()
        or not work.name.startswith(WORK_PREFIX)
        or work.parent.resolve() != Path(tmp).resolve()
    ):
        LOGGER.warning("refusing to remove %s", work)
        return False
    shutil.rmtree(work)
    return True


def resolve_cache_path(cache: Path | None, data: Path, work: Path) -> Path:
    """Returns the annotation cache file; falls back to a per-run file.

    The cross-run cache is disabled (with a warning) if its location is
    not writable.
    """
    path = cache or data / CACHE_NAME
    writable = os.access(path.parent, os.W_OK) and (
        not path.exists() or os.access(path, os.W_OK)
    )
    if writable:
        return path
    LOGGER.warning(
        "annotation cache %s is not writable; cache disabled for this run",
        path,
    )
    return work / CACHE_NAME


# --- task helpers ----------------------------------------------------------


def _enc(obj: Any, path_fields: Sequence[str]) -> dict | None:
    """Encodes a dataclass with Path fields as a JSON-able dict."""
    if obj is None:
        return None
    out = dataclasses.asdict(obj)
    for name in path_fields:
        out[name] = str(out[name])
    return out


def _dec(cls: type, path_fields: Sequence[str]) -> Any:
    """Returns a decoder for :func:`_enc`."""

    def decode(value: dict | None) -> Any:
        if value is None:
            return None
        return cls(
            **{k: Path(v) if k in path_fields else v for k, v in value.items()}
        )

    return decode


def _enc_jobs(jobs: Sequence[split.HaplotypeJob]) -> list[dict]:
    return [_enc(j, _JOB_PATHS) for j in jobs]


def _dec_jobs(values: Sequence[dict]) -> list[split.HaplotypeJob]:
    decode = _dec(split.HaplotypeJob, _JOB_PATHS)
    return [decode(v) for v in values]


def _paths_of(value: Any) -> list[Path]:
    """Flattens None, a path or a list of paths."""
    if value is None:
        return []
    if isinstance(value, (str, os.PathLike)):
        return [Path(value)]
    return [Path(v) for v in value]


def _str_or_list(value: Any) -> Any:
    if value is None:
        return None
    if isinstance(value, (str, os.PathLike)):
        return str(value)
    return [str(v) for v in value]


def _path_or_list(value: Any) -> Any:
    if value is None:
        return None
    if isinstance(value, str):
        return Path(value)
    return [Path(v) for v in value]


def _data_padding(data: Path) -> int | None:
    """Padding (bp) the data directory annotates around each locus."""
    try:
        return DataLayout(data).read_manifest().get("padding")
    except (OSError, ValueError):
        return None


def hap_tasks(
    jobs: Sequence[split.HaplotypeJob],
    settings: config.RunSettings,
    work: Path,
) -> list[scheduler.Task]:
    """Builds the prepare and call tasks of every haplotype.

    Priority is the source BAM size so that the longest jobs start first.
    """
    reference = settings.reference
    threads = HAP_THREADS
    max_flank = _data_padding(settings.data)
    tasks: list[scheduler.Task] = []
    for job in jobs:
        prep_key, call_key = f"prepare:{job.hp}", f"call:{job.hp}"
        size = float(job.source_bam.stat().st_size)
        tasks.append(
            scheduler.Task(
                key=prep_key,
                fn=lambda deps, job=job: realign.prepare_haplotype(
                    job, reference, work, threads, max_flank
                ),
                threads=threads,
                priority=size,
                group=job.hp,
                memory_gb=MEMORY_GB["prepare"],
                inputs=[job.source_bam, reference],
                outputs=lambda r: (
                    [] if r is None else [r.bam, r.bed, r.expected_cn]
                ),
                params={
                    "job": _enc(job, _JOB_PATHS),
                    "max_flank": max_flank,
                },
                encode=lambda r: _enc(r, _PREPARED_PATHS),
                decode=_dec(realign.PreparedHaplotype, _PREPARED_PATHS),
            )
        )
        tasks.append(
            scheduler.Task(
                key=call_key,
                fn=lambda deps, pk=prep_key: _call_one(
                    deps[pk], reference, threads
                ),
                deps=[prep_key],
                threads=threads,
                priority=size,
                group=job.hp,
                memory_gb=MEMORY_GB["call"],
                inputs=lambda deps, pk=prep_key: (
                    [reference]
                    if deps[pk] is None
                    else [
                        reference,
                        deps[pk].bam,
                        deps[pk].bed,
                        deps[pk].expected_cn,
                    ]
                ),
                outputs=_paths_of,
                params={"model": call.model_checksums(call.CLAIR3_MODEL)},
                encode=_str_or_list,
                decode=_path_or_list,
            )
        )
    return tasks


def _call_one(
    prepared: realign.PreparedHaplotype | None,
    reference: Path,
    threads: int,
) -> Path | None:
    """Calls one prepared haplotype; None if it had no reads."""
    if prepared is None:
        return None
    return call.call_haplotype(
        prepared.bam.parent,
        reference,
        threads=threads,
        model_path=call.CLAIR3_MODEL,
    )


def _merge_annotate_paths(deps: Mapping[str, Any]) -> list[Path]:
    """Returns the merged VCF and annotation table of the dependencies."""
    return [deps["merge"], deps["annotate"]]


def _output_task(
    key: str,
    fn: Any,
    deps: Sequence[str],
    inputs: Any,
    *,
    group: str,
    threads: int = 1,
    params: Mapping[str, Any] | None = None,
) -> scheduler.Task:
    """Builds a non-fatal output task whose result is a path or list."""
    return scheduler.Task(
        key=key,
        fn=fn,
        deps=deps,
        threads=threads,
        group=group,
        memory_gb=MEMORY_GB["output"],
        inputs=inputs,
        outputs=_paths_of,
        params=params or {},
        encode=_str_or_list,
        decode=_path_or_list,
    )


# --- the run ---------------------------------------------------------------


def run_pipeline(settings: config.RunSettings) -> RunResult:
    """Runs the whole workflow for one Paraphase result.

    Args:
        settings: Validated run settings.

    Returns:
        The result; ``ok`` is False if any haplotype or output step failed.

    Raises:
        PipelineError: On invalid inputs, a locked work directory, a failed
            global step, or when no haplotype produced calls.
    """
    info = check_inputs(settings)
    sample = info.sample
    out_dir = settings.outdir / sample
    if settings.dry_run:
        for line in plan_lines(settings, info):
            print(line)
        return RunResult(sample, out_dir, dry_run=True)

    work = make_work_dir(settings, sample)
    created = work
    result = RunResult(sample, out_dir, work_dir=work)
    old_tmpdir = os.environ.get("TMPDIR")
    try:
        with scheduler.work_lock(work):
            LOGGER.info(
                "budget: %d threads, %.0f GB",
                settings.threads,
                settings.memory_gb,
            )
            (work / "tmp").mkdir(exist_ok=True)
            os.environ["TMPDIR"] = str(work / "tmp")
            out_dir.mkdir(parents=True, exist_ok=True)
            started = time.monotonic()
            pool = scheduler.Pool(
                settings.threads,
                work / "markers",
                max_memory_gb=settings.memory_gb,
                resume=settings.resume,
            )
            ctx: dict[str, Any] = {}
            if "summary" in settings.formats:
                # A page left by an earlier run must not pass for this one.
                (out_dir / summary.SUMMARY_NAME).unlink(missing_ok=True)
                (out_dir / f"{summary.SUMMARY_NAME}.part").unlink(
                    missing_ok=True
                )
            try:
                _execute(settings, info, work, out_dir, pool, result, ctx)
            except Exception as err:
                result.failures.setdefault("run", str(err))
                raise
            finally:
                if "summary" in settings.formats:
                    page = _write_summary(
                        settings, info, work, out_dir, result, ctx
                    )
                    if page is not None:
                        result.outputs.append(page)
                    elif ctx.get("figures") is not None:
                        # No page: keep the figures as files instead.
                        _publish_figures(
                            ctx["figures_stage"], ctx["figures_final"]
                        )
                total_wall = time.monotonic() - started
                # Kept in the work directory for profiling, not published.
                scheduler.write_timings(
                    work / "timings.tsv", pool.timings, total_wall
                )
                timing = scheduler.summarize_timings(pool.timings, total_wall)
                LOGGER.info(
                    "run summary: wall %.1f s; CPU %.1f s executed, "
                    "%.1f s cached (stored from earlier runs)",
                    timing["total_wall_s"],
                    timing["cpu_s_executed"],
                    timing["cpu_s_cached"],
                )
    except scheduler.LockedError as err:
        raise PipelineError(str(err)) from err
    except scheduler.GlobalStepError as err:
        raise PipelineError(f"{err} (work directory kept: {work})") from err
    finally:
        if old_tmpdir is None:
            os.environ.pop("TMPDIR", None)
        else:
            os.environ["TMPDIR"] = old_tmpdir
    if settings.tmp is not None and settings.work_dir is None and result.ok:
        result.work_removed = remove_work_dir(work, settings.tmp, created)
    return result


def _split_task(
    settings: config.RunSettings, info: Inputs, work: Path
) -> scheduler.Task:
    """Builds step [1], the global split into haplotype jobs."""
    reference = settings.reference
    threads = min(4, settings.threads)
    return scheduler.Task(
        key="split",
        fn=lambda deps: split.split_haplotypes(
            settings.bam,
            info.json_path,
            reference,
            work,
            threads,
            sample=info.sample,
        ),
        threads=threads,
        memory_gb=MEMORY_GB["split"],
        inputs=[
            settings.bam,
            Path(f"{settings.bam}.bai"),
            info.json_path,
            reference,
        ],
        outputs=lambda jobs: [j.source_bam for j in jobs],
        encode=_enc_jobs,
        decode=_dec_jobs,
    )


def _execute(
    settings: config.RunSettings,
    info: Inputs,
    work: Path,
    out_dir: Path,
    pool: scheduler.Pool,
    result: RunResult,
    ctx: dict[str, Any],
) -> None:
    """Steps [1] to [7] under the pool; fills ``result``.

    ``ctx`` receives ``jobs``, ``results`` and ``figures`` for the summary,
    also when a later step raises.
    """
    jobs = pool.run([_split_task(settings, info, work)]).results["split"]
    ctx["jobs"] = jobs
    if not jobs:
        raise PipelineError("no haplotype with confident reads in the JSON")

    # [2-4] per haplotype.
    report = pool.run(hap_tasks(jobs, settings, work))
    result.failures.update(
        {
            f"haplotype {key.split(':', 1)[1]}": msg
            for key, msg in report.failures.items()
            if not msg.startswith("cancelled")
        }
    )
    good = [j for j in jobs if report.results.get(f"call:{j.hp}") is not None]
    for job in jobs:
        key = f"prepare:{job.hp}"
        if key in report.results and report.results[key] is None:
            LOGGER.warning("%s: no read overlaps the locus; skipped", job.hp)
    if not good:
        raise PipelineError(
            "no haplotype produced calls: "
            + "; ".join(f"{k}: {v}" for k, v in result.failures.items())
        )

    # [5]-[7] global graph.
    # Figures are drawn in the work directory (so a resume can reuse them),
    # copied to a sibling of ``figures/`` and swapped in after success, so
    # no figure of an earlier run survives.
    fig_final = out_dir / "figures"
    fig_dir = work / "figures"
    final = pool.run(
        _output_tasks(settings, info, work, out_dir, jobs, good, fig_dir)
    )
    ctx["results"] = final.results
    for key, message in final.failures.items():
        if not message.startswith("cancelled"):
            result.failures[key] = message
    figs_ok = settings.figures and not any(
        k in final.failures for k in ("figures", "rccx")
    )
    # With the summary page the figures live in the page (each has a
    # Download link): no ``figures/`` folder, and none left by an earlier
    # run. The page reads them from the work directory.
    in_page = "summary" in settings.formats
    swapped = figs_ok and not in_page and _publish_figures(fig_dir, fig_final)
    if in_page:
        _remove_figures(fig_final)
    ctx["figures"] = (
        fig_final
        if swapped
        else fig_dir
        if figs_ok and in_page and fig_dir.is_dir()
        else None
    )
    ctx["figures_final"] = fig_final
    ctx["figures_stage"] = fig_dir
    for key in ("vcf", "mviewer", "figures", "rccx"):
        paths = _paths_of(final.results.get(key))
        if key in ("figures", "rccx") and key in final.results:
            if not swapped:
                continue
            paths = [fig_final / p.relative_to(fig_dir) for p in paths]
        result.outputs += paths


def _remove_figures(final: Path) -> None:
    """Deletes the JPEGs of an earlier run and the folder if it is empty."""
    if not final.is_dir():
        return
    for path in [*final.glob("*.jpeg"), *final.glob("*.part")]:
        path.unlink(missing_ok=True)
    try:
        final.rmdir()
    except OSError:
        pass


def rccx_symbols(info: Inputs) -> dict[str, str]:
    """Returns HP -> CYP21 gene of the RCCX haplotypes (empty if none)."""
    if not info.has_rccx:
        return {}
    try:
        return rccx.cyp21_symbols(info.json_path)
    except (KeyError, ValueError, OSError) as err:
        LOGGER.warning("RCCX gene symbols unavailable: %s", err)
        return {}


def tnx_symbols(
    info: Inputs, layout: DataLayout, merged: Path | None, hap_bam_dir: Path
) -> dict[str, str]:
    """Returns HP -> TNXA or TNXB for the TNXB region (empty if unknown)."""
    if merged is None or not layout.rccx_sites.is_file():
        return {}
    try:
        return rccx.classify_tnx_haplotypes(
            info.json_path,
            merged,
            layout.rccx_sites,
            hap_bam_dir=hap_bam_dir,
        )
    except (KeyError, ValueError, OSError) as err:
        LOGGER.warning("TNXB region gene symbols unavailable: %s", err)
        return {}


def _rccx_roles(symbols: Mapping[str, str]) -> dict[str, str]:
    """HP -> "gene" (distal copy) or "pseudogene" (proximal copy)."""
    table = rccx.load_gene_table()
    proximal = {table.pair(f)[0].name for f in ("CYP21", "TNX")}
    return {
        hp: "pseudogene" if sym in proximal else "gene"
        for hp, sym in symbols.items()
    }


def _write_summary(
    settings: config.RunSettings,
    info: Inputs,
    work: Path,
    out_dir: Path,
    result: RunResult,
    ctx: Mapping[str, Any],
) -> Path | None:
    """Writes the HTML summary; a problem is logged, never raised."""
    try:
        jobs = ctx.get("jobs") or []
        results = ctx.get("results") or {}
        fig_dir = ctx.get("figures")
        merged, annotated = results.get("merge"), results.get("annotate")
        table = out_dir / f"{info.sample}.paralog.tsv.gz"
        if "mviewer" not in settings.formats or "mviewer" not in results:
            table = None
            if merged is not None and annotated is not None:
                table = work / "summary_table.tsv.gz"
                mviewer_table.write_mviewer_table(
                    annotated, merged, info.sample, info.sex, table
                )
        region_of = {j.hp: j.region for j in jobs}
        roles = plot.load_roles()
        hap_symbol = {}
        hap_role = {}
        for j in jobs:
            role = plot.resolve_role(j.hp, roles, j.region)
            hap_role[j.hp] = role.role
            if role.role == "gene":
                hap_symbol[j.hp] = role.symbol
        hits: list[tuple] = []
        known = pseudogene_path(settings)
        if merged is not None and known.is_file():
            symbols = rccx_symbols(info) | tnx_symbols(
                info, DataLayout(settings.data), merged, work / "hap"
            )
            hap_symbol.update(symbols)
            hap_role.update(_rccx_roles(symbols))
            hits = plot.known_hits(
                plot.read_known(known),
                plot.read_alleles(merged, set()),
                region_of,
                roles,
                symbols,
            )
        return summary.write_summary(
            out_dir / summary.SUMMARY_NAME,
            sample=info.sample,
            version=__version__,
            work_dir=work,
            failures=result.failures,
            table_path=table,
            paraphase_json=info.json_path,
            known_hits=hits,
            hap_symbol=hap_symbol,
            hap_role=hap_role,
            figures_dir=fig_dir,
        )
    except Exception as err:  # noqa: BLE001 - the summary never fails a run
        LOGGER.warning("summary page not written: %s", err)
        return None


def _output_tasks(
    settings: config.RunSettings,
    info: Inputs,
    work: Path,
    out_dir: Path,
    jobs: Sequence[split.HaplotypeJob],
    good: Sequence[split.HaplotypeJob],
    fig_dir: Path,
) -> list[scheduler.Task]:
    """Builds steps [5] to [7]: merge, region check, annotation, outputs."""
    sample = info.sample
    reference = settings.reference
    layout = DataLayout(settings.data)
    budget = settings.threads
    hap_dirs = {j.hp: work / "hap" / j.hp for j in good}
    dir_list = list(hap_dirs.values())
    cache_path = resolve_cache_path(settings.cache, settings.data, work)
    tasks = [
        scheduler.Task(
            key="merge",
            fn=lambda deps: merge.merge_haplotype_calls(
                dir_list, reference, work
            ),
            memory_gb=MEMORY_GB["merge"],
            inputs=[
                p
                for d in dir_list
                for p in (d / "calls.vcf.gz", d / "reads.bam", d / "region.bed")
            ]
            + [reference],
            outputs=lambda p: [p],
            encode=str,
            decode=Path,
        ),
        scheduler.Task(
            key="verify_regions",
            fn=lambda deps: _verify_regions(settings.data, dir_list),
            deps=["merge"],
            inputs=[d / "region.bed" for d in dir_list]
            + [layout.regions_bed, layout.manifest],
        ),
        scheduler.Task(
            key="annotate",
            fn=lambda deps: annotate.annotate_alleles(
                deps["merge"],
                reference,
                layout,
                work,
                cache_path,
                annotate.VEP_FORK,
            ),
            deps=["merge", "verify_regions"],
            threads=min(budget, annotate.VEP_FORK + annotate.ANNOVAR_THREADS),
            priority=1.0,
            memory_gb=MEMORY_GB["annotate"],
            inputs=lambda deps: [deps["merge"], layout.manifest],
            outputs=lambda p: [p],
            encode=str,
            decode=Path,
        ),
    ]
    if "vcf" in settings.formats:
        out = out_dir / f"{sample}.paralog.vcf.gz"
        tasks.append(
            _output_task(
                "vcf",
                lambda deps: vcf_out.write_annotated_vcf(
                    deps["merge"],
                    deps["annotate"],
                    out,
                    mode=settings.vcf_annotation,
                ),
                ["merge", "annotate"],
                _merge_annotate_paths,
                group="vcf",
                params={"mode": settings.vcf_annotation},
            )
        )
    if "mviewer" in settings.formats:
        table = out_dir / f"{sample}.paralog.tsv.gz"

        def write_table(deps: Mapping[str, Any]) -> Path:
            mviewer_table.write_mviewer_table(
                deps["annotate"], deps["merge"], sample, info.sex, table
            )
            return table

        tasks.append(
            _output_task(
                "mviewer",
                write_table,
                ["merge", "annotate"],
                _merge_annotate_paths,
                group="mviewer",
                params={"sex": info.sex, "sample": sample},
            )
        )
    if settings.figures:
        tasks.append(
            _figures_task(settings, info, layout, fig_dir, jobs, hap_dirs)
        )
        if info.has_rccx:
            tasks.append(_rccx_task(info, layout, work, fig_dir))
    return tasks


def _figures_task(
    settings: config.RunSettings,
    info: Inputs,
    layout: DataLayout,
    fig_dir: Path,
    jobs: Sequence[split.HaplotypeJob],
    hap_dirs: Mapping[str, Path],
) -> scheduler.Task:
    """Builds the read-pileup figure task."""
    region_of = {j.hp: j.region for j in jobs}
    threads = max(1, min(settings.threads, len(set(region_of.values()))))
    pseudo = pseudogene_path(settings)
    gene_models = layout.annovar_humandb / "hg38_refGene.txt"
    work_hap = next(iter(hap_dirs.values())).parent if hap_dirs else fig_dir

    def draw(deps: Mapping[str, Any]) -> list[Path]:
        _clear_stale_figures(fig_dir)
        return plot.plot_regions(
            deps["merge"],
            hap_dirs,
            annotate.pathogenic_allele_keys(deps["annotate"]),
            pseudo,
            settings.reference,
            fig_dir,
            workers=threads,
            region_of=region_of,
            gene_models=gene_models,
            default_transcripts=layout.default_transcript,
            hap_info=info.json_path,
            hap_symbols=rccx_symbols(info)
            | tnx_symbols(info, layout, deps["merge"], work_hap),
        )

    return _output_task(
        "figures",
        draw,
        ["merge", "annotate"],
        lambda deps: (
            _merge_annotate_paths(deps)
            + [pseudo, gene_models, layout.default_transcript, info.json_path]
            + [d / "reads.bam" for d in hap_dirs.values()]
        ),
        group="figures",
        threads=threads,
    )


def _rccx_task(
    info: Inputs, layout: DataLayout, work: Path, fig_dir: Path
) -> scheduler.Task:
    """Builds the RCCX module figure task."""

    def draw_rccx(deps: Mapping[str, Any]) -> Path | None:
        (fig_dir / "RCCX_module.jpeg").unlink(missing_ok=True)
        rows = rccx.classify_rccx(
            info.json_path,
            deps["merge"],
            layout.rccx_sites,
            hap_bam_dir=work / "hap",
        )
        if not rows:
            LOGGER.warning("RCCX: no haplotype rows; no figure")
            return None
        rccx.write_rccx_table(rows, work / "rccx_calls.tsv")
        return rccx.plot_rccx_module(rows, fig_dir / "RCCX_module.jpeg")

    return _output_task(
        "rccx",
        draw_rccx,
        ["merge"],
        lambda deps: [deps["merge"], info.json_path, layout.rccx_sites],
        group="rccx",
    )


def _clear_stale_figures(stage: Path) -> None:
    """Removes read-pileup figures left in a staging directory."""
    if stage.is_dir():
        for path in stage.glob("*.jpeg"):
            if path.name != "RCCX_module.jpeg":
                path.unlink()


def _publish_figures(stage: Path, final: Path) -> bool:
    """Makes ``final`` hold exactly the figures in ``stage``.

    Each file is copied to a ``.part`` name and moved into place with
    ``os.replace``, and figures absent from ``stage`` are removed. No
    directory is renamed, because a directory swap fails on network file
    systems that release a deleted directory late.

    Returns:
        True if the staging directory existed and was published.
    """
    if not stage.is_dir():
        return False
    legacy = final.with_name(final.name + ".tmp")
    if legacy.is_dir():
        shutil.rmtree(legacy)
    final.mkdir(parents=True, exist_ok=True)
    names = {path.name for path in stage.glob("*.jpeg")}
    for path in final.glob("*.jpeg"):
        if path.name not in names:
            path.unlink()
    for name in sorted(names):
        part = final / f"{name}.part"
        shutil.copyfile(stage / name, part)
        os.replace(part, final / name)
    return True


def _verify_regions(data: Path, hap_dirs: Sequence[Path]) -> None:
    """Fails if a calling BED lies outside the data regions."""
    beds = [d / "region.bed" for d in hap_dirs]
    problems = data_problems(data, beds)
    if problems:
        raise PipelineError(
            "calling regions outside the data directory: "
            + "; ".join(problems[:5])
        )

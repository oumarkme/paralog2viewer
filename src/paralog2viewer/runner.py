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
    vcf_out,
)
from paralog2viewer.datalayout import DataLayout

LOGGER = logging.getLogger(__name__)

SAMPLE_PATTERN = re.compile(r"[A-Za-z0-9_-]+")
WORK_PREFIX = "paralog2viewer."
CACHE_NAME = "annotation_cache.sqlite"
PSEUDOGENE_FILE = "known_gene_pseudogene_variants.tsv"
RCCX_KEY = "rccx"
# Memory estimates per task in GB (rough guesses, not measurements).
MEMORY_GB = {
    "split": 4.0,
    "prepare": 4.0,
    "call": 12.0,
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
    """

    sample: str
    sex: str
    json_path: Path
    genes: tuple[str, ...]

    @property
    def has_rccx(self) -> bool:
        """True when the JSON has an RCCX region."""
        return RCCX_KEY in self.genes


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

    Raises:
        PipelineError: If the sex is missing, unknown or inconsistent.
    """
    sexes = {
        v.get("sample_sex")
        for v in data.values()
        if isinstance(v, dict) and v.get("sample_sex") is not None
    }
    if len(sexes) != 1 or next(iter(sexes)) not in ("male", "female"):
        raise PipelineError(
            f"Paraphase JSON has missing or inconsistent sample sex: "
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
    sex = read_sex(data)
    check_bam_dictionary(bam, settings.reference)
    problems = data_problems(settings.data, checksums=settings.verify_checksums)
    if problems:
        raise PipelineError("data directory: " + "; ".join(problems[:5]))
    return Inputs(sample, sex, json_path, tuple(data))


def plan_lines(settings: config.RunSettings, info: Inputs) -> list[str]:
    """Returns the planned steps for ``--dry-run``."""
    sample = info.sample
    out = settings.outdir / sample
    work = (
        settings.tmp / f"{WORK_PREFIX}XXXXXX"
        if settings.tmp
        else Path("work") / sample
    )
    pin = (
        f"pinned to {len(settings.cpus)} CPUs"
        if settings.cpus
        else "affinity unchanged"
    )
    lines = [
        f"sample {sample} (sex {info.sex}); {len(info.genes)} genes in JSON",
        f"budget: {settings.threads} threads ({pin}), "
        f"{settings.max_memory_gb:g} GB; resume "
        f"{'on' if settings.resume else 'off'}",
        f"work directory: {work}",
        "[0] check inputs: done",
        "[1] split_haplotypes (global)",
        f"[2-4] per haplotype, longest first: prepare_haplotype "
        f"({settings.hap_threads} threads) -> call_haplotype "
        f"(Clair3 + Sawfish, gpu={settings.gpu})",
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
            f"[7] write_mviewer_table -> {out / f'{sample}.paralog.mviewer.gz'}"
        )
    if settings.figures:
        lines.append(f"[7] plot_regions -> {out / 'figures'}/")
        if info.has_rccx:
            lines.append(
                f"[7] classify_rccx + plot_rccx_module -> "
                f"{out / 'figures' / 'RCCX_module.jpeg'}"
            )
    lines.append(f"timings -> {out / 'timings.tsv'}")
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


def hap_tasks(
    jobs: Sequence[split.HaplotypeJob],
    settings: config.RunSettings,
    work: Path,
) -> list[scheduler.Task]:
    """Builds the prepare and call tasks of every haplotype.

    Priority is the source BAM size so that the longest jobs start first.
    """
    reference = settings.reference
    threads = min(settings.hap_threads, settings.threads)
    tasks: list[scheduler.Task] = []
    for job in jobs:
        prep_key, call_key = f"prepare:{job.hp}", f"call:{job.hp}"
        size = float(job.source_bam.stat().st_size)
        tasks.append(
            scheduler.Task(
                key=prep_key,
                fn=lambda deps, job=job: realign.prepare_haplotype(
                    job, reference, work, threads
                ),
                threads=threads,
                priority=size,
                group=job.hp,
                memory_gb=MEMORY_GB["prepare"],
                inputs=[job.source_bam, reference],
                outputs=lambda r: (
                    [] if r is None else [r.bam, r.bed, r.expected_cn]
                ),
                params={"job": _enc(job, _JOB_PATHS)},
                encode=lambda r: _enc(r, _PREPARED_PATHS),
                decode=_dec(realign.PreparedHaplotype, _PREPARED_PATHS),
            )
        )
        tasks.append(
            scheduler.Task(
                key=call_key,
                fn=lambda deps, pk=prep_key: _call_one(
                    deps[pk], reference, threads, settings.gpu
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
                params={"gpu": settings.gpu, "model": call.CLAIR3_MODEL},
                encode=_str_or_list,
                decode=_path_or_list,
            )
        )
    return tasks


def _call_one(
    prepared: realign.PreparedHaplotype | None,
    reference: Path,
    threads: int,
    gpu: bool,
) -> Path | None:
    """Calls one prepared haplotype; None if it had no reads."""
    if prepared is None:
        return None
    return call.call_haplotype(
        prepared.bam.parent, reference, threads=threads, gpu=gpu
    )


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
            cpus = scheduler.apply_affinity(settings.cpus)
            LOGGER.info("CPU affinity: %d CPUs", len(cpus))
            if len(cpus) < settings.threads:
                LOGGER.warning(
                    "threads %d exceed the %d CPUs available",
                    settings.threads,
                    len(cpus),
                )
            (work / "tmp").mkdir(exist_ok=True)
            os.environ["TMPDIR"] = str(work / "tmp")
            out_dir.mkdir(parents=True, exist_ok=True)
            started = time.monotonic()
            pool = scheduler.Pool(
                settings.threads,
                work / "markers",
                max_memory_gb=settings.max_memory_gb,
                resume=settings.resume,
            )
            try:
                _execute(settings, info, work, out_dir, pool, result)
            finally:
                total_wall = time.monotonic() - started
                scheduler.write_timings(
                    out_dir / "timings.tsv", pool.timings, total_wall
                )
                summary = scheduler.summarize_timings(pool.timings, total_wall)
                LOGGER.info(
                    "run summary: wall %.1f s; CPU %.1f s executed, "
                    "%.1f s cached (stored from earlier runs)",
                    summary["total_wall_s"],
                    summary["cpu_s_executed"],
                    summary["cpu_s_cached"],
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


def _execute(
    settings: config.RunSettings,
    info: Inputs,
    work: Path,
    out_dir: Path,
    pool: scheduler.Pool,
    result: RunResult,
) -> None:
    """Steps [1] to [7] under the pool; fills ``result``."""
    sample = info.sample
    reference = settings.reference
    layout = DataLayout(settings.data)
    budget = settings.threads

    # [1] split (global step).
    split_threads = min(4, budget)
    split_task = scheduler.Task(
        key="split",
        fn=lambda deps: split.split_haplotypes(
            settings.bam,
            info.json_path,
            reference,
            work,
            split_threads,
            sample=sample,
        ),
        threads=split_threads,
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
    jobs = pool.run([split_task]).results["split"]
    if not jobs:
        raise PipelineError("no haplotype with confident reads in the JSON")
    region_of = {j.hp: j.region for j in jobs}

    # [2-4] per haplotype.
    report = pool.run(hap_tasks(jobs, settings, work))
    failed = {
        key.split(":", 1)[1]: msg
        for key, msg in report.failures.items()
        if not msg.startswith("cancelled")
    }
    result.failures.update({f"haplotype {hp}": m for hp, m in failed.items()})
    good = [
        j
        for j in jobs
        if f"call:{j.hp}" in report.results
        and report.results[f"call:{j.hp}"] is not None
    ]
    for job in jobs:
        if (
            f"prepare:{job.hp}" in report.results
            and report.results[f"prepare:{job.hp}"] is None
        ):
            LOGGER.warning("%s: no read overlaps the locus; skipped", job.hp)
    if not good:
        raise PipelineError(
            "no haplotype produced calls: "
            + "; ".join(f"{k}: {v}" for k, v in result.failures.items())
        )
    hap_dirs = {j.hp: work / "hap" / j.hp for j in good}
    dir_list = list(hap_dirs.values())

    # [5]-[7] global graph.
    annotate_threads = max(1, min(8, budget // 2))
    # Figures are drawn in the work directory (so a resume can reuse them),
    # copied to a sibling of ``figures/`` and swapped in after success, so
    # no figure of an earlier run survives.
    fig_final = out_dir / "figures"
    fig_dir = work / "figures"
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
                annotate_threads,
            ),
            deps=["merge", "verify_regions"],
            threads=min(budget, 2 * annotate_threads),
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
                lambda deps, out=out: vcf_out.write_annotated_vcf(
                    deps["merge"],
                    deps["annotate"],
                    out,
                    mode=settings.vcf_annotation,
                ),
                ["merge", "annotate"],
                lambda deps: [deps["merge"], deps["annotate"]],
                group="vcf",
                params={"mode": settings.vcf_annotation},
            )
        )
    if "mviewer" in settings.formats:
        out = out_dir / f"{sample}.paralog.mviewer.gz"

        def write_table(deps: Mapping[str, Any], out: Path = out) -> Path:
            mviewer_table.write_mviewer_table(
                deps["annotate"], deps["merge"], sample, info.sex, out
            )
            return out

        tasks.append(
            _output_task(
                "mviewer",
                write_table,
                ["merge", "annotate"],
                lambda deps: [deps["merge"], deps["annotate"]],
                group="mviewer",
                params={"sex": info.sex, "sample": sample},
            )
        )
    if settings.figures:
        fig_threads = min(4, budget)
        pseudo = pseudogene_path(settings)

        def draw(deps: Mapping[str, Any]) -> list[Path]:
            _clear_stale_figures(fig_dir)
            return plot.plot_regions(
                deps["merge"],
                hap_dirs,
                annotate.pathogenic_allele_keys(deps["annotate"]),
                pseudo,
                reference,
                fig_dir,
                workers=fig_threads,
                region_of=region_of,
            )

        tasks.append(
            _output_task(
                "figures",
                draw,
                ["merge", "annotate"],
                lambda deps: (
                    [deps["merge"], deps["annotate"], pseudo]
                    + [d / "reads.bam" for d in dir_list]
                ),
                group="figures",
                threads=fig_threads,
            )
        )
        if info.has_rccx:

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

            tasks.append(
                _output_task(
                    "rccx",
                    draw_rccx,
                    ["merge"],
                    lambda deps: [
                        deps["merge"],
                        info.json_path,
                        layout.rccx_sites,
                    ],
                    group="rccx",
                )
            )
    final = pool.run(tasks)
    for key, message in final.failures.items():
        if not message.startswith("cancelled"):
            result.failures[key] = message
    figure_keys = [k for k in ("figures", "rccx") if k in final.results]
    swapped = (
        settings.figures
        and not any(k in final.failures for k in ("figures", "rccx"))
        and _publish_figures(fig_dir, fig_final)
    )
    for key in ("vcf", "mviewer", "figures", "rccx"):
        paths = _paths_of(final.results.get(key))
        if swapped and key in figure_keys:
            paths = [fig_final / p.relative_to(fig_dir) for p in paths]
        result.outputs += paths


def _clear_stale_figures(stage: Path) -> None:
    """Removes read-pileup figures left in a staging directory."""
    if stage.is_dir():
        for path in stage.glob("*.jpeg"):
            if path.name != "RCCX_module.jpeg":
                path.unlink()


def _publish_figures(stage: Path, final: Path) -> bool:
    """Replaces ``final`` by a copy of ``stage`` (copy, then swap).

    Returns:
        True if the staging directory existed and was published.
    """
    if not stage.is_dir():
        return False
    temporary = final.with_name(final.name + ".tmp")
    if temporary.exists():
        shutil.rmtree(temporary)
    shutil.copytree(stage, temporary)
    if final.exists():
        shutil.rmtree(final)
    temporary.rename(final)
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

# SPDX-License-Identifier: GPL-3.0-or-later
"""Dependency-aware worker pool with a thread and memory budget.

Every task runs in its own forked process, so its wall time, CPU time and
peak resident memory can be measured per task (``RUSAGE_SELF`` plus
``RUSAGE_CHILDREN`` of the worker cover the external tools it starts). The
parent process stays single-threaded, which keeps ``fork`` safe.

Resource limits:

* CPU: the sum of ``Task.threads`` of running tasks never exceeds the pool
  budget (a counting semaphore of thread weights, kept by the parent).
  Children inherit the CPU affinity of the process; nothing is pinned here.
* Memory: the sum of ``Task.memory_gb`` estimates of running tasks never
  exceeds the budget. The numbers are estimates, not measurements.

Resume: a finished task leaves ``<markers>/<key>.json`` with a hash of its
name, input fingerprints (content hash below 1 MiB, otherwise size and mtime),
parameters and the hashes of its dependencies. A rerun skips a task whose
hash is unchanged and whose outputs still exist.

Failure policy: a failing task with a ``group`` (for example one haplotype)
is recorded, its dependents are cancelled, other groups continue. A failing
task without a group (a global step) terminates the pool and raises
:class:`GlobalStepError`.
"""

from __future__ import annotations

import contextlib
import dataclasses
import datetime
import fcntl
import hashlib
import json
import logging
import multiprocessing
import multiprocessing.connection
import os
import re
import resource
import signal
import time
import traceback
from collections.abc import Callable, Iterator, Mapping, Sequence
from pathlib import Path
from typing import Any

from paralog2viewer import __version__

LOGGER = logging.getLogger(__name__)

# Files up to this size are fingerprinted by content, larger ones by
# size and modification time.
SMALL_FILE_BYTES = 1 << 20
TIMING_COLUMNS = (
    "step",
    "group",
    "status",
    "threads",
    "wall_s",
    "cpu_s",
    "peak_rss_mb",
    "finished",
)

PathsOrFn = Sequence[Path] | Callable[[Mapping[str, Any]], Sequence[Path]]


class SchedulerError(RuntimeError):
    """Invalid task graph or an unusable pool setting."""


class LockedError(RuntimeError):
    """Another run holds the lock of the work directory."""


@dataclasses.dataclass
class Task:
    """One unit of work.

    Attributes:
        key: Unique step name, for example ``call:HP1``.
        fn: Called in the worker as ``fn(dep_results)`` with the results of
            the dependencies by key. Must return a picklable value.
        deps: Keys of tasks that must finish first.
        threads: Thread weight counted against the CPU budget (>= 1).
        priority: Larger runs earlier among ready tasks (longest first).
        group: Failure group (for example the HP value); None marks a global
            step whose failure stops the run.
        memory_gb: Estimated peak memory, counted against the budget.
        inputs: Files that define the resume hash, or a function of the
            dependency results returning them. Directories are walked.
        outputs: Files that must exist for a marker to be honored, or a
            function of the (decoded) task result.
        params: JSON-serializable parameters that enter the resume hash.
        encode: Turns the result into JSON for the marker (default: as is).
        decode: Inverse of ``encode`` (default: as is).
    """

    key: str
    fn: Callable[[Mapping[str, Any]], Any]
    deps: Sequence[str] = ()
    threads: int = 1
    priority: float = 0.0
    group: str | None = None
    memory_gb: float = 0.0
    inputs: PathsOrFn = ()
    outputs: PathsOrFn = ()
    params: Mapping[str, Any] = dataclasses.field(default_factory=dict)
    encode: Callable[[Any], Any] = lambda x: x
    decode: Callable[[Any], Any] = lambda x: x


@dataclasses.dataclass(frozen=True)
class Timing:
    """Resources used by one task.

    Attributes:
        step: Task key.
        group: Failure group or empty.
        status: ``done`` (executed in this run), ``cached`` (resumed from a
            marker; values are those stored by the run that executed it),
            ``failed`` or ``cancelled``.
        threads: Thread weight.
        wall_s: Wall-clock seconds; ``None`` if unknown (old marker).
        cpu_s: User plus system CPU seconds of the worker and its children;
            ``None`` if unknown.
        peak_rss_mb: Largest resident set of the worker or of its largest
            child, in MiB (not the sum over children); ``None`` if unknown.
        finished: UTC finish time (ISO 8601, seconds) of the executing run,
            or empty if unknown.
    """

    step: str
    group: str
    status: str
    threads: int
    wall_s: float | None = 0.0
    cpu_s: float | None = 0.0
    peak_rss_mb: float | None = 0.0
    finished: str = ""


@dataclasses.dataclass
class Report:
    """Outcome of one :meth:`Pool.run` call.

    Attributes:
        results: Result per finished or skipped task key.
        failures: Error text per failed key; cancelled dependents are
            listed with the reason.
        skipped: Keys skipped through resume markers.
    """

    results: dict[str, Any] = dataclasses.field(default_factory=dict)
    failures: dict[str, str] = dataclasses.field(default_factory=dict)
    skipped: list[str] = dataclasses.field(default_factory=list)


class GlobalStepError(RuntimeError):
    """A task without a failure group failed.

    Attributes:
        key: The failed task.
        report: Partial report of the run.
    """

    def __init__(self, key: str, message: str, report: Report) -> None:
        """Stores the failed key, its message and the partial report."""
        super().__init__(f"{key} failed: {message}")
        self.key = key
        self.report = report


@contextlib.contextmanager
def work_lock(work_dir: Path) -> Iterator[Path]:
    """Holds ``<work_dir>/run.lock`` so two runs cannot share a work dir.

    Raises:
        LockedError: If another process holds the lock.
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    path = work_dir / "run.lock"
    with open(path, "a+", encoding="utf-8") as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as err:
            raise LockedError(
                f"work directory {work_dir} is in use (lock {path})"
            ) from err
        handle.seek(0)
        handle.truncate()
        handle.write(f"{os.getpid()}\n")
        handle.flush()
        try:
            yield path
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def fingerprint(path: Path) -> Any:
    """Returns a JSON-able fingerprint of a file or directory tree."""
    path = Path(path)
    if path.is_dir():
        return {
            p.relative_to(path).as_posix(): fingerprint(p)
            for p in sorted(path.rglob("*"))
            if p.is_file()
        }
    if not path.exists():
        return "missing"
    stat = path.stat()
    if stat.st_size <= SMALL_FILE_BYTES:
        return "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()
    return f"size:{stat.st_size}:mtime_ns:{stat.st_mtime_ns}"


def _fmt(value: float | None, digits: int) -> str:
    """Formats a number, or an empty string when unknown."""
    return "" if value is None else f"{value:.{digits}f}"


def _stored_usage(
    usage: Mapping[str, Any],
) -> tuple[float | None, float | None, float | None, str]:
    """Reads marker usage into Timing fields; missing values become None."""

    def number(name: str) -> float | None:
        value = usage.get(name)
        ok = isinstance(value, (int, float)) and not isinstance(value, bool)
        return float(value) if ok else None

    finished = usage.get("finished")
    return (
        number("wall_s"),
        number("cpu_s"),
        number("peak_rss_mb"),
        finished if isinstance(finished, str) else "",
    )


def _timing_row(timing: Timing) -> str:
    """Formats one ``timings.tsv`` line."""
    return "\t".join(
        [
            timing.step,
            timing.group,
            timing.status,
            str(timing.threads),
            _fmt(timing.wall_s, 2),
            _fmt(timing.cpu_s, 2),
            _fmt(timing.peak_rss_mb, 1),
            timing.finished,
        ]
    )


def summarize_timings(
    timings: Sequence[Timing], total_wall_s: float
) -> dict[str, float]:
    """Sums CPU seconds by origin.

    Args:
        timings: Rows of one invocation.
        total_wall_s: Wall-clock seconds of this invocation.

    Returns:
        ``total_wall_s``, ``cpu_s_executed`` (status ``done``, this run) and
        ``cpu_s_cached`` (status ``cached``, stored by earlier runs; unknown
        values count as 0).
    """
    return {
        "total_wall_s": total_wall_s,
        "cpu_s_executed": sum(
            t.cpu_s or 0.0 for t in timings if t.status == "done"
        ),
        "cpu_s_cached": sum(
            t.cpu_s or 0.0 for t in timings if t.status == "cached"
        ),
    }


def write_timings(
    path: Path,
    timings: Sequence[Timing],
    total_wall_s: float | None = None,
) -> None:
    """Writes ``timings.tsv``: header, one row per task, then summary lines.

    Args:
        path: Output file.
        timings: One row per task.
        total_wall_s: Wall-clock seconds of this invocation; when given,
            ``#``-prefixed, tab-separated summary lines are appended (see
            :func:`summarize_timings`).
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(TIMING_COLUMNS)] + [_timing_row(t) for t in timings]
    if total_wall_s is not None:
        summary = summarize_timings(timings, total_wall_s)
        lines += [f"# {key}\t{value:.2f}" for key, value in summary.items()]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _cpu_seconds() -> float:
    """Returns CPU seconds of this process plus waited-for children."""
    own = resource.getrusage(resource.RUSAGE_SELF)
    kids = resource.getrusage(resource.RUSAGE_CHILDREN)
    return own.ru_utime + own.ru_stime + kids.ru_utime + kids.ru_stime


def _peak_rss_mb() -> float:
    """Returns the larger peak RSS of this process and its children (MiB)."""
    own = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    kids = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    return max(own, kids) / 1024.0  # Linux reports KiB


def _worker(
    conn: multiprocessing.connection.Connection,
    fn: Callable[[Mapping[str, Any]], Any],
    deps: Mapping[str, Any],
) -> None:
    """Runs one task in the worker process and reports through ``conn``."""
    # Own session and process group, so the parent can signal the worker and
    # every tool it starts (Clair3, minimap2, ...) with one ``killpg``.
    with contextlib.suppress(OSError):
        os.setsid()
    # The parent's handlers must not run in the worker.
    for signum in (signal.SIGTERM, signal.SIGINT):
        signal.signal(signum, signal.SIG_DFL)
    start = time.monotonic()
    cpu0 = _cpu_seconds()
    try:
        value = fn(deps)
        message: tuple = ("ok", value)
    except BaseException as err:  # noqa: BLE001 - everything goes to parent
        message = (
            "error",
            f"{type(err).__name__}: {err}",
            traceback.format_exc(),
        )
    usage = (
        time.monotonic() - start,
        _cpu_seconds() - cpu0,
        _peak_rss_mb(),
    )
    try:
        conn.send((message, usage))
    except Exception as err:  # noqa: BLE001 - unpicklable result
        conn.send((("error", f"result not sent: {err}", ""), usage))
    conn.close()


@dataclasses.dataclass
class _Running:
    """Bookkeeping of one running worker."""

    task: Task
    process: multiprocessing.process.BaseProcess
    conn: multiprocessing.connection.Connection
    task_hash: str


@dataclasses.dataclass
class _State:
    """Mutable state of one :meth:`Pool.run` call."""

    pending: dict[str, Task]
    report: Report = dataclasses.field(default_factory=Report)
    hashes: dict[str, str] = dataclasses.field(default_factory=dict)
    queued: dict[str, str] = dataclasses.field(default_factory=dict)
    running: dict[str, _Running] = dataclasses.field(default_factory=dict)
    used_threads: int = 0
    used_mem: float = 0.0
    fatal: tuple[str, str] | None = None


class Pool:
    """Runs task graphs under a thread and memory budget.

    Attributes:
        timings: Timing rows of all :meth:`run` calls, in finishing order.
    """

    def __init__(
        self,
        threads: int,
        markers_dir: Path,
        *,
        max_memory_gb: float = 800.0,
        resume: bool = True,
    ) -> None:
        """Creates a pool.

        Args:
            threads: CPU budget in thread weights (>= 1).
            markers_dir: Directory for done-markers (created).
            max_memory_gb: Memory budget for the task estimates.
            resume: Honor existing markers; markers are written either way.

        Raises:
            SchedulerError: On a budget below 1.
        """
        if threads < 1 or max_memory_gb <= 0:
            raise SchedulerError("threads and memory budget must be positive")
        self.threads = threads
        self.max_memory_gb = max_memory_gb
        self.resume = resume
        self.markers_dir = Path(markers_dir)
        self.markers_dir.mkdir(parents=True, exist_ok=True)
        self.timings: list[Timing] = []
        self._ctx = multiprocessing.get_context("fork")

    def marker_path(self, key: str) -> Path:
        """Returns the marker file of a task key."""
        safe = re.sub(r"[^A-Za-z0-9_.-]", "_", key)[:60]
        digest = hashlib.sha1(key.encode()).hexdigest()[:8]
        return self.markers_dir / f"{safe}.{digest}.json"

    @staticmethod
    def _paths(spec: PathsOrFn, arg: Any) -> list[Path]:
        """Evaluates an inputs or outputs specification."""
        return [Path(p) for p in (spec(arg) if callable(spec) else spec)]

    def _task_hash(
        self, task: Task, deps: Mapping[str, Any], dep_hashes: list[str]
    ) -> str:
        """Hashes name, inputs, parameters and dependency hashes."""
        payload = {
            "version": __version__,
            "step": task.key,
            "params": task.params,
            "inputs": {
                str(p): fingerprint(p) for p in self._paths(task.inputs, deps)
            },
            "deps": dep_hashes,
        }
        text = json.dumps(payload, sort_keys=True, default=str)
        return hashlib.sha256(text.encode()).hexdigest()

    def _try_resume(
        self, task: Task, task_hash: str
    ) -> tuple[bool, Any, Mapping[str, Any]]:
        """Returns ``(True, result, usage)`` if a valid marker exists.

        ``usage`` holds the stored ``wall_s``, ``cpu_s``, ``peak_rss_mb`` and
        ``finished``; it is empty for markers written before these existed.
        """
        if not self.resume:
            return False, None, {}
        path = self.marker_path(task.key)
        try:
            marker = json.loads(path.read_text(encoding="utf-8"))
            if marker["hash"] != task_hash:
                return False, None, {}
            result = task.decode(marker["result"])
            if any(not p.exists() for p in self._paths(task.outputs, result)):
                return False, None, {}
            usage = marker.get("usage")
        except (OSError, ValueError, KeyError, TypeError):
            return False, None, {}
        return True, result, usage if isinstance(usage, dict) else {}

    def _write_marker(
        self,
        task: Task,
        task_hash: str,
        result: Any,
        usage: Mapping[str, Any],
    ) -> None:
        """Writes the done-marker atomically, with the task's real usage."""
        path = self.marker_path(task.key)
        tmp = path.with_suffix(".tmp")
        body = {
            "key": task.key,
            "hash": task_hash,
            "result": task.encode(result),
            "usage": dict(usage),
        }
        tmp.write_text(json.dumps(body, sort_keys=True), encoding="utf-8")
        os.replace(tmp, path)

    def _validate(self, tasks: Sequence[Task]) -> dict[str, Task]:
        """Checks keys, dependencies, weights and cycles."""
        by_key: dict[str, Task] = {}
        for task in tasks:
            if task.key in by_key:
                raise SchedulerError(f"duplicate task {task.key}")
            if not 1 <= task.threads <= self.threads:
                raise SchedulerError(
                    f"{task.key}: threads {task.threads} outside "
                    f"1..{self.threads}"
                )
            if task.memory_gb > self.max_memory_gb:
                raise SchedulerError(
                    f"{task.key}: memory {task.memory_gb} GB above budget"
                )
            by_key[task.key] = task
        for task in tasks:
            for dep in task.deps:
                if dep not in by_key:
                    raise SchedulerError(f"{task.key}: unknown dep {dep}")
        state: dict[str, int] = {}

        def visit(key: str) -> None:
            if state.get(key) == 2:
                return
            if state.get(key) == 1:
                raise SchedulerError(f"dependency cycle at {key}")
            state[key] = 1
            for dep in by_key[key].deps:
                visit(dep)
            state[key] = 2

        for key in by_key:
            visit(key)
        return by_key

    def run(self, tasks: Sequence[Task]) -> Report:
        """Runs a task graph to completion.

        Args:
            tasks: Tasks with unique keys; dependencies stay inside this
                list.

        Returns:
            The report with results, failures and resumed keys.

        Raises:
            SchedulerError: On an invalid graph.
            GlobalStepError: If a task without a group failed. Running
                workers are terminated first.
        """
        state = _State(pending=dict(self._validate(tasks)))
        previous = self._install_signal_handlers()
        try:
            while state.pending or state.running:
                progressed = self._cancel_blocked(state)
                progressed |= self._start_ready(state)
                if state.fatal is not None:
                    break
                if not state.running:
                    if state.pending and not progressed:
                        raise SchedulerError(
                            f"no runnable task among {sorted(state.pending)}"
                        )
                    continue
                waitables: list[Any] = []
                for run_ in state.running.values():
                    waitables += [run_.conn, run_.process.sentinel]
                multiprocessing.connection.wait(waitables, timeout=1.0)
                self._reap(state)
                if state.fatal is not None:
                    break
        finally:
            for run_ in state.running.values():
                self._terminate(run_)
            self._restore_signal_handlers(previous)
        if state.fatal is not None:
            for key in state.pending:
                state.report.failures.setdefault(key, "cancelled: run stopped")
            raise GlobalStepError(*state.fatal, state.report)
        return state.report

    def _fail(
        self, state: _State, task: Task, message: str, wall: float = 0.0
    ) -> None:
        """Records a failed task."""
        state.report.failures[task.key] = message
        self.timings.append(
            Timing(task.key, task.group or "", "failed", task.threads, wall)
        )
        LOGGER.error("%s failed: %s", task.key, message)

    def _cancel_blocked(self, state: _State) -> bool:
        """Cancels pending tasks whose dependency failed or was cancelled."""
        failures = state.report.failures
        progressed = False
        for key in list(state.pending):
            task = state.pending[key]
            bad = [d for d in task.deps if d in failures]
            if bad:
                del state.pending[key]
                failures[key] = f"cancelled: {bad[0]} failed"
                self.timings.append(
                    Timing(key, task.group or "", "cancelled", task.threads)
                )
                progressed = True
        return progressed

    def _start_ready(self, state: _State) -> bool:
        """Resumes or starts ready tasks, highest priority first.

        Returns:
            True if any task changed state. A failing global step sets
            ``state.fatal`` and ends the pass.
        """
        report = state.report
        ready = sorted(
            (
                t
                for t in state.pending.values()
                if all(d in report.results for d in t.deps)
            ),
            key=lambda t: (-t.priority, t.key),
        )
        progressed = False
        for task in ready:
            dep_results = {d: report.results[d] for d in task.deps}
            task_hash = state.queued.get(task.key)
            if task_hash is None:
                # Inputs and dependency results are fixed once a task is
                # ready, so a task waiting for budget is hashed only once.
                try:
                    task_hash = self._task_hash(
                        task, dep_results, [state.hashes[d] for d in task.deps]
                    )
                except Exception as err:  # noqa: BLE001
                    del state.pending[task.key]
                    self._fail(state, task, f"hashing inputs: {err}")
                    progressed = True
                    if task.group is None:
                        state.fatal = (task.key, report.failures[task.key])
                        break
                    continue
                hit, result, stored = self._try_resume(task, task_hash)
                if hit:
                    del state.pending[task.key]
                    state.hashes[task.key] = task_hash
                    report.results[task.key] = result
                    report.skipped.append(task.key)
                    self.timings.append(
                        Timing(
                            task.key,
                            task.group or "",
                            "cached",
                            task.threads,
                            *_stored_usage(stored),
                        )
                    )
                    progressed = True
                    continue
                state.queued[task.key] = task_hash
            if (
                state.used_threads + task.threads > self.threads
                or state.used_mem + task.memory_gb > self.max_memory_gb
            ):
                continue
            self.marker_path(task.key).unlink(missing_ok=True)
            parent, child = self._ctx.Pipe(duplex=False)
            proc = self._ctx.Process(
                target=_worker,
                args=(child, task.fn, dep_results),
                name=task.key,
            )
            proc.start()
            child.close()
            state.running[task.key] = _Running(task, proc, parent, task_hash)
            del state.pending[task.key]
            state.used_threads += task.threads
            state.used_mem += task.memory_gb
            progressed = True
        return progressed

    def _reap(self, state: _State) -> None:
        """Collects finished workers; a global failure sets ``state.fatal``."""
        report = state.report
        for key in list(state.running):
            run_ = state.running[key]
            task = run_.task
            message = None
            done = run_.conn.poll()
            if not done:
                if run_.process.is_alive():
                    continue
                # The worker may have sent its result and exited between the
                # two checks; look again before declaring it dead.
                done = run_.conn.poll()
            if done:
                try:
                    message = run_.conn.recv()
                except EOFError:
                    message = None
            run_.process.join()
            run_.conn.close()
            del state.running[key]
            state.used_threads -= task.threads
            state.used_mem -= task.memory_gb
            if message is None:
                self._fail(
                    state, task, f"worker died (exit {run_.process.exitcode})"
                )
            else:
                outcome, (wall, cpu, rss) = message
                if outcome[0] == "ok":
                    report.results[key] = outcome[1]
                    state.hashes[key] = run_.task_hash
                    finished = datetime.datetime.now(datetime.UTC).strftime(
                        "%Y-%m-%dT%H:%M:%SZ"
                    )
                    usage = {
                        "wall_s": wall,
                        "cpu_s": cpu,
                        "peak_rss_mb": rss,
                        "finished": finished,
                    }
                    self._write_marker(task, run_.task_hash, outcome[1], usage)
                    self.timings.append(
                        Timing(
                            key,
                            task.group or "",
                            "done",
                            task.threads,
                            wall,
                            cpu,
                            rss,
                            finished,
                        )
                    )
                else:
                    self._fail(state, task, outcome[1], wall)
                    self.timings[-1] = dataclasses.replace(
                        self.timings[-1], cpu_s=cpu, peak_rss_mb=rss
                    )
                    LOGGER.debug("%s traceback:\n%s", key, outcome[2])
            if key in report.failures and task.group is None:
                state.fatal = (key, report.failures[key])

    @staticmethod
    def _install_signal_handlers() -> dict[int, Any]:
        """Turns SIGTERM and SIGINT into exceptions inside ``run``.

        The exception unwinds through ``run``'s ``finally`` (workers and
        their process groups are terminated) and through ``work_lock`` (the
        lock is released). Returns the previous handlers; empty when not in
        the main thread, where handlers cannot be set.
        """

        def interrupt(signum: int, _frame: Any) -> None:
            if signum == signal.SIGINT:
                raise KeyboardInterrupt
            raise SystemExit(128 + signum)

        previous: dict[int, Any] = {}
        for signum in (signal.SIGTERM, signal.SIGINT):
            try:
                previous[signum] = signal.signal(signum, interrupt)
            except ValueError:  # not the main thread
                return {}
        return previous

    @staticmethod
    def _restore_signal_handlers(previous: Mapping[int, Any]) -> None:
        """Puts back the handlers saved by ``_install_signal_handlers``."""
        for signum, handler in previous.items():
            signal.signal(signum, handler)

    @staticmethod
    def _signal_group(run_: _Running, signum: int) -> None:
        """Signals the worker's process group (the worker alone as fallback)."""
        pid = run_.process.pid
        try:
            os.killpg(pid, signum)
        except ProcessLookupError:
            # No group yet (worker not past setsid) or all members gone; the
            # pid is only signalled while the worker is not yet reaped.
            if run_.process.exitcode is None:
                with contextlib.suppress(ProcessLookupError):
                    os.kill(pid, signum)

    @classmethod
    def _terminate(cls, run_: _Running) -> None:
        """Stops a worker and every process in its group."""
        cls._signal_group(run_, signal.SIGTERM)
        run_.process.join(timeout=10)
        # Children may outlive the worker; the group is killed regardless.
        cls._signal_group(run_, signal.SIGKILL)
        run_.process.join()
        run_.conn.close()

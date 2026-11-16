# SPDX-License-Identifier: GPL-3.0-or-later
"""Scheduler tests with fake tasks.

Workers are forked, so tasks log to a file instead of changing memory.
"""

import json
import os
import time

import pytest

from paralog2viewer import scheduler
from paralog2viewer.scheduler import (
    GlobalStepError,
    Pool,
    SchedulerError,
    Task,
)


def _logger(log):
    def work(key, sleep=0.0, fail=False, value=None):
        def fn(deps):
            with open(log, "a", encoding="utf-8") as handle:
                handle.write(f"start {key} {time.monotonic():.4f}\n")
            time.sleep(sleep)
            if fail:
                raise RuntimeError("boom")
            with open(log, "a", encoding="utf-8") as handle:
                handle.write(f"end {key} {time.monotonic():.4f}\n")
            return value if value is not None else key

        return fn

    return work


def _events(log):
    if not log.exists():
        return []
    rows = [line.split() for line in log.read_text().splitlines()]
    return [(kind, key, float(t)) for kind, key, t in rows]


def _starts(log):
    return [k for kind, k, _ in _events(log) if kind == "start"]


def _pool(tmp_path, threads=4, **kwargs):
    return Pool(threads, tmp_path / "markers", **kwargs)


def test_dependency_order_and_dep_results(tmp_path):
    log = tmp_path / "log"
    work = _logger(log)
    seen = {}

    def last(deps):
        seen.update(deps)
        (tmp_path / "deps.txt").write_text(",".join(sorted(deps.items())[0]))
        return "c"

    tasks = [
        Task("c", last, deps=["a", "b"]),
        Task("b", work("b", value="vb"), deps=["a"]),
        Task("a", work("a", value="va")),
    ]
    report = _pool(tmp_path).run(tasks)
    assert _starts(log) == ["a", "b"]
    assert report.results == {"a": "va", "b": "vb", "c": "c"}
    assert (tmp_path / "deps.txt").read_text() == "a,va"


def test_thread_budget_never_exceeded_and_used(tmp_path):
    log = tmp_path / "log"
    work = _logger(log)
    tasks = [Task(f"t{i}", work(f"t{i}", 0.3), threads=2) for i in range(6)]
    _pool(tmp_path, threads=4).run(tasks)
    events = sorted(_events(log), key=lambda e: e[2])
    live = peak = 0
    for kind, _, _ in events:
        live += 2 if kind == "start" else -2
        peak = max(peak, live)
    assert peak == 4


def test_memory_guard_serializes_heavy_tasks(tmp_path):
    log = tmp_path / "log"
    work = _logger(log)
    tasks = [Task(f"m{i}", work(f"m{i}", 0.2), memory_gb=3.0) for i in range(3)]
    _pool(tmp_path, threads=4, max_memory_gb=5.0).run(tasks)
    events = sorted(_events(log), key=lambda e: e[2])
    live = peak = 0
    for kind, _, _ in events:
        live += 1 if kind == "start" else -1
        peak = max(peak, live)
    assert peak == 1


def test_task_waiting_for_budget_is_hashed_once(tmp_path, monkeypatch):
    calls = []
    real = Pool._task_hash

    def counting(self, task, deps, dep_hashes):
        calls.append(task.key)
        return real(self, task, deps, dep_hashes)

    monkeypatch.setattr(Pool, "_task_hash", counting)
    work = _logger(tmp_path / "log")
    tasks = [Task(f"w{i}", work(f"w{i}", 0.3), threads=2) for i in range(3)]
    _pool(tmp_path, threads=2).run(tasks)
    assert sorted(calls) == ["w0", "w1", "w2"]


def test_longest_first_priority(tmp_path):
    log = tmp_path / "log"
    work = _logger(log)
    tasks = [
        Task("small", work("small"), priority=1),
        Task("big", work("big"), priority=100),
        Task("mid", work("mid"), priority=10),
    ]
    _pool(tmp_path, threads=1).run(tasks)
    assert _starts(log) == ["big", "mid", "small"]


def test_resume_skips_and_reruns_only_dependents(tmp_path):
    log = tmp_path / "log"
    work = _logger(log)
    data = tmp_path / "in.txt"
    data.write_text("one")

    def tasks():
        return [
            Task("a", work("a"), inputs=[data]),
            Task("b", work("b"), deps=["a"]),
            Task("other", work("other")),
        ]

    _pool(tmp_path).run(tasks())
    assert sorted(_starts(log)) == ["a", "b", "other"]
    log.unlink()
    report = _pool(tmp_path).run(tasks())
    assert _starts(log) == []
    assert sorted(report.skipped) == ["a", "b", "other"]

    data.write_text("two")
    report = _pool(tmp_path).run(tasks())
    assert sorted(_starts(log)) == ["a", "b"]
    assert report.skipped == ["other"]


def test_cached_rows_keep_stored_usage_and_summary(tmp_path):
    log = tmp_path / "log"
    work = _logger(log)

    def tasks():
        return [Task("a", work("a", 0.2)), Task("b", work("b"))]

    first = _pool(tmp_path)
    first.run(tasks())
    real = {t.step: t for t in first.timings}
    assert all(t.status == "done" and t.finished for t in real.values())

    second = _pool(tmp_path)
    second.run(tasks())
    cached = {t.step: t for t in second.timings}
    assert {t.status for t in cached.values()} == {"cached"}
    assert cached["a"].wall_s == real["a"].wall_s >= 0.2
    assert cached["a"].cpu_s == real["a"].cpu_s
    assert cached["a"].peak_rss_mb == real["a"].peak_rss_mb > 0
    assert cached["a"].finished == real["a"].finished

    mixed = _pool(tmp_path)
    mixed.run(tasks() + [Task("c", work("c"))])
    path = tmp_path / "out" / "timings.tsv"
    scheduler.write_timings(path, mixed.timings, 12.5)
    lines = path.read_text().splitlines()
    assert lines[0].split("\t") == list(scheduler.TIMING_COLUMNS)
    summary = {
        k[2:]: float(v) for k, v in (ln.split("\t") for ln in lines[-3:])
    }
    done = {t.step: t for t in mixed.timings}
    assert summary["total_wall_s"] == 12.5
    assert summary["cpu_s_executed"] == pytest.approx(done["c"].cpu_s, abs=0.01)
    assert summary["cpu_s_cached"] == pytest.approx(
        real["a"].cpu_s + real["b"].cpu_s, abs=0.01
    )


def test_old_marker_without_usage_is_cached_with_empty_values(tmp_path):
    work = _logger(tmp_path / "log")
    pool = _pool(tmp_path)
    pool.run([Task("a", work("a"))])
    marker = pool.marker_path("a")
    body = json.loads(marker.read_text())
    del body["usage"]
    marker.write_text(json.dumps(body))

    again = _pool(tmp_path)
    report = again.run([Task("a", work("a"))])
    assert report.skipped == ["a"]
    row = again.timings[0]
    assert row.status == "cached"
    assert (row.wall_s, row.cpu_s, row.peak_rss_mb, row.finished) == (
        None,
        None,
        None,
        "",
    )
    path = tmp_path / "t.tsv"
    scheduler.write_timings(path, again.timings, 1.0)
    cells = path.read_text().splitlines()[1].split("\t")
    assert cells[2:] == ["cached", "1", "", "", "", ""]


def test_no_resume_reruns_everything(tmp_path):
    log = tmp_path / "log"
    work = _logger(log)
    _pool(tmp_path).run([Task("a", work("a"))])
    _pool(tmp_path, resume=False).run([Task("a", work("a"))])
    assert _starts(log) == ["a", "a"]


def test_resume_reruns_when_output_missing(tmp_path):
    log = tmp_path / "log"
    work = _logger(log)
    out = tmp_path / "out.txt"
    out.write_text("x")
    task = Task("a", work("a"), outputs=[out])
    _pool(tmp_path).run([task])
    out.unlink()
    _pool(tmp_path).run([task])
    assert _starts(log) == ["a", "a"]


def test_group_failure_is_isolated(tmp_path):
    log = tmp_path / "log"
    work = _logger(log)
    tasks = [
        Task("p:h1", work("p:h1", fail=True), group="h1"),
        Task("c:h1", work("c:h1"), deps=["p:h1"], group="h1"),
        Task("p:h2", work("p:h2"), group="h2"),
        Task("c:h2", work("c:h2"), deps=["p:h2"], group="h2"),
    ]
    report = _pool(tmp_path).run(tasks)
    assert "boom" in report.failures["p:h1"]
    assert report.failures["c:h1"].startswith("cancelled")
    assert report.results["c:h2"] == "c:h2"
    assert "c:h1" not in _starts(log)


def test_global_failure_stops_run(tmp_path):
    log = tmp_path / "log"
    work = _logger(log)
    tasks = [
        Task("g", work("g", fail=True)),
        Task("after", work("after"), deps=["g"]),
    ]
    with pytest.raises(GlobalStepError) as exc:
        _pool(tmp_path).run(tasks)
    assert exc.value.key == "g"
    assert "after" not in _starts(log)
    assert "after" in exc.value.report.failures


def test_worker_death_is_a_failure(tmp_path):
    def die(deps):
        os._exit(3)

    report = _pool(tmp_path).run([Task("d", die, group="x")])
    assert "worker died" in report.failures["d"]


def test_timings_file(tmp_path):
    log = tmp_path / "log"
    work = _logger(log)
    pool = _pool(tmp_path)
    pool.run(
        [Task("a", work("a", 0.2)), Task("f", work("f", fail=True), group="g")]
    )
    path = tmp_path / "out" / "timings.tsv"
    scheduler.write_timings(path, pool.timings)
    rows = [line.split("\t") for line in path.read_text().splitlines()]
    assert rows[0] == list(scheduler.TIMING_COLUMNS)
    by_step = {r[0]: r for r in rows[1:]}
    assert by_step["a"][2] == "done" and float(by_step["a"][4]) >= 0.2
    assert float(by_step["a"][6]) > 0
    assert by_step["f"][2] == "failed"


def test_invalid_graphs(tmp_path):
    fn = _logger(tmp_path / "log")
    pool = _pool(tmp_path)
    with pytest.raises(SchedulerError, match="unknown dep"):
        pool.run([Task("a", fn("a"), deps=["zz"])])
    with pytest.raises(SchedulerError, match="cycle"):
        pool.run(
            [Task("a", fn("a"), deps=["b"]), Task("b", fn("b"), deps=["a"])]
        )
    with pytest.raises(SchedulerError, match="threads"):
        pool.run([Task("a", fn("a"), threads=9)])
    with pytest.raises(SchedulerError, match="duplicate"):
        pool.run([Task("a", fn("a")), Task("a", fn("a"))])


def test_work_lock_blocks_second_holder(tmp_path):
    with scheduler.work_lock(tmp_path / "w"):
        with pytest.raises(scheduler.LockedError):
            with scheduler.work_lock(tmp_path / "w"):
                pass
    with scheduler.work_lock(tmp_path / "w"):
        pass


def test_result_sent_just_before_exit_is_not_a_death(tmp_path, monkeypatch):
    """poll() False, worker exits after sending: the result must be kept."""
    pool = _pool(tmp_path)
    real_poll = scheduler.multiprocessing.connection.Connection.poll
    calls = {"n": 0}

    def late_poll(self, *args, **kwargs):
        # The first poll misses the message; later polls see it.
        calls["n"] += 1
        if calls["n"] == 1:
            time.sleep(0.5)  # let the worker send and exit
            return False
        return real_poll(self, *args, **kwargs)

    monkeypatch.setattr(
        scheduler.multiprocessing.connection.Connection, "poll", late_poll
    )
    report = pool.run([Task("quick", lambda deps: "done", group="x")])
    assert report.failures == {}
    assert report.results["quick"] == "done"


def test_children_die_with_the_worker_group(tmp_path):
    """terminate() kills grandchildren started by a task."""
    pid_file = tmp_path / "child.pid"

    def spawn(deps):
        import subprocess

        child = subprocess.Popen(["sleep", "60"])
        pid_file.write_text(str(child.pid))
        time.sleep(60)

    pool = _pool(tmp_path)
    ctx = pool._ctx
    parent, child_conn = ctx.Pipe(duplex=False)
    proc = ctx.Process(target=scheduler._worker, args=(child_conn, spawn, {}))
    proc.start()
    child_conn.close()
    run_ = scheduler._Running(Task("s", spawn), proc, parent, "h")
    for _ in range(100):
        if pid_file.exists() and pid_file.read_text():
            break
        time.sleep(0.05)
    grandchild = int(pid_file.read_text())
    os.kill(grandchild, 0)  # alive
    pool._terminate(run_)
    time.sleep(0.2)
    with pytest.raises(ProcessLookupError):
        os.kill(grandchild, 0)


def test_sigterm_in_run_stops_workers_and_releases_lock(tmp_path):
    """SIGTERM while a task runs unwinds run(), so the lock is released."""
    import signal
    import threading

    from paralog2viewer.scheduler import work_lock

    work = tmp_path / "w"
    timer = threading.Timer(0.5, os.kill, args=(os.getpid(), signal.SIGTERM))
    timer.start()
    try:
        with pytest.raises(SystemExit), work_lock(work):
            _pool(tmp_path).run([Task("slow", lambda d: time.sleep(30))])
    finally:
        timer.cancel()
    with work_lock(work):  # lock is free again
        pass
    assert signal.getsignal(signal.SIGTERM) is not None


def test_version_change_invalidates_resume(tmp_path, monkeypatch):
    log = tmp_path / "log"
    work = _logger(log)
    _pool(tmp_path).run([Task("a", work("a"))])
    monkeypatch.setattr(scheduler, "__version__", "99.0.0")
    report = _pool(tmp_path).run([Task("a", work("a"))])
    assert report.skipped == []

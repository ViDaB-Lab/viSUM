"""Bounded scheduling, spawn safety, and deterministic candidate evaluation."""

from argparse import Namespace
from concurrent.futures.process import BrokenProcessPool
import copy
import multiprocessing
import os
from pathlib import Path
import sys
import time

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
import run_viharmony as harmony


def context():
    return harmony.CandidateContext("sample", "dna", 1000, False, False, False, True)


def catalog():
    return harmony.make_msl_catalog([("Riboviria", "", "", "", "", "", "", "")])


def task(index):
    final_id = f"c{index}"
    return harmony.CandidateTask(
        final_id, 100, dict(record_type="input_contig", coordinates="", original_length="100"),
        final_id, dict(original_id="same-name", original_header="same-name comment"), "viral",
        [dict(sequence_id=final_id, tool="genomad", classification="virus", evidence_strength="strong", r__Realm="Riboviria")],
        [], (),
    )


class FakeFuture:
    def __init__(self, values, fail=False):
        self.values, self.fail, self.cancelled = values, fail, False

    def result(self):
        if self.fail:
            raise ValueError("sentinel worker failure")
        return self.values

    def cancel(self):
        self.cancelled = True


class FakeExecutor:
    def __init__(self, fail=False):
        self.futures, self.fail = [], fail

    def submit(self, function, batch):
        future = FakeFuture(batch, self.fail and not self.futures)
        self.futures.append(future)
        return future


def test_queue_is_lazy_bounded_and_closes_pending_work():
    consumed = []
    def inputs():
        for index in range(10000):
            consumed.append(index)
            yield index
    executor = FakeExecutor()
    results = harmony._ordered_candidate_batches(executor, inputs(), 3, 4)
    assert next(results) == 0
    assert consumed == list(range(12))
    assert next(results) == 1
    assert consumed == list(range(12))
    results.close()
    assert all(f.cancelled for f in executor.futures[1:])
    assert len(consumed) == 12  # Did not enqueue the rest of the dataset.


def test_batched_results_remain_ordered_at_boundaries():
    assert list(harmony._ordered_candidate_batches(FakeExecutor(), range(101), 3, 4)) == list(range(101))


def test_failed_future_cancels_pending_work_and_propagates():
    executor = FakeExecutor(fail=True)
    with pytest.raises(ValueError, match="sentinel worker failure"):
        list(harmony._ordered_candidate_batches(executor, range(100), 1, 4))
    assert len(executor.futures) == 4
    assert all(f.cancelled for f in executor.futures[1:])


@pytest.mark.parametrize("workers", [0, -1, True, 1.5, "2", None])
def test_invalid_worker_counts_fail_before_loading_inputs(workers):
    with pytest.raises(ValueError, match="positive integer"):
        harmony.run(Namespace(workers=workers))


def test_serial_mode_never_creates_a_pool(monkeypatch):
    def forbidden(**kwargs):
        raise AssertionError("Serial path created a process pool")
    monkeypatch.setattr(harmony, "ProcessPoolExecutor", forbidden)
    with harmony.candidate_results([task(0)], context(), catalog(), 1) as results:
        assert [r.metadata["final_sequence_id"] for r in results] == ["c0"]


def test_candidate_evaluation_does_not_mutate_its_evidence():
    candidate = task(0)
    original = copy.deepcopy(candidate)
    harmony.evaluate_candidate(candidate, context(), catalog())
    assert candidate == original


@pytest.mark.parametrize("workers", [2, 4, 12])
def test_spawned_candidate_results_match_serial_and_workers_exit(workers):
    original_children = {child.pid for child in multiprocessing.active_children()}
    tasks = [task(i) for i in range(48)]
    expected = [harmony.evaluate_candidate(item, context(), catalog()) for item in tasks]
    with harmony.candidate_results(iter(tasks), context(), catalog(), workers, 2) as results:
        assert list(results) == expected
    assert {child.pid for child in multiprocessing.active_children()} <= original_children


def _delayed_echo(batch):
    time.sleep(0.2 if batch[0] == 0 else 0.01)
    return [(value, os.getpid()) for value in batch]


def _crash(batch):
    os._exit(7)


def test_out_of_order_completion_still_returns_input_order(monkeypatch):
    monkeypatch.setattr(harmony, "_evaluate_candidate_batch", _delayed_echo)
    with harmony.candidate_results(range(40), context(), catalog(), 4, 1) as results:
        values = list(results)
    assert [value for value, _ in values] == list(range(40))
    assert 1 < len({pid for _, pid in values}) <= 4


def test_worker_exception_reaches_parent_and_pool_is_reaped():
    original_children = {child.pid for child in multiprocessing.active_children()}
    invalid = harmony.CandidateTask("bad", 100, {}, "bad", {}, "viral", [], [], ())
    with pytest.raises(KeyError, match="record_type"):
        with harmony.candidate_results([invalid, task(1)], context(), catalog(), 2, 1) as results:
            list(results)
    assert {child.pid for child in multiprocessing.active_children()} <= original_children


def test_abrupt_worker_exit_does_not_hang_or_succeed(monkeypatch):
    original_children = {child.pid for child in multiprocessing.active_children()}
    monkeypatch.setattr(harmony, "_evaluate_candidate_batch", _crash)
    with pytest.raises(BrokenProcessPool):
        with harmony.candidate_results(range(4), context(), catalog(), 2, 1) as results:
            list(results)
    assert {child.pid for child in multiprocessing.active_children()} <= original_children

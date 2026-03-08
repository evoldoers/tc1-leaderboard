"""Integration tests: run full validation on test-entry."""

import json
import pytest
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from validate import validate_entry, update_leaderboard


TEST_ENTRY = Path(__file__).parent.parent.parent / "entries" / "test-entry"


def test_validate_test_entry_passes():
    report = validate_entry(TEST_ENTRY, commit="abc1234")
    assert report["team"] == "test-entry"
    assert report["n_sequences"] == 3
    assert report["n_families"] == 1
    assert report["scores"]["total"] > 0

    # All checks should pass
    for cat, check in report["checks"].items():
        assert check["status"] != "fail", f"{cat} failed: {check['messages']}"


def test_validate_test_entry_scores_reasonable():
    report = validate_entry(TEST_ENTRY, commit="abc1234")
    scores = report["scores"]
    assert 0 < scores["total"] <= 100
    assert 0 <= scores["diversity"] <= 1
    assert 0 <= scores["annotation"] <= 1
    assert 0 <= scores["functionality"] <= 1
    assert 0 <= scores["size"] <= 1


def test_validate_missing_dir(tmp_path):
    report = validate_entry(tmp_path / "nonexistent")
    assert report["scores"]["total"] == 0.0


def test_validate_empty_dir(tmp_path):
    entry = tmp_path / "empty-team"
    entry.mkdir()
    report = validate_entry(entry)
    assert report["scores"]["total"] == 0.0
    assert report["checks"]["format"]["status"] == "fail"


def test_update_leaderboard(tmp_path):
    scores_dir = tmp_path / "scores"
    scores_dir.mkdir()

    # Write two score files
    for team, total in [("alpha", 80.0), ("beta", 60.0)]:
        (scores_dir / f"{team}.json").write_text(json.dumps({
            "team": team,
            "scores": {"total": total, "diversity": 0.5, "annotation": 0.9,
                       "functionality": 0.3, "size": 0.8},
            "n_sequences": 10,
            "n_families": 2,
            "commit": "abc",
        }))

    lb_path = tmp_path / "leaderboard.json"
    update_leaderboard(scores_dir, lb_path)

    data = json.loads(lb_path.read_text())
    assert len(data["entries"]) == 2
    assert data["entries"][0]["team"] == "alpha"  # higher score first
    assert data["entries"][1]["team"] == "beta"


def test_deleted_entry_removed_from_leaderboard(tmp_path):
    """When a team's score file is deleted, leaderboard should no longer include them."""
    scores_dir = tmp_path / "scores"
    scores_dir.mkdir()

    for team, total in [("alive", 80.0), ("doomed", 60.0)]:
        (scores_dir / f"{team}.json").write_text(json.dumps({
            "team": team,
            "scores": {"total": total, "diversity": 0.5, "annotation": 0.9,
                       "functionality": 0.3, "size": 0.8},
            "n_sequences": 10,
            "n_families": 2,
            "commit": "abc",
        }))

    lb_path = tmp_path / "leaderboard.json"
    update_leaderboard(scores_dir, lb_path)
    data = json.loads(lb_path.read_text())
    assert len(data["entries"]) == 2

    # Delete doomed team's score
    (scores_dir / "doomed.json").unlink()
    update_leaderboard(scores_dir, lb_path)

    data = json.loads(lb_path.read_text())
    assert len(data["entries"]) == 1
    assert data["entries"][0]["team"] == "alive"

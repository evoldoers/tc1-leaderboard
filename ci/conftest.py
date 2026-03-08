"""Pytest configuration and shared fixtures for CI tests."""

from pathlib import Path

import pytest

FIXTURES_DIR = Path(__file__).parent / "tests" / "fixtures"
TEST_ENTRY_DIR = Path(__file__).parent.parent / "entries" / "test-entry"


@pytest.fixture
def test_entry_dir():
    """Path to the test-entry used for integration tests."""
    return TEST_ENTRY_DIR


@pytest.fixture
def fixtures_dir():
    """Path to the test fixtures directory."""
    return FIXTURES_DIR

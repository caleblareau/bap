"""Shared pytest fixtures and helpers for the bap test suite."""

import os

import pytest

# Absolute path to the tests/ directory and the bundled data fixtures.
TESTS_DIR = os.path.dirname(os.path.realpath(__file__))
DATA_DIR = os.path.join(TESTS_DIR, "data")
FASTQ_DIR = os.path.join(TESTS_DIR, "fastq_br")


@pytest.fixture
def data_dir():
    return DATA_DIR


@pytest.fixture
def small_bam():
    """A small, indexed hg19 bam bundled with the repo."""
    return os.path.join(DATA_DIR, "jaccardPairsForIGV.bam")

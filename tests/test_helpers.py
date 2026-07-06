"""Unit tests for the pure helper functions in bap.bapHelp.

These require no external tools and run in well under a second. They also lock
in fixes for functions that were previously broken under Python 3 (e.g.
``string_hamming_distance`` used the Python-2-only ``itertools.imap``).
"""

import os
import sys

import pytest

from bap import bapHelp
from bap.bapHelp import (
    bap_version,
    run_cmd,
    string_hamming_distance,
    intersection,
    findIdx,
    file_len,
    verify_file,
    get_software_path,
    rev_comp,
    available_cpu_count,
    gettime,
)
from bap.bap2ProjectClass import mitoChr


# ---------------------------------------------------------------------------
# Sequence helpers
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("seq,expected", [
    ("ACGT", "ACGT"),   # reverse complement is a palindrome here
    ("AAAA", "TTTT"),
    ("ATCG", "CGAT"),
    ("N", "N"),
    ("", ""),
])
def test_rev_comp(seq, expected):
    assert rev_comp(seq) == expected


@pytest.mark.parametrize("a,b,expected", [
    ("karolin", "kathrin", 3),
    ("AAAA", "AAAA", 0),
    ("ACGT", "TGCA", 4),
    ("", "", 0),
])
def test_string_hamming_distance(a, b, expected):
    assert string_hamming_distance(a, b) == expected


# ---------------------------------------------------------------------------
# List helpers
# ---------------------------------------------------------------------------

def test_intersection():
    assert intersection([1, 2, 3], [2, 3, 4]) == [2, 3]
    assert intersection(["chr1", "chr2"], ["chrX"]) == []


def test_findIdx():
    assert findIdx(["a", "b", "c"], ["b"]) == [1]
    assert findIdx(["a", "b", "c"], ["a", "c"]) == [0, 2]


# ---------------------------------------------------------------------------
# File helpers
# ---------------------------------------------------------------------------

def test_file_len(tmp_path):
    p = tmp_path / "lines.txt"
    p.write_text("a\nb\nc\n")
    assert file_len(str(p)) == 3


def test_verify_file_ok(tmp_path):
    p = tmp_path / "ok.txt"
    p.write_text("hello\n")
    assert verify_file(str(p)) == str(p)


def test_verify_file_missing_exits(tmp_path):
    with pytest.raises(SystemExit):
        verify_file(str(tmp_path / "does_not_exist.txt"))


# ---------------------------------------------------------------------------
# Tool discovery
# ---------------------------------------------------------------------------

def test_get_software_path_found():
    # The Python interpreter is always discoverable via its basename on PATH
    # is not guaranteed, so use an absolute-path override instead.
    assert get_software_path("definitely-not-a-real-tool", sys.executable) == sys.executable


def test_get_software_path_missing_exits():
    with pytest.raises(SystemExit):
        get_software_path("definitely-not-a-real-tool-xyz", "")


# ---------------------------------------------------------------------------
# Reference-genome mitochondrial chromosome naming
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("ref,expected", [
    ("hg19", "chrM"),
    ("mm10", "chrM"),
    ("GRCh38", "MT"),
    ("hg19_mm10_c", "humanM"),
    ("some_unknown_ref", "hg19_chrM"),
])
def test_mitoChr_defaults(ref, expected):
    assert mitoChr(ref, "default") == expected


def test_mitoChr_explicit_override():
    assert mitoChr("hg19", "myMito") == "myMito"


# ---------------------------------------------------------------------------
# Misc helpers
# ---------------------------------------------------------------------------

def test_available_cpu_count():
    n = available_cpu_count()
    assert isinstance(n, int)
    assert n >= 1


def test_bap_version_is_string():
    v = bap_version()
    assert isinstance(v, str)
    assert v != ""


def test_gettime_format():
    # Ends with the ": " separator used throughout the log messages
    assert gettime().endswith(": ")


# ---------------------------------------------------------------------------
# run_cmd: the subprocess wrapper that replaced os.system(<string>)
# ---------------------------------------------------------------------------

def test_run_cmd_success():
    result = run_cmd([sys.executable, "-c", "pass"])
    assert result.returncode == 0


def test_run_cmd_failure_exits():
    with pytest.raises(SystemExit):
        run_cmd([sys.executable, "-c", "import sys; sys.exit(3)"])


def test_run_cmd_failure_no_check():
    result = run_cmd([sys.executable, "-c", "import sys; sys.exit(3)"], check=False)
    assert result.returncode == 3


def test_run_cmd_captures_combined_log(tmp_path):
    log = tmp_path / "out.log"
    run_cmd(
        [sys.executable, "-c", "import sys; print('to_out'); print('to_err', file=sys.stderr)"],
        log_file=str(log),
    )
    contents = log.read_text()
    assert "to_out" in contents
    assert "to_err" in contents


def test_run_cmd_stderr_only_log(tmp_path):
    log = tmp_path / "err.log"
    run_cmd(
        [sys.executable, "-c", "import sys; print('to_out'); print('to_err', file=sys.stderr)"],
        log_file=str(log),
        log_stderr_only=True,
    )
    contents = log.read_text()
    assert "to_err" in contents
    assert "to_out" not in contents

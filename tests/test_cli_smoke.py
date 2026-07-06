"""Fast, hermetic smoke tests for the click CLIs.

These invoke the command objects in-process with click's ``CliRunner`` and do
not require samtools/bedtools/R/snakemake. They verify that every active entry
point imports cleanly, renders ``--help``, and reports a version.
"""

import pytest
from click.testing import CliRunner

from bap.cli_bap2 import main as bap2_main
from bap.cli_bap_frag import main as bap_frag_main
from bap.cli_bap_bulk_frag import main as bap_bulk_frag_main
from bap.cli_reanno import main as bap_reanno_main
from bap.barcode.cli_barcode import main as bap_barcode_main

ALL_COMMANDS = [
    bap2_main,
    bap_frag_main,
    bap_bulk_frag_main,
    bap_reanno_main,
    bap_barcode_main,
]


@pytest.mark.parametrize("command", ALL_COMMANDS)
def test_help_renders(command):
    result = CliRunner().invoke(command, ["--help"])
    assert result.exit_code == 0
    assert "Usage:" in result.output


@pytest.mark.parametrize("command", ALL_COMMANDS)
def test_version_reports(command):
    result = CliRunner().invoke(command, ["--version"])
    assert result.exit_code == 0
    assert "version" in result.output.lower()


def test_bap2_support_lists_builtin_genomes():
    # `support` prints the genome list and then sys.exit()s with a message,
    # so the exit code is non-zero but the output holds the genome names.
    result = CliRunner().invoke(bap2_main, ["support"])
    assert "hg19" in result.output
    assert "GRCh38" in result.output

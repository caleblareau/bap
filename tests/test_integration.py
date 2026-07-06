"""End-to-end integration tests for the bap2 pipeline.

These exercise the real workflow and therefore require the external tools the
pipeline orchestrates (samtools, bedtools, R + pipeline R packages, snakemake).
They are automatically skipped when those tools are not on PATH, so the default
``pytest`` run stays green on a bare machine. Run them explicitly with:

    pytest -m integration
"""

import shutil

import pytest
from click.testing import CliRunner

from bap.cli_bap2 import main as bap2_main

# External tools the full pipeline shells out to; skip when any are missing.
REQUIRED_TOOLS = ("samtools", "bedtools", "R", "snakemake")
_missing = [t for t in REQUIRED_TOOLS if shutil.which(t) is None]
requires_external_tools = pytest.mark.skipif(
    bool(_missing), reason="external tools not on PATH: " + ", ".join(_missing)
)


@pytest.mark.integration
@requires_external_tools
def test_bap2_bam_end_to_end(small_bam, tmp_path):
    outdir = tmp_path / "bap2_out"
    result = CliRunner().invoke(
        bap2_main,
        ["bam", "-i", small_bam, "-bt", "XB", "-r", "hg19", "-z", "-o", str(outdir)],
    )

    # Surface pipeline output if something went wrong, to aid debugging.
    assert result.exit_code == 0, result.output

    name = "jaccardPairsForIGV"  # derived from the input .bam prefix
    final_dir = outdir / "final"
    assert (final_dir / (name + ".bap.bam")).exists()
    assert (final_dir / (name + ".barcodeTranslate.tsv")).exists()
    assert (final_dir / (name + ".fragments.tsv.gz")).exists()

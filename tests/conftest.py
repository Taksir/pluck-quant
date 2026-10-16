import shutil
from pathlib import Path
import pytest

DATA = Path(__file__).parent / "data"


@pytest.fixture(scope="session")
def panel_dir(tmp_path_factory):
    from pluck_quant import panel
    out = tmp_path_factory.mktemp("panel")
    panel.compile_panel(DATA / "mini_design" / "inputs", DATA / "mini_design" / "inputs" / "probes_s26.npz",
                        DATA / "genes.json", 50, out, "reverse")
    return out


@pytest.fixture(scope="session")
def counts_file(panel_dir, tmp_path_factory):
    from pluck_quant import count
    out = tmp_path_factory.mktemp("counts") / "cc021.counts.txt"
    count.count_stream(panel_dir / "panel.txt", DATA / "reads" / "CC021_first500k.fastq.gz", out, threads=4)
    return out

"""
Unit tests for bin/make_samplesheet.py
"""
import gzip, csv, sys, tempfile
from pathlib import Path

# Add bin/ to path so we can import the script
sys.path.insert(0, str(Path(__file__).parent.parent / 'bin'))
from make_samplesheet import detect_marker, parse_fastq_dir, load_metadata, MARKER_ALIASES


# ── detect_marker ─────────────────────────────────────────────────────────────

def test_detect_marker_dash_separator():
    sid, marker = detect_marker("H4D-EA100_2024-18S")
    assert sid == "H4D-EA100_2024"
    assert marker == "18S"

def test_detect_marker_underscore_separator():
    sid, marker = detect_marker("SAMPLE_16S")
    assert sid == "SAMPLE"
    assert marker == "16S"

def test_detect_marker_coi_alias():
    sid, marker = detect_marker("POND_A-COI")
    assert marker == "CO1"

def test_detect_marker_its1_alias():
    sid, marker = detect_marker("RIVER_B-ITS1")
    assert marker == "ITS"

def test_detect_marker_rbcl_lowercase():
    sid, marker = detect_marker("SAMPLE-rbcL")
    assert marker == "RBCL"

def test_detect_marker_unknown_returns_none():
    sid, marker = detect_marker("SAMPLE_UNKNOWN")
    assert sid is None
    assert marker is None

def test_detect_marker_no_separator():
    sid, marker = detect_marker("NOSEPARATOR")
    assert sid is None
    assert marker is None


# ── parse_fastq_dir ────────────────────────────────────────────────────────────

def _make_fastq(path: Path):
    """Create a minimal valid gzipped FASTQ file."""
    with gzip.open(path, 'wt') as f:
        f.write("@READ1\nACGT\n+\nIIII\n")


def test_parse_fastq_dir_paired():
    with tempfile.TemporaryDirectory() as tmpdir:
        d = Path(tmpdir)
        _make_fastq(d / "SAMPLE_A-16S_S1_R1_001.fastq.gz")
        _make_fastq(d / "SAMPLE_A-16S_S1_R2_001.fastq.gz")
        samples, skipped = parse_fastq_dir(d, absolute=False)
        assert len(samples) == 1
        assert ('SAMPLE_A', '16S') in samples
        assert 'R1' in samples[('SAMPLE_A', '16S')]
        assert 'R2' in samples[('SAMPLE_A', '16S')]
        assert skipped == []


def test_parse_fastq_dir_multi_marker():
    with tempfile.TemporaryDirectory() as tmpdir:
        d = Path(tmpdir)
        for marker in ['16S', '18S', 'ITS']:
            _make_fastq(d / f"SAMP-{marker}_S1_R1_001.fastq.gz")
            _make_fastq(d / f"SAMP-{marker}_S1_R2_001.fastq.gz")
        samples, _ = parse_fastq_dir(d, absolute=False)
        assert len(samples) == 3
        markers_found = {m for _, m in samples}
        assert markers_found == {'16S', '18S', 'ITS'}


def test_parse_fastq_dir_skips_unrecognised():
    with tempfile.TemporaryDirectory() as tmpdir:
        d = Path(tmpdir)
        _make_fastq(d / "SAMPLE_A-16S_S1_R1_001.fastq.gz")
        _make_fastq(d / "SAMPLE_A-16S_S1_R2_001.fastq.gz")
        (d / "random_file.txt").write_text("not a fastq")
        (d / "bad_name_R1_001.fastq.gz").write_bytes(b"")
        samples, skipped = parse_fastq_dir(d, absolute=False)
        assert len(samples) == 1
        assert len(skipped) == 1  # only bad_name; .txt is not .fastq.gz


def test_parse_fastq_dir_single_end():
    with tempfile.TemporaryDirectory() as tmpdir:
        d = Path(tmpdir)
        _make_fastq(d / "SAMPLE_A-CO1_S1_R1_001.fastq.gz")
        samples, _ = parse_fastq_dir(d, absolute=False)
        key = ('SAMPLE_A', 'CO1')
        assert key in samples
        assert 'R1' in samples[key]
        assert 'R2' not in samples[key]


# ── load_metadata ──────────────────────────────────────────────────────────────

def test_load_metadata_tsv():
    with tempfile.NamedTemporaryFile(suffix='.tsv', mode='w', delete=False) as f:
        f.write("sample_id\thabitat\tsite\n")
        f.write("SAMPLE_A\tpond\tNorth\n")
        f.write("SAMPLE_B\triver\tSouth\n")
        name = f.name
    id_col, header, rows = load_metadata(Path(name))
    assert id_col == 'sample_id'
    assert 'SAMPLE_A' in rows
    assert rows['SAMPLE_A']['habitat'] == 'pond'


def test_load_metadata_csv():
    with tempfile.NamedTemporaryFile(suffix='.csv', mode='w', delete=False) as f:
        f.write("sample,group\n")
        f.write("RIVER_C,upland\n")
        name = f.name
    id_col, header, rows = load_metadata(Path(name))
    assert 'RIVER_C' in rows
    assert rows['RIVER_C']['group'] == 'upland'


def test_marker_aliases_canonical():
    """All alias keys should map to canonical marker names."""
    canonical = {'16S', '18S', 'ITS', 'CO1', '12S', 'RBCL'}
    for alias, canon in MARKER_ALIASES.items():
        assert canon in canonical, f"{alias} maps to {canon} which is not canonical"

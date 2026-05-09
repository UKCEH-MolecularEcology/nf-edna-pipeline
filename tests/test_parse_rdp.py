"""
Unit tests for bin/parse_rdp.py
"""
import sys, subprocess, tempfile
from pathlib import Path

BIN = Path(__file__).parent.parent / 'bin' / 'parse_rdp.py'
sys.path.insert(0, str(Path(__file__).parent.parent / 'bin'))
from parse_rdp import parse_fixrank_line, RANK_LOOKUP, OUT_RANKS


# ── parse_fixrank_line ────────────────────────────────────────────────────────

def _rdp_cols(seq_id, ranks):
    """Build a cols list as RDP fixrank produces: [seqid, orientation, rank, name, conf, ...]"""
    cols = [seq_id, '+']
    for rank, name, conf in ranks:
        cols += [rank, name, str(conf)]
    return cols


def test_basic_ranks_parsed():
    cols = _rdp_cols('ASV1', [
        ('superkingdom', 'Plantae',         0.95),
        ('phylum',       'Tracheophyta',    0.92),
        ('class',        'Magnoliopsida',   0.90),
        ('order',        'Lamiales',        0.88),
        ('family',       'Plantaginaceae',  0.85),
        ('genus',        'Plantago',        0.80),
        ('species',      'Plantago lanceolata', 0.75),
    ])
    rank_val, rank_boot = parse_fixrank_line(cols)
    assert rank_val['genus'] == 'Plantago'
    assert rank_boot['genus'] == 0.80
    assert rank_val['phylum'] == 'Tracheophyta'


def test_confidence_values_stored():
    cols = _rdp_cols('ASV2', [
        ('genus',   'Rosa',  0.63),
        ('species', 'Rosa canina', 0.55),
    ])
    _, rank_boot = parse_fixrank_line(cols)
    assert abs(rank_boot['genus'] - 0.63) < 1e-9
    assert abs(rank_boot['species'] - 0.55) < 1e-9


def test_empty_cols_returns_empty_dicts():
    rank_val, rank_boot = parse_fixrank_line(['ASV1', '+'])
    assert rank_val == {}
    assert rank_boot == {}


# ── RANK_LOOKUP coverage ──────────────────────────────────────────────────────

def test_rank_lookup_covers_all_output_ranks():
    for rank in OUT_RANKS:
        assert rank in RANK_LOOKUP, f"Missing output rank in RANK_LOOKUP: {rank}"


def test_rank_lookup_superkingdom_or_kingdom():
    """Kingdom column should accept both 'superkingdom' and 'kingdom'."""
    assert 'superkingdom' in RANK_LOOKUP['Kingdom']
    assert 'kingdom' in RANK_LOOKUP['Kingdom']


# ── End-to-end via subprocess ─────────────────────────────────────────────────

def _run(in_text, cutoff='0.70'):
    with tempfile.TemporaryDirectory() as tmpdir:
        inp = Path(tmpdir) / 'rdp.txt'
        out = Path(tmpdir) / 'taxonomy.tsv'
        inp.write_text(in_text)
        result = subprocess.run(
            [sys.executable, str(BIN), str(inp), str(out), cutoff],
            capture_output=True, text=True
        )
        assert result.returncode == 0, result.stderr
        return out.read_text().splitlines()


def test_e2e_basic():
    lines = _run(
        'ASV1\t+\t'
        'superkingdom\tPlantae\t0.95\t'
        'phylum\tTracheophyta\t0.92\t'
        'class\tMagnoliopsida\t0.90\t'
        'order\tLamiales\t0.88\t'
        'family\tPlantaginaceae\t0.85\t'
        'genus\tPlantago\t0.80\t'
        'species\tPlantago lanceolata\t0.75\n'
    )
    header = lines[0].split('\t')
    data   = lines[1].split('\t')
    assert data[0] == 'ASV1'
    assert data[header.index('Genus')] == 'Plantago'
    assert data[header.index('Species')] == 'Plantago lanceolata'


def test_e2e_below_cutoff_blanked():
    lines = _run(
        'ASV1\t+\t'
        'superkingdom\tPlantae\t0.95\t'
        'phylum\tTracheophyta\t0.92\t'
        'class\tMagnoliopsida\t0.90\t'
        'order\tLamiales\t0.88\t'
        'family\tPlantaginaceae\t0.50\t'    # below cutoff
        'genus\tPlantago\t0.40\t'           # below cutoff
        'species\tPlantago lanceolata\t0.30\n'  # below cutoff
    )
    header = lines[0].split('\t')
    data   = lines[1].split('\t')
    assert data[header.index('Family')] == ''
    assert data[header.index('Genus')]  == ''
    assert data[header.index('Order')]  == 'Lamiales'


def test_e2e_empty_input_writes_header_only():
    lines = _run('')
    assert len(lines) == 1
    assert lines[0].startswith('asv_id\t')

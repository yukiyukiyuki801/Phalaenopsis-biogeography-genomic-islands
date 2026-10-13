#!/usr/bin/env python3
"""Refresh unconstrained tree inputs without changing the reported geography."""
from pathlib import Path
import argparse
import csv
import io
import re
from Bio import Phylo

ROOT = Path(__file__).resolve().parents[2]

def prepare(calibration, output_root=None):
    base = ROOT / 'orchid_paper2'
    summary = base / '16_rad_beast/final_unconstrained_60m_summary'
    inp = (output_root or base / '19_biogeobears/unconstrained_29tip') / calibration / 'max_range_3/inputs'
    # Refuse empty directories rather than silently changing the fitted model.
    for name in ('geography.data', 'area_order.tsv', 'allowed_states.tsv', 'graph_adjacency.tsv'):
        if not (inp / name).is_file():
            raise FileNotFoundError(f'Retain the supplied configuration first: {inp / name}')
    areas = list(csv.DictReader((inp / 'area_order.tsv').open(), delimiter='\t'))
    assert [r['code'] for r in areas] == list('ABCEFGHIJK')
    states = [r['state_code'] for r in csv.DictReader((inp / 'allowed_states.tsv').open(), delimiter='\t')]
    assert len(states) == 45 and max(len(s.split('+')) for s in states) == 3
    geography = (inp / 'geography.data').read_text().splitlines()
    tips = {l.split()[0] for l in geography[1:] if l.strip()}
    assert len(tips) == 29
    mcc = (summary / f'{calibration}_unconstrained_final_mcc_29tip.newick').read_text().strip()
    posterior = [l.strip() for l in (summary / f'{calibration}_unconstrained_100_posterior_29tip.newick').read_text().splitlines() if l.strip()]
    assert len(posterior) == 100
    for tree in [mcc] + posterior:
        assert {t.name for t in Phylo.read(io.StringIO(tree), 'newick').get_terminals()} == tips
    # Pruning the outgroup leaves a root stem in the upstream export. It is
    # outside the ingroup crown and was excluded from all reported range fits.
    def crown_only(tree):
        return re.sub(r':[0-9.eE+\-]+;$', ';', tree)
    mcc=crown_only(mcc)
    posterior=[crown_only(t) for t in posterior]
    (inp / 'mcc_myr.newick').write_text(mcc + '\n')
    for i, tree in enumerate(posterior, 1):
        (inp / f'posterior_{i:03d}_myr.newick').write_text(tree + '\n')
    (inp / 'max_range_size.txt').write_text('3\n')
    print(f'{calibration}: 29 tips, ten areas, 45 allowed states, MCC + 100 posterior trees')

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--calibration', choices=['givnish', 'zhang', 'both'], default='both')
    p.add_argument('--output-root', type=Path, help='Testing copy containing the supplied configuration')
    args = p.parse_args()
    for calibration in (('givnish', 'zhang') if args.calibration == 'both' else (args.calibration,)):
        prepare(calibration, args.output_root)

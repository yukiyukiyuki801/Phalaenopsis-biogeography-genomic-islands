#!/usr/bin/env python3
"""Read-only integrity, syntax, and key dataset checks for the manuscript export."""
from pathlib import Path
import ast
import hashlib
import subprocess
import xml.etree.ElementTree as ET
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
for p in (ROOT/'orchid_paper2/scripts').glob('*.py'):ast.parse(p.read_text(),filename=str(p))
for p in (ROOT/'orchid_paper2/scripts').iterdir():
    if p.suffix in ('.sh','.sbatch'):subprocess.run(['bash','-n',str(p)],check=True)
for cal in ('givnish','zhang'):
    p=ROOT/f'orchid_paper2/16_rad_beast/xml/rad29_peq_{cal}_combined_rateprior_unconstrained_60m.xml'
    text=p.read_text();ET.fromstring(text)
    assert 'pamMonophyly' not in text and 'papMonophyly' not in text
    assert 'chainLength="60000000"' in text
    assert 'ingroupCalibration' in text and 'rootCalibration' in text
v=ROOT/'phalaenopsis_mac20260730/Phalaenopsis_eq_LD_filter.vcf'
sites=0
for line in v.open():
    if line.startswith('#CHROM'):assert len(line.rstrip().split('\t')[9:])==170
    elif not line.startswith('#'):sites+=1
assert sites==6596,sites
counts=pd.read_csv(ROOT/'tables/Table_S3.tsv',sep='\t')
assert len(counts)==18 and counts.n_individuals.sum()==170
g=pd.read_csv(ROOT/'orchid_paper2/19_local_genealogy/results/genealogy_class_permutation_tests.tsv',sep='\t')
assert len(g)==12 and set(g.analysis)=={'all_650_phase_free'}
assert len(list((ROOT/'orchid_paper2/19_local_genealogy/trees_all650_diploid_ml').glob('*.treefile')))==650
models=pd.read_csv(ROOT/'orchid_paper2/22_diversification/summary/rpanda_model_comparison.tsv',sep='\t')
assert set(models.model)=={'N1','E1','E2','E3','E4','T1','T2','T3','T4'}
if (ROOT/'checksums.sha256').exists():
    for line in (ROOT/'checksums.sha256').read_text().splitlines():
        checksum,path=line.split('  ',1)
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==checksum,path
print('PASS: Python/shell syntax, final XMLs, 170 samples/6596 SNPs, 18 units, 650 trees, 12 genealogy tests, nine RPANDA models, and available checksums.')

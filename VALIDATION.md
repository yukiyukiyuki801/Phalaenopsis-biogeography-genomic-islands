# Release validation

This is a curated research-workflow export, not a newly executed analysis.

Completed checks:

- Python AST parsing and Bash syntax checks for packaged analysis scripts.
- Well-formed final BEAST XMLs: 60 million states, two age calibrations,
  broad clock prior, and no PAM or PAP monophyly constraint.
- Published VCF: 170 individuals and 6,596 LD-pruned SNPs.
- Migration assignments: 18 units totaling 170 individuals.
- Genealogy archive: 650 ML trees and the original 12-contrast permutation family.
- RPANDA comparison: nine retained models, N1, E1–E4 and T1–T4.
- Regenerated individual network distances agree with the saved matrix to
  maximum absolute difference 1.12 × 10⁻¹⁶.
- The packaged migration converter reproduces 6,596 SNPs, 170 individuals,
  and 18 populations without additional frequency filtering.
- Source/adaptation hashes and release-file SHA-256 manifest.
- Manuscript-scope checks: no Herodotools raster/inference scripts or old
  constrained-tree/SNAP(P)ER input paths in executable workflows. Shared
  bioregion polygons and assignments are retained under BioGeoBEARS.
- RAD sample assignment: 94 selected individuals and 29 population consensuses.
- BioGeoBEARS inputs: ten sampled areas, 45 allowed states including null,
  maximum range size three. Input-refresh helper checked in a temporary copy;
  it retains the fitted geography and uses only the unconstrained trees.

Run `python scripts/check_package.py` to repeat the integrity checks. The
checker does not submit jobs, install software, or fit models.

Limitations:

- Long MCMC chains, orthology inference, alignments, full variant calling and
  statistical model fits were not rerun during packaging.
- Batch examples retain institution-specific modules and scheduler settings;
  adapt these and install dependencies before execution.
- Large sequencing inputs and external software are not shipped. Public
  accessions and preparation scripts identify their sources.
- The migration summary uses approximate lag-one pooled ESS and ordinary
  R-hat, not rank-normalized bulk ESS. These summaries do not guarantee mixing.
- Saved final summaries and manuscript table inputs are supplied to support
  inspection; the export has not been tested end-to-end in a fresh environment.
- Manuscript figure assembly and annotated Word files are intentionally not
  uploaded. Analysis plotting scripts are supplied where retained; no claim
  is made that all historical plotting edits can be recreated in one command.

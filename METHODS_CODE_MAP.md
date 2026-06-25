# Methods-to-Code Map

This file maps the manuscript Methods sections to the curated scripts in this repository. It is intended to make the public code package transparent without redistributing source datasets or large intermediate files.

## Study datasets and preprocessing

Manuscript methods:

- Study design and participants
- DNAm sample and measurement
- Quality control and normalization

Code:

- `scripts/01_discovery_preprocessing/01_sgpd_load_idat_original.R`
- `scripts/01_discovery_preprocessing/02_sgpd_champ_hc.R`
- `scripts/01_discovery_preprocessing/03_sgpd_champ_pd.R`
- `scripts/01_discovery_preprocessing/05_peg1_load_idat_original.R`
- `scripts/01_discovery_preprocessing/06_peg1_champ_hc.R`
- `scripts/01_discovery_preprocessing/07_peg1_champ_pd.R`
- `scripts/01_discovery_preprocessing/10_peg2_saliva_champ_hc.R`
- `scripts/01_discovery_preprocessing/11_peg2_saliva_champ_pd.R`
- `scripts/05_reviewer_requested_analyses/03_r2_9_reestimate_cell_fractions.R`

Covered steps include raw IDAT loading, sample-level QC, sex prediction with `minfi::getSex`, probe filtering, BMIQ normalization, ComBat batch correction, and ChAMP Refbase cell-composition adjustment.

## Discovery-stage sex-DMP analysis and inflation correction

Manuscript methods:

- Sex difference analysis of individual datasets
- Bacon correction
- QQ plots and genomic inflation factors

Code:

- `scripts/01_discovery_preprocessing/02_sgpd_champ_hc.R`
- `scripts/01_discovery_preprocessing/03_sgpd_champ_pd.R`
- `scripts/01_discovery_preprocessing/06_peg1_champ_hc.R`
- `scripts/01_discovery_preprocessing/07_peg1_champ_pd.R`
- `scripts/01_discovery_preprocessing/10_peg2_saliva_champ_hc.R`
- `scripts/01_discovery_preprocessing/11_peg2_saliva_champ_pd.R`

Covered steps include `champ.DMP`, empirical-null correction with `bacon`, and output of bacon-corrected effect sizes, standard errors, p-values, and FDR values.

## Meta-analysis

Manuscript methods:

- Meta-analysis for sex-DMPs

Code:

- `scripts/02_meta_analysis/01_prepare_metal_inputs.R`
- `scripts/02_meta_analysis/02_postprocess_metal_outputs.R`
- `scripts/02_meta_analysis/03_metal_pd_protocol.txt`
- `scripts/02_meta_analysis/04_metal_hc_protocol.txt`
- `scripts/02_meta_analysis/05_ppmi_saliva_meta_postprocess.R`

The METAL executable is not redistributed. Users should install METAL separately and run the provided protocol files.

## DMR analysis

Manuscript methods:

- Identification of sex-differential methylation regions
- DMRcate analysis
- comb-p analysis
- Union-based definition of final PD-unique DMRs

Code:

- `scripts/01_discovery_preprocessing/04_sgpd_dmr.R`
- `scripts/01_discovery_preprocessing/08_peg1_dmr.R`
- `scripts/01_discovery_preprocessing/12_peg2_saliva_dmr.R`
- `scripts/03_dmp_dmr_downstream/03_prepare_bed_for_combp.R`
- `scripts/03_dmp_dmr_downstream/04_run_combp.py`
- `scripts/03_dmp_dmr_downstream/05_dmr_combp.R`
- `scripts/03_dmp_dmr_downstream/06_dmr_dmp_dmr_cate4_downstream.R`
- `scripts/03_dmp_dmr_downstream/07_dmr_bumphunter_dmr_cate7_downstream.R`

The comb-p package is not redistributed. Users should install comb-p separately.

## Functional enrichment and network analysis

Manuscript methods:

- GO enrichment analysis
- DMP-set background comparison with missMethyl
- PPI and hub-gene analysis

Code:

- `scripts/03_dmp_dmr_downstream/01_dmp_downstream.R`
- `scripts/03_dmp_dmr_downstream/02_dmp_kegg_go_up_down.R`
- `scripts/03_dmp_dmr_downstream/05_dmr_combp.R`
- `scripts/03_dmp_dmr_downstream/06_dmr_dmp_dmr_cate4_downstream.R`
- `scripts/05_reviewer_requested_analyses/01_r2_1_dmp_background_comparison.R`

STRING and Cytoscape/cytoHubba analyses require use of external web or desktop tools and are therefore not fully automated in this repository. Exported PPI and hub-gene result files are not redistributed here.

## PPMI validation and longitudinal analyses

Manuscript methods:

- Validation in the PPMI cohort
- Longitudinal DNAm and RNA time-series analyses
- Sex-based EWAS and TWAS in PPMI

Code:

- `scripts/04_validation_and_integration/00_ppmi_bl_champ_pd_pipeline.R`
- `scripts/04_validation_and_integration/01_ppmi_time_series_dnam_2199cpg.R`
- `scripts/04_validation_and_integration/02_ppmi_time_series_dnam_536cpg.R`
- `scripts/04_validation_and_integration/03_rna_time_series_459gene.R`
- `scripts/04_validation_and_integration/04_rna_time_series_2199rna.R`

PPMI source data must be downloaded from the PPMI data portal under its access policy. The baseline DNAm preprocessing and EWAS pipeline is provided as a representative script; other PPMI methylation batches were processed using the same pipeline with visit-specific sample sheets and input directories. The longitudinal validation scripts use the resulting PPMI EWAS and TWAS output files as inputs for validation and clustering.

## eQTM and CIT analyses

Manuscript methods:

- MatrixEQTL cis- and trans-eQTM analyses
- Age-adjusted eQTM analysis
- Age plus blood cell-composition-adjusted sensitivity eQTM analysis
- Causal inference test

Code:

- `scripts/04_validation_and_integration/05a_prepare_dnam_inputs_for_eqtm.R`
- `scripts/04_validation_and_integration/05_eqtm_matrixeqtl.R`
- `scripts/04_validation_and_integration/06_znf727_plots.R`
- `scripts/04_validation_and_integration/07_cit_example.R`
- `scripts/04_validation_and_integration/08_cit_clinical_traits.R`
- `scripts/05_reviewer_requested_analyses/02_r2_9_check_cell_fraction_absence.R`
- `scripts/05_reviewer_requested_analyses/04_r2_9_eqtm_age_cell_sensitivity.R`

## Cross-tissue validation

Manuscript methods:

- Saliva-based replication
- In silico blood-brain methylation concordance analysis

Code:

- `scripts/01_discovery_preprocessing/10_peg2_saliva_champ_hc.R`
- `scripts/01_discovery_preprocessing/11_peg2_saliva_champ_pd.R`
- `scripts/01_discovery_preprocessing/13_peg2_saliva_cross_tissue.R`
- `scripts/02_meta_analysis/05_ppmi_saliva_meta_postprocess.R`

The blood-brain methylation concordance checks used external online resources: the Blood Brain DNA Methylation Comparison Tool and BECon. These web-tool queries are not redistributed as executable code.

## Peer-review sensitivity analyses

Code:

- `scripts/05_reviewer_requested_analyses/01_r2_1_dmp_background_comparison.R`
- `scripts/05_reviewer_requested_analyses/02_r2_9_check_cell_fraction_absence.R`
- `scripts/05_reviewer_requested_analyses/03_r2_9_reestimate_cell_fractions.R`
- `scripts/05_reviewer_requested_analyses/04_r2_9_eqtm_age_cell_sensitivity.R`

These scripts correspond to analyses added during revision in response to reviewer comments.

# Epigenome-wide sex-differential DNA methylation analysis in Parkinson's disease

This repository contains curated analysis scripts for the manuscript:

Epigenome-wide meta-analysis of sex-differential DNA methylation in Parkinson's disease.

The scripts were curated from the authors' analysis folder for transparency and reproducibility. Public source data are not redistributed here. Users should obtain the underlying methylation, transcriptomic, and clinical data from the original cohort repositories, including GEO records for SGPD, PEG1, and PEG2, and the PPMI data portal under the applicable data-use terms.

## Repository contents

- `scripts/01_discovery_preprocessing/`: methylation array preprocessing, quality control, normalization, batch correction, cell-composition correction, sex-DMP analysis, and DMR analysis for SGPD, PEG1, and PEG2.
- `scripts/02_meta_analysis/`: input preparation and output post-processing for METAL meta-analysis, plus METAL protocol files.
- `scripts/03_dmp_dmr_downstream/`: DMP/DMR downstream annotation, enrichment, comb-p input preparation, and comb-p wrapper script.
- `scripts/04_validation_and_integration/`: PPMI baseline DNAm preprocessing/EWAS pipeline, longitudinal DNAm validation, RNA time-series analyses, eQTM input preparation and analysis, ZNF727 visualization, and CIT analyses.
- `scripts/05_reviewer_requested_analyses/`: additional analyses performed during peer review, including DMP background comparisons and age plus cell-composition-adjusted eQTM sensitivity analyses.
- `METHODS_CODE_MAP.md`: mapping between manuscript Methods sections and the corresponding scripts.
- `CODE_MANIFEST.tsv`: public script inventory with workflow labels, SHA-256 checksums, and file sizes.
- `PUBLIC_CODE_CHECKSUMS.tsv`: checksums for all public files in this repository.

## Data and external software

Input data and large intermediate files are not included. The scripts expect the user to place downloaded and preprocessed data files in the working directory or adjust file paths accordingly.

The PPMI baseline DNAm preprocessing and EWAS pipeline is provided as a representative script. Other PPMI methylation batches were processed using the same pipeline with visit-specific sample sheets and input directories. The longitudinal validation scripts use the resulting PPMI EWAS output files as inputs.

External tools used by the workflow include:

- R and Bioconductor packages including `ChAMP`, `minfi`, `limma`, `sva`, `bacon`, `DMRcate`, `missMethyl`, `clusterProfiler`, `MatrixEQTL`, `CIT`, and plotting/data-processing packages.
- METAL for inverse-variance fixed-effect meta-analysis. The METAL executable is not redistributed in this repository.
- comb-p for spatially correlated p-value region analysis. Install comb-p separately before running the comb-p wrapper script.
- Cytoscape/cytoHubba and STRING were used for PPI and hub-gene analyses; these steps may require manual export from those tools.

## Suggested workflow

1. Obtain raw IDAT files and phenotype/sample-sheet files from the original public cohorts.
2. Run the relevant scripts in `scripts/01_discovery_preprocessing/` for each cohort and subgroup.
3. Prepare and run METAL using scripts and protocol files in `scripts/02_meta_analysis/`.
4. Run DMP and DMR downstream analyses in `scripts/03_dmp_dmr_downstream/`.
5. Run PPMI preprocessing/EWAS, validation, RNA, eQTM, ZNF727, and CIT scripts in `scripts/04_validation_and_integration/`.
6. Run peer-review sensitivity analyses in `scripts/05_reviewer_requested_analyses/`.

Some scripts retain local file names from the original analysis environment. Before rerunning, edit input/output paths to match the local directory structure and confirm that package versions are compatible with the analysis environment.

## License

No license has been selected for this repository. Please contact the corresponding author for reuse beyond review and reproducibility purposes.

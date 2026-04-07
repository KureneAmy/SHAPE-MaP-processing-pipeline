# SHAPE-MaP Analysis Report

---

**Report Date:** 2026-04-07  
**Analysis Date:** 2026/04/07  
**Institution:** BIG  
**Principal Investigator:** HYJ  
**Project ID:** report-test  
**Output Directory:** `/mnt1/4.NAS2025/zhangam/easyomics_test/report/SHAPE-MaP/output`

---

## Table of Contents

1. [Project Information](#1-project-information)
2. [Analysis Overview](#2-analysis-overview)
3. [Data Quality Assessment](#3-data-quality-assessment)
4. [SHAPE Reactivity Analysis](#4-shape-reactivity-analysis)
5. [ShapeMapper Results Summary](#5-shapemapper-results-summary)
6. [RNA Structure Prediction](#6-rna-structure-prediction)
7. [QC Metrics](#7-qc-metrics)
8. [Methods](#8-methods)
9. [Contact Information](#9-contact-information)

---

## 1. Project Information

| Field              | Value                          |
|--------------------|--------------------------------|
| Report Title       | SHAPE-MaP Analysis Report             |
| Report Date        | 2026-04-07              |
| Analysis Date      | 2026/04/07            |
| Institution        | BIG     |
| Principal Investigator | HYJ     |
| Project ID         | report-test      |
| Number of Samples  | 1             |
| Output Directory   | `/mnt1/4.NAS2025/zhangam/easyomics_test/report/SHAPE-MaP/output`             |
| MultiQC Report     | `/mnt1/4.NAS2025/zhangam/easyomics_test/report/SHAPE-MaP/output/multiqc_report.html`         |

---

## 2. Analysis Overview

This report summarises the results of a **SHAPE-MaP** (Selective 2′-Hydroxyl Acylation
analysed by Primer Extension and Mutational Profiling) experiment processed with
**ShapeMapper2** and the **EasyOmics SHAPE-MaP Processing Pipeline**.

A total of **1** sample(s) were processed. SHAPE reactivity values
were computed and used for RNA secondary-structure prediction with **RNAstructure**.

### Samples Analysed


- **test** — 1 target(s):
    - `Hs_DRAIC_ncRNA`
  


---

## 3. Data Quality Assessment


### Sample: test


#### Target: Hs_DRAIC_ncRNA

**Read Processing Summary**

| Metric | Value |
|--------|-------|
| Overall Alignment Rate | 29.71% |
| Concordantly Mapped Reads (1×) | 22,148 |
| Mean Read Depth | 76114.7 |


**Per-Channel Quality Control**

| Channel | Total Read Depth | Total Mutations | Median Mutation Rate | Mean Read Depth |
|---------|-----------------|-----------------|---------------------|-----------------|
| Modified | 132,059,024 | 1,233,108 | 0.0058 | 76114.7 |
| Untreated | 138,850,709 | 579,155 | 0.0000 | 80029.2 |
| Denatured | 0 | 0 | N/A | 0.0 |





---

## 4. SHAPE Reactivity Analysis


### Sample: test


#### Target: Hs_DRAIC_ncRNA

**Mutation Rates by Channel**

| Channel | Median Mutation Rate | Total Mutations | Total Read Depth |
|---------|---------------------|-----------------|-----------------|
| Modified | 0.0058 | 1,233,108 | 132,059,024 |
| Untreated | 0.0000 | 579,155 | 138,850,709 |
| Denatured | N/A | 0 | 0 |
| High-Quality Profile (%) | 27.7% | — | — |


**Reactivity Summary**


| Metric | Value |
|--------|-------|
| Nucleotides with Valid Reactivity | 480 |
| Mean Reactivity | 0.4494 |
| Median Reactivity | 0.2664 |
| Std Dev | 0.5502 |
| Maximum Reactivity | 3.0408 |
| Reactive Fraction (> 0.4) | 41.7% |





---

## 5. ShapeMapper Results Summary

| Sample | Target | Alignment Rate | Mean Depth | Mut. Rate (mod.) | Mut. Rate (untr.) | Reactive Fraction |
|--------|--------|---------------|------------|------------------|-------------------|-------------------|
| test | Hs_DRAIC_ncRNA | 29.71% | 76114.7 | 0.0058 | 0.0000 | 41.7% |


---

## 6. RNA Structure Prediction

RNA secondary structures were predicted using **RNAstructure Fold** with
SHAPE-directed folding constraints.

### Method

SHAPE reactivity values from each sample were incorporated as pseudo-free-energy
constraints to guide the thermodynamic folding algorithm. The minimum free energy (MFE)
structure is reported.

### Output Files


#### Sample: test


**Target: Hs_DRAIC_ncRNA**

- CT file: `/mnt1/4.NAS2025/zhangam/easyomics_test/report/SHAPE-MaP/output/test/Hs_DRAIC_ncRNA/structure/test_Hs_DRAIC_ncRNA.ct`
- DBN file: `/mnt1/4.NAS2025/zhangam/easyomics_test/report/SHAPE-MaP/output/test/Hs_DRAIC_ncRNA/structure/test_Hs_DRAIC_ncRNA.dbn`
- SVG visualisation: `/mnt1/4.NAS2025/zhangam/easyomics_test/report/SHAPE-MaP/output/test/Hs_DRAIC_ncRNA/structure/test_Hs_DRAIC_ncRNA_folding.svg`




---

## 7. QC Metrics

### MultiQC Report

A comprehensive MultiQC report aggregating FastQC results for all samples is available at:

```
/mnt1/4.NAS2025/zhangam/easyomics_test/report/SHAPE-MaP/output/multiqc_report.html
```

### Sample QC Summary

| Sample | Target | Total Input Pairs | Pairs Merged | Alignment Rate | Mean Read Depth | Status |
|--------|--------|-------------------|--------------|----------------|-----------------|--------|
| test | Hs_DRAIC_ncRNA | N/A | N/A | 29.71% | 76114.7 | LOW |


---

## 8. Methods

### Analysis Pipeline

SHAPE-MaP experiments were processed using the **EasyOmics SHAPE-MaP Processing Pipeline**
based on Snakemake. The pipeline encompasses the following major steps:

1. **Quality Control** – Raw reads were assessed with FastQC and aggregated with MultiQC.
2. **Reference Splitting** – Multi-target reference FASTA files were split into individual sequences.
3. **ShapeMapper2** – Modified, untreated (and optionally denatured) read pairs were jointly
   processed to derive per-nucleotide SHAPE reactivity values.
4. **RNA Structure Prediction** – SHAPE-constrained secondary-structure prediction was performed
   with RNAstructure `Fold`.
5. **Structure Visualisation** – Arc diagrams in SVG format were generated with RNAstructure `draw`.

### Software

| Software | Version | Purpose |
|----------|---------|---------|
| Snakemake | — | Workflow management |
| FastQC | — | Read-level quality control |
| MultiQC | — | QC report aggregation |
| BBMerge | — | Paired-end read merging |
| Bowtie2 | — | Read alignment |
| ShapeMapper2 | 2.3.0 | SHAPE reactivity computation |
| RNAstructure | — | RNA secondary-structure prediction |

*Version numbers are extracted automatically from ShapeMapper log files where available.*

### References

1. Busan S, Weeks KM. Accurate detection of chemical modifications in RNA by mutational profiling
   (MaP) with ShapeMapper 2. *RNA*. 2018;24(2):143–148.
2. Reuter JS, Mathews DH. RNAstructure: software for RNA secondary structure prediction and
   analysis. *BMC Bioinformatics*. 2010;11:129.
3. Ewels P, et al. MultiQC: summarize analysis results for multiple tools and samples in a single
   report. *Bioinformatics*. 2016;32(19):3047–3048.

---

## 9. Contact Information

For questions regarding this analysis, please contact:

- **Institution:** BIG
- **Principal Investigator:** HYJ
- **Pipeline:** EasyOmics SHAPE-MaP Processing Pipeline

---

*This report was automatically generated by `scripts/compile_report.py` on 2026-04-07.*

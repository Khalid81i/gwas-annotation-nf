# gwas-annotation-nf

[![Nextflow](https://img.shields.io/badge/Nextflow-%E2%89%A524.0-23aa62?style=flat-square&logo=nextflow)](https://www.nextflow.io/)
[![Docker](https://img.shields.io/badge/Run_with-Docker-2496ED?style=flat-square&logo=docker)](https://www.docker.com/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg?style=flat-square)](https://opensource.org/licenses/MIT)

A Nextflow (DSL2) pipeline for annotating a sample VCF against the GWAS Catalog,
computing effect-allele-aligned genotype dosages, and summarising results by trait.

Originally developed in Bash during a clinical bioinformatics internship at
**NoorDX** (Jeddah, Saudi Arabia) and ported to Nextflow as an open-source,
reproducible, containerised pipeline.

---

## Overview

Given a single-sample VCF and a pre-downloaded GWAS Catalog TSV, the pipeline:

1. **Normalises** the VCF — splits multiallelic sites, assigns `CHROM:POS:REF:ALT` IDs, indexes output
2. **Filters** the GWAS Catalog to genome-wide significant associations (p ≤ 5×10⁻⁸)
3. **Matches** sample variants to GWAS hits by chromosome + position
4. **Computes** effect-allele-aligned genotype dosage (0/1/2), flipping dosage when the effect allele is REF rather than ALT
5. **Summarises** matched variants by trait with per-trait variant counts, total dosages, and a dosage-weighted mean effect size

## Pipeline architecture

```
  sample.vcf.gz          gwas_catalog.tsv
       |                        |
       v                        v
  NORMALIZE_VCF            FILTER_GWAS
       |                        |
       +----------+-------------+
                  v
           MATCH_VARIANTS
                  |
                  v
           COMPUTE_DOSAGE
                  |
                  v
          SUMMARISE_TRAITS
                  |
                  v
   annotated_variants.tsv + trait_summary.tsv
```

Every process runs inside a pinned Docker container, so results are reproducible
across any machine with Docker and Nextflow installed.

## Quick start

### Requirements

- [Nextflow](https://www.nextflow.io/) >= 24.0
- [Docker](https://www.docker.com/) (or Docker Engine on Linux)

### Clone and run

```bash
git clone https://github.com/<your-username>/gwas-annotation-nf.git
cd gwas-annotation-nf

nextflow run main.nf \
    --vcf  data/test/match_test.vcf.gz \
    --gwas data/test/mock_gwas.tsv
```

On first run, Nextflow will pull the required Docker containers
(`bcftools`, `gawk`, `python`). Subsequent runs reuse the cached images.

### Parameters

| Parameter | Description                                                   | Required |
|-----------|---------------------------------------------------------------|----------|
| `--vcf`   | Input VCF file (`.vcf.gz`)                                    | Yes      |
| `--gwas`  | GWAS Catalog TSV download                                     | Yes      |
| `--outdir`| Results directory (default: `results/`)                       | No       |

## Input formats

### VCF

A standard bgzipped VCF for a single sample. For real data, use human VCFs
aligned to a recent reference build (hg38 recommended). The pipeline assumes
chromosome naming without the `chr` prefix; pre-normalise if your VCF uses
`chr1`, `chr2`, etc.

### GWAS Catalog

The full tab-separated download from the
[NHGRI-EBI GWAS Catalog](https://www.ebi.ac.uk/gwas/docs/file-downloads).
Key columns used: `CHR_ID`, `CHR_POS`, `STRONGEST SNP-RISK ALLELE`, `P-VALUE`,
`OR or BETA`, `DISEASE/TRAIT`. Column indices are detected dynamically so the
pipeline will tolerate column reorderings between releases.

## Outputs

All outputs land in `results/` (configurable via `--outdir`):

```
results/
├── normalised/
│   ├── <sample>.norm.vcf.gz
│   └── <sample>.norm.vcf.gz.tbi
├── gwas_filtered/
│   └── gwas_significant.tsv
├── matched/
│   └── matched_variants.tsv
├── annotated/
│   └── annotated_variants.tsv
├── summary/
│   └── trait_summary.tsv
└── reports/
    └── execution_report.html
```

### `annotated_variants.tsv`

One row per variant that matched a significant GWAS hit.

| Column                 | Description                                              |
|------------------------|----------------------------------------------------------|
| `CHR`                  | Chromosome                                               |
| `POS`                  | Position (1-based)                                       |
| `REF`                  | Reference allele                                         |
| `ALT`                  | Alternate allele                                         |
| `GT`                   | Sample genotype (e.g. `0/1`, `1/1`)                      |
| `SNP_RISK_ALLELE`      | Raw GWAS Catalog field (e.g. `rs123-A`)                  |
| `OR_BETA`              | Effect size from the GWAS Catalog                        |
| `TRAIT`                | Reported disease or trait                                |
| `EFFECT_ALLELE`        | Parsed effect allele letter                              |
| `EFFECT_ALLELE_DOSAGE` | 0, 1, or 2 — number of effect alleles the sample carries |

### `trait_summary.tsv`

One row per trait with a compact summary across all matched variants:

| Column                      | Description                                                       |
|-----------------------------|-------------------------------------------------------------------|
| `TRAIT`                     | Trait name                                                        |
| `N_VARIANTS`                | Number of matched significant variants                            |
| `TOTAL_DOSAGE`              | Sum of effect-allele dosages across matched variants              |
| `DOSAGE_WEIGHTED_MEAN_BETA` | Dosage-weighted mean effect size (beta or ln-OR as reported)      |

## Effect-allele alignment logic

The core conceptual operation is in `COMPUTE_DOSAGE`. The VCF genotype is always
relative to (REF, ALT), but the GWAS Catalog effect size is relative to the
*risk/effect allele*, which may be either REF or ALT. Dosage is aligned so that
the reported value always counts effect alleles:

| GWAS effect allele | VCF genotype | Raw ALT dosage | Effect-allele dosage |
|--------------------|--------------|----------------|----------------------|
| = ALT              | 0/0          | 0              | 0                    |
| = ALT              | 0/1          | 1              | 1                    |
| = ALT              | 1/1          | 2              | 2                    |
| = REF              | 0/0          | 0              | 2                    |
| = REF              | 0/1          | 1              | 1                    |
| = REF              | 1/1          | 2              | 0                    |
| != REF and != ALT  | any          | any            | `NA` (flagged)       |

Without this alignment, downstream analyses (polygenic scoring, trait
prediction) would silently invert the direction of effect for roughly half the
variants in a typical dataset.

## Project structure

```
gwas-annotation-nf/
├── main.nf                        Top-level workflow
├── nextflow.config                Pipeline defaults and Docker enablement
├── modules/
│   └── local/
│       ├── normalize_vcf.nf       bcftools norm + ID assignment
│       ├── filter_gwas.nf         p-value filter + column extraction
│       ├── match_variants.nf      intersect VCF with filtered GWAS
│       ├── compute_dosage.nf      genotype -> effect-allele-aligned dosage
│       └── summarise_traits.nf    aggregate by trait
├── data/
│   └── test/                      Minimal test VCF + mock GWAS TSV
└── results/                       Created on first run
```

## Background

This pipeline reimplements, in production-grade Nextflow DSL2, the variant
annotation and GWAS interpretation logic I developed during a clinical
bioinformatics internship at NoorDX (Jeddah, Saudi Arabia). The original Bash
implementation processed real diagnostic sequencing data; this open-source
version uses only public test data and is safe to share, extend, and deploy.
A detailed technical case study of the original NoorDX work — including the AncestryG81 (population genetics, PLINK 2.0 + ADMIXTURE) and PharmaGB (pharmacogenomics) layers, engineering challenges encountered, and what I would do differently — is available at [`NOORDX_CASE_STUDY.md`](./NOORDX_CASE_STUDY.md).

Key improvements over the original Bash pipeline:

- **Reproducibility** — every tool runs in a pinned container
- **Portability** — identical execution on a laptop, HPC cluster, or cloud
- **Observability** — execution reports, resumability, per-process tagging
- **Modularity** — each step is an independent DSL2 module, testable in isolation

## Future work

- Docker-pinned reference genome handling for strict build verification
- Optional sample-VCF pre-liftover between reference builds
- Integration with [PGS Catalog](https://www.pgscatalog.org/) scoring files for
  direct polygenic score computation
- nf-core-style `schema.json` parameter validation
- CI test with GitHub Actions running `-profile test`

## Author

**Khalid Dallol**
BEng (Hons) Electrical and Electronic Engineering, University of Greenwich.
Final-year student with clinical bioinformatics experience at NoorDX.

## License

MIT — see [LICENSE](LICENSE) for details.

# gwas-annotation-nf

A Nextflow (DSL2) pipeline for annotating sample VCF files against the GWAS Catalog.

Ported from a clinical bash pipeline originally developed at NoorDX (Jeddah, Saudi Arabia).

## Overview

Given a single-sample VCF and a GWAS Catalog download, this pipeline:

1. Normalises the VCF (splits multiallelics, standardises chromosome naming, assigns `CHROM:POS:REF:ALT` IDs)
2. Filters the GWAS Catalog to genome-wide significant hits (p ≤ 5×10⁻⁸)
3. Intersects sample variants with GWAS hits on chromosome + position
4. Computes allele dosage with effect-allele alignment
5. Produces a per-variant annotated TSV and a trait-level summary

## Status

In development.

## Author

Khalid Dallol

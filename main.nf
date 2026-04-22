#!/usr/bin/env nextflow

include { NORMALIZE_VCF }  from './modules/local/normalize_vcf.nf'
include { FILTER_GWAS }    from './modules/local/filter_gwas.nf'
include { MATCH_VARIANTS } from './modules/local/match_variants.nf'

workflow {

    // Validate inputs
    if (!params.vcf)  error "Please provide a VCF file with --vcf"
    if (!params.gwas) error "Please provide a GWAS Catalog TSV with --gwas"

    // Input channels
    vcf_ch  = Channel.fromPath(params.vcf,  checkIfExists: true)
    gwas_ch = Channel.fromPath(params.gwas, checkIfExists: true)

    // Step 1: normalise VCF
    NORMALIZE_VCF(vcf_ch)

    // Step 2: filter GWAS Catalog
    FILTER_GWAS(gwas_ch)

    // Step 3: intersect VCF with filtered GWAS
    MATCH_VARIANTS(
        NORMALIZE_VCF.out.vcf,
        NORMALIZE_VCF.out.tbi,
        FILTER_GWAS.out.tsv
    )

    // View outputs
    NORMALIZE_VCF.out.vcf.view  { "Normalised VCF: $it" }
    FILTER_GWAS.out.tsv.view    { "Filtered GWAS:  $it" }
    MATCH_VARIANTS.out.tsv.view { "Matched:        $it" }
}

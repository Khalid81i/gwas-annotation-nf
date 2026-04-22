#!/usr/bin/env nextflow

include { NORMALIZE_VCF } from './modules/local/normalize_vcf.nf'
include { FILTER_GWAS }   from './modules/local/filter_gwas.nf'

workflow {

    // Validate inputs
    if (!params.vcf)  error "Please provide a VCF file with --vcf"
    if (!params.gwas) error "Please provide a GWAS Catalog TSV with --gwas"

    // Create channels
    vcf_ch  = Channel.fromPath(params.vcf,  checkIfExists: true)
    gwas_ch = Channel.fromPath(params.gwas, checkIfExists: true)

    // Run processes
    NORMALIZE_VCF(vcf_ch)
    FILTER_GWAS(gwas_ch)

    // View outputs
    NORMALIZE_VCF.out.vcf.view { "Normalised VCF: $it" }
    FILTER_GWAS.out.tsv.view   { "Filtered GWAS: $it" }
}

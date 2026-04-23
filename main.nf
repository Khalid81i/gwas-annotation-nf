#!/usr/bin/env nextflow

include { NORMALIZE_VCF }  from './modules/local/normalize_vcf.nf'
include { FILTER_GWAS }    from './modules/local/filter_gwas.nf'
include { MATCH_VARIANTS } from './modules/local/match_variants.nf'
include { COMPUTE_DOSAGE } from './modules/local/compute_dosage.nf'

workflow {

    if (!params.vcf)  error "Please provide a VCF file with --vcf"
    if (!params.gwas) error "Please provide a GWAS Catalog TSV with --gwas"

    vcf_ch  = Channel.fromPath(params.vcf,  checkIfExists: true)
    gwas_ch = Channel.fromPath(params.gwas, checkIfExists: true)

    NORMALIZE_VCF(vcf_ch)
    FILTER_GWAS(gwas_ch)
    MATCH_VARIANTS(
        NORMALIZE_VCF.out.vcf,
        NORMALIZE_VCF.out.tbi,
        FILTER_GWAS.out.tsv
    )
    COMPUTE_DOSAGE(MATCH_VARIANTS.out.tsv)

    NORMALIZE_VCF.out.vcf.view   { "Normalised VCF: $it" }
    FILTER_GWAS.out.tsv.view     { "Filtered GWAS:  $it" }
    MATCH_VARIANTS.out.tsv.view  { "Matched:        $it" }
    COMPUTE_DOSAGE.out.tsv.view  { "Annotated:      $it" }
}

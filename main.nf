#!/usr/bin/env nextflow

include { NORMALIZE_VCF } from './modules/local/normalize_vcf.nf'

workflow {

    // Validate required inputs
    if (!params.vcf) {
        error "Please provide a VCF file with --vcf"
    }

    // Create a channel from the input VCF
    vcf_ch = Channel.fromPath(params.vcf, checkIfExists: true)

    // Run normalisation
    NORMALIZE_VCF(vcf_ch)

    // Print output path on completion
    NORMALIZE_VCF.out.vcf.view { "Normalised VCF: $it" }
}

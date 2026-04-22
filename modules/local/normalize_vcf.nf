process NORMALIZE_VCF {

    tag "$vcf.baseName"
    publishDir "${params.outdir}/normalised", mode: 'copy'

    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_0'

    input:
    path vcf

    output:
    path "${vcf.baseName}.norm.vcf.gz",     emit: vcf
    path "${vcf.baseName}.norm.vcf.gz.tbi", emit: tbi

    script:
    """
    bcftools norm \\
        -m -any \\
        -Oz \\
        -o ${vcf.baseName}.norm.vcf.gz \\
        ${vcf}

    # Assign CHROM:POS:REF:ALT style IDs
    bcftools annotate \\
        --set-id '%CHROM:%POS:%REF:%ALT' \\
        -Oz \\
        -o ${vcf.baseName}.norm.vcf.gz.tmp \\
        ${vcf.baseName}.norm.vcf.gz

    mv ${vcf.baseName}.norm.vcf.gz.tmp ${vcf.baseName}.norm.vcf.gz
    tabix -p vcf ${vcf.baseName}.norm.vcf.gz
    """
}

process MATCH_VARIANTS {

    tag "${vcf.baseName}"
    publishDir "${params.outdir}/matched", mode: 'copy'

    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_0'

    input:
    path vcf
    path tbi
    path gwas_tsv

    output:
    path "matched_variants.tsv", emit: tsv

    script:
    """
    # Extract variants from VCF into a simple TSV: CHROM, POS, REF, ALT, genotype
    bcftools query \\
        -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' \\
        ${vcf} > vcf_variants.tsv

    # Write matched-variants header
    echo -e "CHR\\tPOS\\tREF\\tALT\\tGT\\tSNP_RISK_ALLELE\\tOR_BETA\\tTRAIT" > matched_variants.tsv

    # Load VCF variants into awk, keyed by CHR:POS, then match against GWAS rows
    awk -F'\\t' '
        NR==FNR {
            # First file: VCF variants
            key = \$1 ":" \$2
            vcf[key] = \$3 "\\t" \$4 "\\t" \$5
            next
        }
        FNR > 1 {
            # Second file: filtered GWAS (skip header)
            # Columns: CHR, POS, SNP_RISK_ALLELE, P_VALUE, OR_BETA, TRAIT
            key = \$1 ":" \$2
            if (key in vcf) {
                print \$1 "\\t" \$2 "\\t" vcf[key] "\\t" \$3 "\\t" \$5 "\\t" \$6
            }
        }
    ' vcf_variants.tsv ${gwas_tsv} >> matched_variants.tsv
    """
}

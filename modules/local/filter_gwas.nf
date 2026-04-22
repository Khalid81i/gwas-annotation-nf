process FILTER_GWAS {

    tag "$gwas.baseName"
    publishDir "${params.outdir}/gwas_filtered", mode: 'copy'

    container 'quay.io/biocontainers/gawk:5.1.0'

    input:
    path gwas

    output:
    path "gwas_significant.tsv", emit: tsv

    script:
    """
    # Find column indices dynamically (GWAS Catalog columns can shift between releases)
    header=\$(head -n 1 ${gwas})

    idx_chr=\$(echo "\$header" | awk -F'\\t' '{for(i=1;i<=NF;i++) if(\$i=="CHR_ID") print i}')
    idx_pos=\$(echo "\$header" | awk -F'\\t' '{for(i=1;i<=NF;i++) if(\$i=="CHR_POS") print i}')
    idx_snp=\$(echo "\$header" | awk -F'\\t' '{for(i=1;i<=NF;i++) if(\$i=="STRONGEST SNP-RISK ALLELE") print i}')
    idx_p=\$(echo   "\$header" | awk -F'\\t' '{for(i=1;i<=NF;i++) if(\$i=="P-VALUE") print i}')
    idx_or=\$(echo  "\$header" | awk -F'\\t' '{for(i=1;i<=NF;i++) if(\$i=="OR or BETA") print i}')
    idx_trait=\$(echo "\$header" | awk -F'\\t' '{for(i=1;i<=NF;i++) if(\$i=="DISEASE/TRAIT") print i}')

    # Write new header
    echo -e "CHR\\tPOS\\tSNP_RISK_ALLELE\\tP_VALUE\\tOR_BETA\\tTRAIT" > gwas_significant.tsv

    # Filter and extract
    awk -F'\\t' -v c=\$idx_chr -v p=\$idx_pos -v s=\$idx_snp -v pv=\$idx_p -v ob=\$idx_or -v t=\$idx_trait '
        NR > 1 && \$c != "" && \$p != "" && \$pv != "" && \$pv+0 <= 5e-8 {
            print \$c "\\t" \$p "\\t" \$s "\\t" \$pv "\\t" \$ob "\\t" \$t
        }
    ' ${gwas} >> gwas_significant.tsv
    """
}

process COMPUTE_DOSAGE {

    tag "compute_dosage"
    publishDir "${params.outdir}/annotated", mode: 'copy'

    container 'quay.io/biocontainers/python:3.11'

    input:
    path matched_tsv

    output:
    path "annotated_variants.tsv", emit: tsv

    script:
    """
    python3 <<'PYEOF'

def parse_genotype(gt):
    '''Return ALT allele count (0/1/2) from a VCF genotype string, or None if unparseable.'''
    if gt in ('./.', '.', ''):
        return None
    # Handle phased (|) or unphased (/) genotypes
    alleles = gt.replace('|', '/').split('/')
    try:
        alleles = [int(a) for a in alleles]
    except ValueError:
        return None
    # Dosage = number of ALT alleles (allele index > 0)
    return sum(1 for a in alleles if a > 0)

def parse_effect_allele(snp_risk_allele):
    '''Extract the effect allele letter from a GWAS Catalog "SNP-risk allele" string like "rs123-A".'''
    if '-' not in snp_risk_allele:
        return None
    return snp_risk_allele.rsplit('-', 1)[1].strip().upper()

def align_dosage(alt_dosage, ref, alt, effect_allele):
    '''Return effect-allele dosage (count of EA), aligning direction.'''
    if alt_dosage is None or effect_allele is None:
        return 'NA'
    if effect_allele == alt.upper():
        # EA is ALT, dosage already counts EA
        return str(alt_dosage)
    elif effect_allele == ref.upper():
        # EA is REF, flip
        return str(2 - alt_dosage)
    else:
        # EA matches neither (likely strand flip or multiallelic weirdness)
        return 'NA'

with open('${matched_tsv}') as f_in, open('annotated_variants.tsv', 'w') as f_out:
    header = f_in.readline().rstrip('\\n').split('\\t')
    f_out.write('\\t'.join(header + ['EFFECT_ALLELE', 'EFFECT_ALLELE_DOSAGE']) + '\\n')

    # Expected columns: CHR, POS, REF, ALT, GT, SNP_RISK_ALLELE, OR_BETA, TRAIT
    col = {name: i for i, name in enumerate(header)}

    for line in f_in:
        fields = line.rstrip('\\n').split('\\t')
        gt    = fields[col['GT']]
        ref   = fields[col['REF']]
        alt   = fields[col['ALT']]
        sra   = fields[col['SNP_RISK_ALLELE']]

        alt_dosage    = parse_genotype(gt)
        effect_allele = parse_effect_allele(sra)
        aligned       = align_dosage(alt_dosage, ref, alt, effect_allele or '')

        out_fields = fields + [effect_allele or 'NA', aligned]
        f_out.write('\\t'.join(out_fields) + '\\n')

PYEOF
    """
}

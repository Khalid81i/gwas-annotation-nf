process SUMMARISE_TRAITS {

    tag "summarise_traits"
    publishDir "${params.outdir}/summary", mode: 'copy'

    container 'quay.io/biocontainers/python:3.11'

    input:
    path annotated_tsv

    output:
    path "trait_summary.tsv", emit: tsv

    script:
    """
    python3 <<'PYEOF'
from collections import defaultdict

traits = defaultdict(lambda: {'n_variants': 0, 'total_dosage': 0, 'weighted_beta_sum': 0.0, 'weight_sum': 0})

with open('${annotated_tsv}') as f:
    header = f.readline().rstrip('\\n').split('\\t')
    col = {name: i for i, name in enumerate(header)}

    for line in f:
        fields = line.rstrip('\\n').split('\\t')
        if not fields or fields[0] == '':
            continue
        trait   = fields[col['TRAIT']]
        dosage  = fields[col['EFFECT_ALLELE_DOSAGE']]
        or_beta = fields[col['OR_BETA']]

        if dosage == 'NA' or dosage == '':
            continue
        try:
            dosage = int(dosage)
            beta   = float(or_beta) if or_beta not in ('', 'NA') else None
        except ValueError:
            continue

        traits[trait]['n_variants']   += 1
        traits[trait]['total_dosage'] += dosage
        if beta is not None:
            traits[trait]['weighted_beta_sum'] += dosage * beta
            traits[trait]['weight_sum']        += dosage

with open('trait_summary.tsv', 'w') as f_out:
    f_out.write('TRAIT\\tN_VARIANTS\\tTOTAL_DOSAGE\\tDOSAGE_WEIGHTED_MEAN_BETA\\n')
    for trait in sorted(traits.keys()):
        t = traits[trait]
        mean_beta = (t['weighted_beta_sum'] / t['weight_sum']) if t['weight_sum'] > 0 else 'NA'
        if mean_beta != 'NA':
            mean_beta = f"{mean_beta:.4f}"
        f_out.write(f"{trait}\\t{t['n_variants']}\\t{t['total_dosage']}\\t{mean_beta}\\n")

PYEOF
    """
}

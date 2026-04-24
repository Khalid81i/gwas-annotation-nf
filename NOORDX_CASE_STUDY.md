# Case study: clinical bioinformatics at NoorDX

*A technical writeup of the variant annotation, ancestry, and pharmacogenomics work I did during a six-week bioinformatics internship at NoorDX (Jeddah, Saudi Arabia), August 2025, as a final-year Electrical and Electronic Engineering student at the University of Greenwich.*

---

## Context

NoorDX is a diagnostics and discovery laboratory based in Saudi Arabia. I joined for a six-week intensive placement with no prior professional bioinformatics background — I came in with solid programming and Linux from my Electrical and Electronic Engineering degree, and had to pick up the domain knowledge in parallel with doing the work.

The brief was open. I was given access to genotype data and asked to make something useful out of it: first to produce interpretable outputs from a single sample's VCF, and then to scale the approach into something the lab could actually use across multiple samples. Over six weeks that grew into two internally-named products — **AncestryG81** (a population genetics layer) and **PharmaGB** (a pharmacogenomics layer) — both built on the same underlying variant-annotation foundation.

This writeup covers how that foundation was built, the engineering problems that came up along the way, and what I'd change if I were starting again now.

---

## The core problem

A sample VCF file contains thousands to millions of variants. The clinically or biologically interesting information about those variants lives elsewhere — in public resources like the GWAS Catalog, PharmGKB, and 1000 Genomes reference panels. The job of a variant annotation pipeline is to stitch these together: take a VCF, pull out each variant, look up what's known about it, align directions correctly, and hand back a structured, interpretable table.

On paper this sounds trivial. In practice it's where most of the time goes, because real data doesn't behave. Chromosome naming conventions differ between files (`chr1` vs `1`), reference builds don't match, multiallelic sites complicate position lookup, effect alleles don't always match the VCF's ALT allele, rsIDs are unreliable, and whichever public resource you're pulling from will inevitably change its column layout between releases.

## What I built

The pipeline was written in Bash, using `bcftools`, `awk`, `tabix`, and `bgzip` as the core tool stack, with PLINK 2.0 handling the population genetics side and ADMIXTURE for ancestry fractions.

### Variant annotation and GWAS interpretation

**Normalisation.** The incoming VCF was first normalised with `bcftools norm -m -any` to split multiallelic sites into single-allele lines. Variant IDs were then reassigned to a consistent `CHROM:POS:REF:ALT` format using `bcftools annotate --set-id`, because relying on rsIDs proved unreliable — they were sometimes missing, sometimes outdated, and sometimes didn't match between files. Position-and-allele identity is the most robust key.

**GWAS Catalog parsing.** The GWAS Catalog is distributed as a wide TSV with a shifting column layout between releases. Rather than hard-coding column indices, the parser reads the header row and looks up the positions of the columns it needs by name — `CHR_ID`, `CHR_POS`, `STRONGEST SNP-RISK ALLELE`, `P-VALUE`, `OR or BETA`, `DISEASE/TRAIT`. Rows were filtered to genome-wide significance (*p* ≤ 5 × 10⁻⁸) and projected to a compact working table.

**Variant matching.** Sample variants were intersected with significant GWAS hits on the `CHROM:POS` key using an `awk` hash join, pulling across the trait, effect size, and effect allele for every match.

**Effect-allele-aligned dosage computation.** This is where the work actually requires some care. The VCF genotype (`0/0`, `0/1`, `1/1`) is relative to the REF/ALT pair. The GWAS effect size is relative to the *risk allele*, which may be either REF or ALT. If you naively report ALT dosage, you've got the direction of effect wrong for roughly half of all variants — an error that silently inverts downstream interpretation. The dosage pass parses the effect allele from the GWAS "strongest SNP-risk allele" string, compares it to both REF and ALT, and either keeps or flips the dosage so that the reported value always counts copies of the effect allele. Edge cases (effect allele matching neither REF nor ALT, indicating a possible strand flip or multiallelic issue) were flagged as `NA` rather than silently zeroed.

The output of this stage was a per-variant annotated table suitable for Excel-ready reporting and downstream polygenic-style aggregation, covering traits from metabolic (obesity, diabetes) and cardiovascular categories through to cognition, mood, sleep, and autoimmune risk.

### AncestryG81 — population genetics layer

On top of the annotation foundation, AncestryG81 placed a sample against global reference populations from the **1000 Genomes Project**. The workflow:

- Convert the sample VCF to PLINK binary format (`.bed`/`.bim`/`.fam`) with PLINK 2.0.
- Filter SNPs, prune for linkage disequilibrium, and merge with the 1000 Genomes reference panel. This step was where most of the debugging time went — chromosome naming mismatches between the sample and reference had to be resolved, and `--allow-extra-chr` handling was necessary for non-standard contigs that would otherwise break the merge.
- Run PCA with PLINK 2.0 on the merged panel and project the sample's coordinates onto the global population space.
- Generate 2D and 3D PCA plots with custom colour maps for 10+ reference populations (YRI, CEU, GBR, IBS, TSI, CHB, JPT, GIH, PUR, GWD and others), clustering by continental group.
- Run **ADMIXTURE** with supervised and unsupervised K-population ancestry modelling to produce ancestry fractions for the sample.
- Overlay ancestry-informative markers onto the GWAS annotation output, linking population-level variation back to phenotype associations.

The end result was a personal-ancestry report: global population coordinates, closest-population matching, regional affinity, and ancestry-weighted trait context.

### PharmaGB — pharmacogenomics layer

Parallel to AncestryG81, PharmaGB extended the annotation pipeline into drug-response interpretation. The genotype data was interpreted against a curated set of pharmacogenomically-relevant genes:

- **CYP2D6** — antidepressants, opioids, beta-blockers
- **CYP2C19** — clopidogrel, PPIs, SSRIs
- **CYP2C9** and **VKORC1** — warfarin
- **SLCO1B1** — statin myopathy risk
- **DPYD** — fluorouracil toxicity
- **TPMT** / **NUDT15** — thiopurines
- **COMT** / **DRD2** — neuropsychiatric response patterns

A secondary nutrition and supplement layer covered MTHFR/folate, vitamin D response, omega-3, caffeine metabolism, and lactose tolerance.

The output was a structured interpretation report flagging likely dose-adjustment candidates, elevated side-effect risks, metaboliser status, and alternative-medication suggestions — organised as pharma hit tables in CSV and Excel formats. Separate filters were built for the most common drug classes (statins, warfarin, clopidogrel, antidepressants, metformin, opioids).

### Reusable pipeline scripts

The work was packaged into reusable Bash pipelines:

- **`gwas_effect_pipeline.sh`** — variant annotation and dosage
- **`gwas_trait_pharm_pipeline.sh`** — combined trait and pharmacogenomics reporting

Both supported multi-sample processing across the different VCFs I worked with (KHLD, KhT, Ash, FTH, 716-1, 723-1, 529-3 and others) with fallback logic for `chr` vs no-`chr` chromosome naming, primary-chromosome fallback, automatic normalisation, and graceful handling of missing reference contigs. A run on a new sample required no manual intervention.

## Engineering problems that actually took time

Six things consumed disproportionate amounts of the six weeks, and are worth naming because they're the parts you only learn by doing:

- **Chromosome naming.** Different sources use `chr1` or `1`. A single mismatch silently drops every match. The pipeline had to detect and normalise both directions, with fallback logic.
- **Reference build mismatches.** Coordinates in hg19 versus hg38 aren't interchangeable. Getting this wrong is catastrophic and silent. Every input file's build had to be verified, not assumed.
- **Missing or alternate contigs.** VCFs occasionally reference decoy or alt contigs (`KI270706.1` and friends) that aren't in the GWAS Catalog or the 1000 Genomes reference. These had to be handled gracefully rather than crashing the merge or join.
- **Effect-allele direction.** The dosage alignment described above. Getting this wrong isn't a software bug — the pipeline runs fine and produces numbers — it's a *scientific* bug that silently inverts half your downstream signal.
- **Empty intersections.** Early versions of the pipeline were correct but produced empty output because of subtle key mismatches. Debugging an empty output is much harder than debugging a wrong output, because there's nothing to inspect.
- **Dependency fragility.** Early designs reached out to dbSNP for rsID resolution; the dependency turned out to be unreliable enough that I refactored to remove it entirely and key on position instead. Broken GWAS Catalog downloads, openpyxl edge cases, and similar ecosystem failures each consumed an afternoon at some point.

None of these are in a textbook. They're the boring, specific, real things that consume most of the wall-clock time in clinical pipeline work.

## What I delivered

By the end of the placement:

- A modular Bash pipeline (`gwas_effect_pipeline.sh`, `gwas_trait_pharm_pipeline.sh`) that took a sample VCF and produced an annotated variant table aligned to the GWAS Catalog, with correct effect-allele dosage — run across seven or more real samples.
- AncestryG81: a population genetics and ancestry reporting system built on PLINK 2.0 and ADMIXTURE, with 2D and 3D PCA visualisations and custom population colour maps.
- PharmaGB: a pharmacogenomics interpretation layer covering a curated panel of drug-response genes and a secondary nutrition layer.
- A custom Excel-format genomics knowledgebase consolidating traits, genes, genotypes, rsIDs, and confidence levels.

## What I'd do differently

Two things, both learned in the six months since finishing:

**1. Build it in Nextflow from the start.** The Bash pipeline was fine for the scale and scope it operated at, but it doesn't travel. Another person can't easily run it on their machine. It doesn't resume from a failed step. There's no built-in execution report. All the things a mature workflow manager gives you for free. I've since re-implemented the GWAS annotation core of this work as an open-source Nextflow DSL2 pipeline, [`gwas-annotation-nf`](https://github.com/Khalid81i/gwas-annotation-nf), which runs every step inside pinned Docker containers and is reproducible on any machine with Nextflow and Docker installed. That's the version I'd hand to a clinical team today.

**2. Use established functional annotation tools rather than rolling custom keyword filters.** For the pharmacogenomics layer especially, the gold-standard resources are PharmGKB (with its curated drug-gene-variant annotations and CPIC guidelines) and VEP with pharmacogenomics plugins. A production pipeline should be built on those, not on keyword matching against trait strings. The custom version worked for the scope we were operating at, but a proper clinical tool would integrate PharmGKB's annotation files directly and apply CPIC-compliant star-allele nomenclature (e.g. CYP2D6\*4) rather than single-SNP keyword filtering.

## What this experience gave me

Six weeks is short, but the density was high. I came out of it able to read a VCF fluently, design and debug a multi-stage Bash pipeline, run a full PLINK + ADMIXTURE ancestry workflow end-to-end, reason about effect-allele alignment without getting it wrong, and — importantly — recognise all the places where a naive approach silently produces wrong answers rather than obviously wrong answers. That last one is the one that matters most in clinical bioinformatics, and it's not something you can pick up from a tutorial.

If you want the cleaned-up, open-source, containerised version of the core annotation pipeline, it's here:
**[github.com/Khalid81i/gwas-annotation-nf](https://github.com/Khalid81i/gwas-annotation-nf)**

— Khalid Dallol, April 2026

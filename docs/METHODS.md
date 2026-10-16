# Methods and scope, version 0.1

## Intended use

This is a research annotation and assay-planning workflow for human chromatin-regulator missense SNVs. It implements the computational entry point to a functional assessment framework, not the editing, sequencing or phenotyping experiments themselves. It does not analyze CUT&Tag FASTQ files or reproduce the study's experimental effect sizes. Functional outcomes from the motivating KMT2C study are not training labels and are not distributed here.

The 14-gene starter panel provides examples across several chromatin functions. Membership is not a claim that every gene causes neurodevelopmental disease. The default neurodevelopmental context adds a conditional neural-phenotyping suggestion; `--context general` removes that suggestion. Gene-to-assay mappings were manually drafted from reviewed human UniProt entries, checked on 2026-10-16, and are starting points rather than exhaustive protocols.

## Input and allele resolution

1. Read uncompressed VCF, CSV/TSV, or one HGVS per text line. Split VCF ALT alleles. Keep invalid rows with an error; retain duplicates as separate input records with unique row IDs.
2. Require GRCh38. Genomic positions are one-based, with forward-strand REF/ALT. A conflicting VCF reference declaration or per-row assembly is rejected. The command assumes GRCh38 when no header declares an assembly; it cannot infer a genome build from coordinates alone.
3. Accept single-base coding substitutions on explicitly versioned `NM_` or `ENST` transcripts. Protein-only HGVS, intronic HGVS, indels and structural variants are outside this release.
4. Resolve each query through Ensembl VEP HGVS. Require one genomic resolution, GRCh38, one position, two unambiguous single-base alleles and a known strand. Complement alleles when the response is on the negative strand. For genomic input, the resolved allele must match the submitted allele.
5. Select the explicit gene if supplied; otherwise prefer genes present in the chosen panel. An explicit transcript must match the versioned VEP `hgvsc` prefix or `mane_select` exactly. Without an explicit transcript, require one gene and one MANE Select transcript, or one canonical transcript if MANE is absent. Ambiguity yields a resolution error rather than an arbitrary selection.
6. Require a missense consequence on that selected transcript and membership in the chosen panel. Returning an `annotated` status means allele/transcript resolution and missense scope succeeded; it does not promise complete scores or a usable second provider. Review `warnings`, `myvariant_status` and score coverage as well.

Variant identity is checked against VEP's mapping and MyVariant's returned identifiers and available per-source alleles/coordinates. There is no separate local reference FASTA comparison. A reference mismatch, obsolete transcript or unresolvable representation should be corrected before interpreting output.

## Evidence extraction

VEP is queried with `AlphaMissense`, `CADD`, `REVEL`, `canonical`, `mane`, `hgvs`, `protein` and `domains` enabled. Transcript-level REVEL, AlphaMissense, SIFT, PolyPhen and CADD fields are read only from the selected consequence. Missing fields remain null. Model prediction labels, including SIFT low-confidence labels when present, are retained in JSON. A requested field is not a guarantee that the service has a value.

MyVariant is queried using the normalized genomic HGVS with `assembly=hg38`. Its allele identifier must match exactly. The workflow reads ClinVar RCV significance, condition, review-status and accession fields; complete assertions remain in the raw response. Different condition-level records can differ without representing the same clinical question. The summary does not automatically resolve them.

The population fields are the overall `gnomad_exome.af.af` and `gnomad_genome.af.af` values, if returned. `gnomad_max_observed_af` is their maximum, **not an ancestry-specific maximum**. Missingness is kept distinct from zero. These services can lag the current gnomAD browser; check source versions in provenance and confirm important population claims directly.

Three genomic-site conservation values are read from dbNSFP via MyVariant: phyloP 100-way vertebrate, GERP++ RS and phastCons 100-way vertebrate. MyVariant transcript-specific score arrays are deliberately not flattened: array lengths and transcript mappings need not align. Ambiguous array-valued conservation fields are withheld rather than averaged.

## Research flags and experimental suggestions

REVEL ≥0.5 and observed AF <0.001 are default configurable research criteria, not estimates of disease probability. The defaults mirror candidate-screening choices in the motivating work. No new classifier is trained, no overall pathogenicity score is calculated and no clinical category is assigned. The AlphaMissense label is a predictor category, distinct from ClinVar clinical assertions.

Rules identify missing population data, AF above a user threshold, REVEL above/below a user threshold, potentially conflicting clinical assertions, and disagreement between the REVEL threshold and the published AlphaMissense likely-benign/likely-pathogenic categories (<0.34 or >0.564 in the current rule). This disagreement flag is a prompt for review, not a statistical independence claim.

Each successfully annotated row receives suggestions for genotype verification, allele-specific modeling, protein abundance/localization and assay controls. The gene panel supplies the primary molecular readout and source. Additional KMT2C/KMT2D suggestions connect H3K4me1 to regulatory chromatin and transcription. Neural phenotyping is conditional on gene-disease relevance. Suggestions do not infer loss versus gain of function, effect size, penetrance or causality from scores alone.

Histone-mark suggestions include total-H3 normalization and documentation of spike-in calibration. Peak overlap with an enhancer annotation is descriptive evidence; TSS distance, promoter overlap, H3K27ac/accessibility and transcription can help contextualize it. A single mutation does not establish that every substitution in the same domain behaves similarly. Independent clones and cultures, controls and rescue/reversion are proposed as future assay-design options; technical aliquots are identified as technical replication.

## Figures and exports

Coverage includes every input allele in the denominator, including errors. The heatmap shows four raw 0–1 predictor values, with SIFT explicitly labeled in the reverse direction. CADD and conservation values remain in the table because their scales differ. Missing cells are gray. Heatmaps paginate after 35 rows. Position plots group by both gene and selected transcript, zoom to the submitted position range, and do not represent a full-length protein or its domain boundaries.

HTML content is escaped. CSV string cells that could be interpreted as formulas are prefixed with an apostrophe; negative numeric values remain numeric. JSON retains the original values. Input IDs should be non-identifying because they are included in local outputs. Genotype/sample columns in VCF are neither used nor sent to the annotation services.

## Provenance and reproducibility

Each successful response is stored with URL, UTC retrieval time and SHA-256 of canonicalized parsed JSON. This hash identifies the parsed response, not the HTTP byte stream. Provenance also records input-file SHA-256, workflow version, parameters, source version metadata when exposed, status counts, request failures and hashes of output files. The gene panel is snapshotted.

The online cache freshness window is seven days. Offline mode allows stale entries and labels them in request metadata; it performs no network calls. HTTP requests have a 30-second timeout, bounded retries and rate spacing. A run's `raw/` directory is a portable cache for replaying those exact requests. Metadata retrieval failures do not fabricate version identifiers or block otherwise usable variant annotation. Plugin dataset versions that the provider does not expose are recorded as unavailable rather than inferred from an Ensembl release.

The CLI returns 0 when all rows reach `annotated` (possibly with incomplete sources), 2 when a report is created with one or more failed/out-of-scope rows, and 1 for a fatal input/output error. CI tests use synthetic data and do not prove that live upstream schemas or every supported gene/allele are stable. Live smoke validation should be reported separately from unit tests.

## Current boundaries

- Human GRCh38 missense SNVs only; no liftover, genotype/phase inference, segregation analysis or inheritance model.
- No automated ACMG/AMP classification, PS3/BS3 assignment, clinical diagnosis or calibrated disease-probability estimate.
- No coverage-aware gnomAD absence inference; no ancestry-specific popmax calculation.
- No automatic literature-mining claim that a variant has been experimentally studied.
- No model that predicts which specific cell phenotype will change.
- No claim that this implementation is benchmarked across all chromatin regulators. Broader validation requires independent variant sets and prespecified performance measures.

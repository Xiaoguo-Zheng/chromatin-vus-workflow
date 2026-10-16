# Chromatin VUS Workflow

A reproducible path from a **variant list** to **traceable database evidence, figures and an experimental follow-up plan** for human chromatin regulators.

[中文说明](README.zh-CN.md) · [Methods](docs/METHODS.md) · [Data sources and terms](docs/DATA_SOURCES.md)

The workflow connects variant selection to endogenous editing, chromatin measurements and cellular phenotyping. It was motivated by a KMT2C functional-assessment study: substitutions in the same domain can have different effects, so candidate selection should lead to allele-specific experiments. This software organizes evidence; it has not been clinically validated and does not assign ACMG classifications.

## What it does

- Accepts **human GRCh38/hg38 VCF**, genomic CSV/TSV, or versioned `NM_`/`ENST` coding HGVS.
- Resolves alleles and selects the requested transcript, otherwise a unique MANE Select or canonical transcript.
- Queries **REVEL, AlphaMissense, SIFT, PolyPhen, CADD**, three conservation measures, **ClinVar** assertions and available **gnomAD** frequencies.
- Produces searchable HTML, CSV/JSON tables, PNG/SVG figures, an assay plan and a provenance record with response hashes.
- Uses an editable **14-gene starter panel**, spanning histone and DNA modification, methylated-DNA binding and chromatin remodeling. It is not an exhaustive chromatin or NDD gene catalogue.

Version 0.1 supports **missense SNVs only**. Indels, CNVs, splice variants and protein-only inputs are reported as unsupported, rather than passed to missense predictors. Upstream services may return missing scores or errors; every input row remains in the output.

## Install and run

Python 3.10 or newer:

```bash
git clone https://github.com/Xiaoguo-Zheng/chromatin-vus-workflow.git
cd chromatin-vus-workflow
python -m venv .venv
```

Activate with `.venv\Scripts\Activate.ps1` on Windows PowerShell, or `source .venv/bin/activate` on macOS/Linux, then:

```bash
python -m pip install -e .
chromatin-vus demo --out outputs/demo
chromatin-vus annotate examples/kmt2c_hgvs.csv --out outputs/kmt2c
```

Open `outputs/kmt2c/report.html` in a browser. The demo uses **invented values** and needs no network. The KMT2C example queries public services and uses no unpublished laboratory results. A typical short list takes a few minutes, depending on service availability.

Optional upload interface:

```bash
python -m pip install -e ".[ui]"
streamlit run app.py --server.address 127.0.0.1
```

## Input formats

Genomic coordinates are **one-based, GRCh38, genomic forward strand**:

```csv
input_id,chrom,pos,ref,alt,gene,transcript,assembly
V4834A,7,152139219,A,G,KMT2C,NM_170606.3,GRCh38
```

Or supply versioned coding HGVS:

```csv
input_id,gene,hgvs
V4834A,KMT2C,NM_170606.3:c.14501T>C
```

The V4834A example illustrates negative-strand conversion: transcript `T>C` corresponds to genomic `A>G`. A plain text file can contain one HGVS per line. VCF multi-ALT records are split into separate alleles; sample genotype fields are not interpreted. VCF FILTER is not used for quality selection: prefilter variants before annotation. Use `.vcf`, not compressed VCF, in this release. No genome-build liftover is performed.

## Outputs

| File | Contents |
| --- | --- |
| `report.html` | Filterable evidence table, plots and rule-based follow-up suggestions |
| `annotations.csv`, `annotations.json` | Scores, selected transcript, consequences, domain IDs, clinical/population evidence, missingness and errors |
| `recommendations.csv`, `recommendations.json` | Proposed experiments, triggering rationale and sources |
| `figures/*.png`, `figures/*.svg` | Coverage, raw predictor heatmaps and transcript-specific position plots |
| `provenance.json` | Parameters, software/data versions when exposed, retrieval times and hashes |
| `raw/*.json` | Exact parsed provider responses with integrity metadata |
| `panel_snapshot.json` | Gene panel used for this run |

Each run requires an empty output directory. Keep the complete output directory to preserve the audit trail. Output, cache and environment folders are excluded from Git by default.

## Reuse and customize

```bash
# Replay the cached responses without network access, including stale entries.
chromatin-vus annotate examples/kmt2c_hgvs.csv --offline --out outputs/replay

# Fetch new annotations; preserve the earlier run for comparison.
chromatin-vus annotate examples/kmt2c_hgvs.csv --refresh --out outputs/refreshed

# Replace the starter panel with a reviewed panel using the same columns.
chromatin-vus annotate variants.csv --panel my_panel.csv --context general --out outputs/custom
```

The cache defaults to `.cache/chromatin-vus` with a seven-day online freshness window. A run's `raw/` directory can itself be used as `--cache` for an offline replay. Default capacity is 200 alleles per run; `--limit` can explicitly increase it, but service rate limits still apply.

Default REVEL `0.5` and observed AF `0.001` thresholds are **editable research filters**, reflecting the motivating study's selection strategy. They are not validated clinical cutoffs. Scores are retained separately; the workflow does not average them or count correlated models as independent votes. SIFT runs in the opposite direction to the other displayed 0–1 predictors. Missing frequency is not zero frequency. AlphaMissense's model label is separate from ClinVar's clinical assertions.

## How follow-up plans work

Rules connect the gene's known function to a candidate assay, then add checks for missing population data, conflicting ClinVar assertions and predictor discordance. KMT2C/KMT2D prompts include H3K4me1, chromatin distribution, transcription, genotype verification and a relevant cellular phenotype. Other gene families receive different readouts. The plan does not predict the direction or size of the experimental effect. Domain membership is returned as annotation, not treated as proof of a shared mechanism.

## Reproducibility and contribution

```bash
python -m pip install -e ".[test]"
python -m pytest -q
```

Tests use synthetic fixtures; CI needs no database access. Contributions should include a source for gene-to-assay mappings, tests for parsing/extraction changes and explicit handling of missing values. See [methods](docs/METHODS.md) for current boundaries.

Online mode sends variant HGVS queries to Ensembl and genomic allele queries to MyVariant.info. Use inputs you are authorized to submit to those services. This code is MIT-licensed; upstream data and scores have separate terms, including noncommercial restrictions for some predictors. See [data sources](docs/DATA_SOURCES.md) before redistributing annotations or using the workflow commercially.

For a publication, cite the exact release or commit, its settings and accessed source versions, plus the original resources listed in [data sources](docs/DATA_SOURCES.md). This project does not yet have a software DOI.

"""Editable assay suggestions, not a trained pathogenicity model."""
import csv
from importlib.resources import files
from .input import InputError

FUNCTIONAL_EVIDENCE = "https://doi.org/10.1186/s13073-019-0690-2"


def load_panel(path=None):
    resource = path or files("chromatin_vus").joinpath("data/panel.csv")
    with open(resource, encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    required = {"gene", "family", "readout", "experiments", "source_url"}
    if not rows or not required.issubset(rows[0]):
        raise InputError("Panel needs gene,family,readout,experiments,source_url columns.")
    panel = {}
    for row in rows:
        gene = row["gene"].strip().upper()
        if not gene or gene in panel or any(not row[k].strip() for k in required):
            raise InputError("Panel has duplicate genes or missing required values.")
        panel[gene] = {k: row[k].strip() for k in required}
    return panel


def recommend(row, panel, revel_threshold=0.5, af_threshold=0.001, context="neurodevelopment"):
    suggestions, flags = [], []

    def add(rule, rationale, action, source=FUNCTIONAL_EVIDENCE):
        suggestions.append({"input_id": row["input_id"], "gene": row.get("gene", ""), "rule": rule,
                            "rationale": rationale, "suggestion": action, "source": source})

    if row["status"] != "annotated":
        add("resolve_input", row.get("error", row["status"]), "Resolve the allele, transcript or scope issue before choosing a functional assay.", "https://rest.ensembl.org/documentation/info/vep_hgvs_get")
        return flags, suggestions
    af = row.get("gnomad_max_observed_af")
    revel = row.get("revel")
    if af is None:
        flags.append("population_data_missing")
        add("population_review", "No usable population AF was returned.", "Check gnomAD directly, including coverage and ancestry-specific counts; missing data cannot establish rarity.", "https://gnomad.broadinstitute.org/")
    elif af >= af_threshold:
        flags.append("above_configured_AF_threshold")
        add("frequency_review", f"Observed AF {af:g} is at least the user threshold {af_threshold:g}.", "Review disease prevalence, inheritance, penetrance and the relevant population before prioritization.", "https://gnomad.broadinstitute.org/")
    else:
        flags.append("below_configured_AF_threshold")
    if revel is None:
        flags.append("REVEL_missing")
    elif revel >= revel_threshold:
        flags.append("REVEL_above_configured_threshold")
    else:
        flags.append("REVEL_below_configured_threshold")
    labels = " | ".join(row.get("clinvar_significance", [])).lower()
    if "conflict" in labels or ("benign" in labels and "pathogenic" in labels):
        flags.append("ClinVar_review_needed")
        add("clinical_conflict", "ClinVar records contain conflicting classifications or different condition-level assertions.", "Inspect the original submissions, conditions, review status and dates before interpreting the discrepancy.", "https://www.ncbi.nlm.nih.gov/clinvar/")
    am = row.get("alphamissense")
    if revel is not None and am is not None and ((revel >= revel_threshold and am < 0.34) or (revel < revel_threshold and am > 0.564)):
        flags.append("predictor_discordance")
        add("prediction_review", "REVEL threshold and AlphaMissense model category disagree.", "Check transcript matching and model scope. Use assays to test a hypothesis; do not count correlated predictors as independent evidence.", "https://doi.org/10.1126/science.adg7492")
    add("allele_model", "The workflow resolves a missense substitution, which may alter several protein functions.", "Confirm endogenous genotype and allele composition; compare matched wild type with the edited allele. Measure protein abundance and localization. A knockout alone does not establish the effect of this substitution.")
    info = panel[row["gene"]]
    add("molecular_readout", f"{row['gene']} belongs to {info['family']}; candidate readout: {info['readout']}.", info["experiments"], info["source_url"])
    add("assay_controls", "Clonal effects and assay variation can mimic variant effects.", "Use independently derived clones and independently initiated cultures, appropriate reference variants and rescue or reversion where feasible. Split aliquots measure technical repeatability. Plan sample size from assay variation.")
    if row["gene"] in {"KMT2C", "KMT2D"}:
        add("chromatin_to_function", "H3K4me1 changes motivate tests of regulatory function; peak location alone does not define enhancer activity.", "Quantify histone marks relative to total H3 and document any spike-in calibration. Profile H3K4me1 with H3K27ac or accessibility and RNA measurements. Report distance to TSS, promoter overlap and candidate distal regulatory regions separately. Test catalytic and scaffolding contributions.", "https://doi.org/10.1016/j.molcel.2024.02.030")
    if context == "neurodevelopment":
        add("cellular_context", "The user selected a neurodevelopmental research context; disease relevance must be established for the gene and allele.", "If supported by gene-disease evidence, test a defined neural differentiation model with staged lineage markers, cell viability and a prespecified phenotype. Compare molecular and cellular effects; a cellular difference alone does not assign clinical pathogenicity.")
    return flags, suggestions

"""Allele- and transcript-aware extraction. Missing values stay missing."""
from urllib.parse import quote, urlencode
import math
import re
from .http import ProviderError
from .input import InputError, genomic_hgvs

ENSEMBL = "https://rest.ensembl.org"
MYVARIANT = "https://myvariant.info/v1"
SCORES = ["revel", "alphamissense", "sift", "polyphen", "cadd_phred", "phylop100way", "gerp_rs", "phastcons100way"]


def number(value):
    if isinstance(value, bool) or isinstance(value, (dict, list)) or value in (None, "", "."):
        return None
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (ValueError, TypeError):
        return None


def nested(obj, path):
    for part in path.split("."):
        if not isinstance(obj, dict):
            return None
        obj = obj.get(part)
    return obj


def as_list(value):
    return value if isinstance(value, list) else [] if value is None else [value]


def normalize_vep_record(record):
    if record.get("assembly_name") != "GRCh38":
        raise InputError("VEP returned a different or unspecified genome assembly.")
    if record.get("start") != record.get("end"):
        raise InputError("Resolved variant is not an SNV; outside v0.1 scope.")
    alleles = record.get("allele_string", "").split("/")
    if len(alleles) != 2 or any(a not in "ACGT" or len(a) != 1 for a in alleles):
        raise InputError("The HGVS resolves to ambiguous or non-SNV alleles.")
    if record.get("strand") == -1:
        alleles = [a.translate(str.maketrans("ACGT", "TGCA")) for a in alleles]
    elif record.get("strand") != 1:
        raise InputError("VEP did not specify the genomic allele strand.")
    return genomic_hgvs(record["seq_region_name"], record["start"], *alleles)


def select_transcript(record, variant, panel):
    transcripts = [t for t in record.get("transcript_consequences", []) if t.get("gene_symbol")]
    if variant.gene:
        transcripts = [t for t in transcripts if t["gene_symbol"] == variant.gene]
        if not transcripts:
            raise InputError("The supplied gene does not match VEP consequences.")
    else:
        panel_matches = [t for t in transcripts if t["gene_symbol"] in panel]
        if panel_matches:
            transcripts = panel_matches
    if variant.transcript:
        requested = variant.transcript
        transcripts = [t for t in transcripts if t.get("mane_select") == requested or t.get("hgvsc", "").split(":")[0] == requested]
        if len(transcripts) != 1:
            raise InputError("The requested transcript version has no unique VEP match. Supply current MANE or a matching versioned ENST identifier.")
        return transcripts[0], "requested_transcript"
    if len({t["gene_symbol"] for t in transcripts}) != 1:
        raise InputError("Multiple or no genes match; supply a gene and transcript.")
    preferred = [t for t in transcripts if t.get("mane_select")]
    rule = "MANE_Select"
    if not preferred:
        preferred = [t for t in transcripts if t.get("canonical")]
        rule = "unique_canonical"
    if len(preferred) != 1:
        raise InputError("No unique MANE Select/canonical transcript; specify a transcript.")
    return preferred[0], rule


def extract_myvariant(obj, genomic):
    expected = "chr" + genomic
    if not isinstance(obj, dict) or obj.get("_id") != expected:
        raise InputError("MyVariant returned an allele identifier inconsistent with the GRCh38 query.")
    result = {"clinvar_significance": [], "clinvar_review_status": [], "clinvar_conditions": [], "clinvar_accessions": [],
              "gnomad_exome_af": None, "gnomad_genome_af": None, "gnomad_max_observed_af": None,
              "phylop100way": None, "gerp_rs": None, "phastcons100way": None, "warnings": []}
    match = re.fullmatch(r"([^:]+):g\.(\d+)([ACGT])>([ACGT])", genomic)
    chrom, pos, ref, alt = match.groups()
    for key in ["clinvar", "dbnsfp", "gnomad_exome", "gnomad_genome"]:
        source = obj.get(key)
        if not isinstance(source, dict):
            continue
        for name, expected_value in [("ref", ref), ("alt", alt), ("chrom", chrom)]:
            if name in source and str(source[name]).removeprefix("chr") != expected_value:
                raise InputError(f"{key} allele identity disagrees with the query.")
        position = nested(source, "hg38.start")
        if position is not None and int(position) != int(pos):
            raise InputError(f"{key} coordinates disagree with GRCh38 query.")
    clinvar = obj.get("clinvar") or {}
    if not isinstance(clinvar, dict):
        raise InputError("Unsupported ClinVar response schema; annotation withheld.")
    for record in as_list(clinvar.get("rcv")):
        if not isinstance(record, dict):
            continue
        for source, dest in [("clinical_significance", "clinvar_significance"), ("review_status", "clinvar_review_status"), ("accession", "clinvar_accessions")]:
            result[dest].extend(str(v) for v in as_list(record.get(source)) if v)
        for condition in as_list(record.get("conditions")):
            if isinstance(condition, dict) and condition.get("name"):
                result["clinvar_conditions"].append(condition["name"])
    for key in ["clinvar_significance", "clinvar_review_status", "clinvar_conditions", "clinvar_accessions"]:
        result[key] = sorted(set(result[key]))
    for dataset in ["gnomad_exome", "gnomad_genome"]:
        result[dataset + "_af"] = number(nested(obj, dataset + ".af.af"))
        if result[dataset + "_af"] is not None and not 0 <= result[dataset + "_af"] <= 1:
            raise InputError("gnomAD returned an AF outside [0, 1].")
    frequencies = [result[d + "_af"] for d in ["gnomad_exome", "gnomad_genome"] if result[d + "_af"] is not None]
    result["gnomad_max_observed_af"] = max(frequencies) if frequencies else None
    for key, path in {"phylop100way": "dbnsfp.phylop.100way_vertebrate.score", "gerp_rs": "dbnsfp.gerp++.rs", "phastcons100way": "dbnsfp.phastcons.100way_vertebrate.score"}.items():
        result[key] = number(nested(obj, path))
    return result


def annotate(variant, client, panel):
    result = variant.to_dict() | {"status": "input_error" if variant.error else "pending", "warnings": [], "score_sources": {}}
    result.update({k: None for k in SCORES})
    if variant.error:
        return result
    options = {"AlphaMissense": 1, "CADD": 1, "REVEL": 1, "canonical": 1, "mane": 1, "hgvs": 1, "protein": 1, "domains": 1, "content-type": "application/json"}
    url = ENSEMBL + "/vep/human/hgvs/" + quote(variant.hgvs, safe="") + "?" + urlencode(options)
    try:
        records = client.get(url)
        if not isinstance(records, list) or len(records) != 1:
            raise InputError("HGVS has no unique genomic resolution; provide an exact genomic allele.")
        record = records[0]
        genomic = normalize_vep_record(record)
        if ":g." in variant.hgvs and genomic != variant.hgvs:
            raise InputError("VEP allele does not match the supplied genomic allele.")
        t, rule = select_transcript(record, variant, panel)
        result.update({"genomic_hgvs": genomic, "gene": t["gene_symbol"], "transcript_id": t.get("transcript_id"),
                       "mane_select": t.get("mane_select"), "selection_rule": rule, "hgvsc": t.get("hgvsc"), "hgvsp": t.get("hgvsp"),
                       "protein_position": t.get("protein_start"), "consequences": t.get("consequence_terms", []),
                       "domains": sorted({str(d.get("db")) + ":" + str(d.get("name")) for d in t.get("domains", [])})})
        if result["gene"] not in panel:
            result.update(status="out_of_scope", error="Gene is outside the selected epigenetic-regulator panel. Provide an explicit custom panel to expand it.")
            return result
        if "missense_variant" not in result["consequences"]:
            result.update(status="out_of_scope", error="This is not a missense consequence on the selected transcript; do not apply missense predictors to this allele.")
            return result
        for dest, path in {"revel": "revel", "alphamissense": "alphamissense.am_pathogenicity", "sift": "sift_score", "polyphen": "polyphen_score", "cadd_phred": "cadd_phred"}.items():
            result[dest] = number(nested(t, path))
            if dest != "cadd_phred" and result[dest] is not None and not 0 <= result[dest] <= 1:
                result[dest] = None
                result["warnings"].append(f"{dest}: provider value outside [0, 1]; withheld")
            if result[dest] is not None:
                result["score_sources"][dest] = {"provider": "Ensembl VEP", "transcript": t.get("transcript_id"), "field": path}
        result["alphamissense_class"] = nested(t, "alphamissense.am_class")
        result["predictor_labels"] = {k: t.get(k) for k in ["sift_prediction", "polyphen_prediction"]}
        for predictor, label in result["predictor_labels"].items():
            if "low_confidence" in str(label).lower():
                result["warnings"].append(f"{predictor}: {label}; review the predictor's confidence before interpretation.")
        result["warnings"].extend(str(w) for w in as_list(record.get("warnings")))
        result["status"] = "annotated"
    except ProviderError as exc:
        result.update(status="provider_error", error=f"Ensembl VEP: {exc}")
        return result
    except (InputError, KeyError, TypeError, ValueError, AttributeError) as exc:
        result.update(status="resolution_error", error=str(exc))
        return result
    url = MYVARIANT + "/variant/" + quote("chr" + genomic, safe="") + "?" + urlencode({"assembly": "hg38", "fields": "dbnsfp,clinvar,gnomad_exome,gnomad_genome,hg38"})
    try:
        extra = extract_myvariant(client.get(url), genomic)
        result["warnings"].extend(extra.pop("warnings"))
        result.update(extra)
        result["myvariant_status"] = "matched"
        for key in ["phylop100way", "gerp_rs", "phastcons100way"]:
            if result[key] is not None:
                result["score_sources"][key] = {"provider": "MyVariant.info / dbNSFP", "scope": "genomic-site conservation"}
    except (ProviderError, InputError, TypeError, ValueError, AttributeError) as exc:
        result["myvariant_status"] = "not_found" if str(exc) == "HTTP 404" else "unavailable"
        result["warnings"].append(f"MyVariant: {exc}; clinical and population annotation may be incomplete")
    if result.get("gnomad_max_observed_af") is None:
        result["warnings"].append("Population AF unavailable; this does not mean absence from gnomAD or AF=0.")
    missing = [key for key in SCORES if result[key] is None]
    if missing:
        result["warnings"].append("Missing scores: " + ", ".join(missing))
    return result

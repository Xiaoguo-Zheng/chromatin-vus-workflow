import copy
import hashlib
import json
from pathlib import Path
import pytest

from chromatin_vus.annotation import annotate, normalize_vep_record, select_transcript, extract_myvariant
from chromatin_vus.http import Client, ProviderError, utc_now
from chromatin_vus.input import Variant, InputError, parse_hgvs, read_variants
from chromatin_vus.pipeline import demo, csv_cell, run
from chromatin_vus.rules import load_panel, recommend


@pytest.fixture
def vep():
    # Entirely synthetic fixture: identifiers/coordinates do not describe a real allele.
    return {"assembly_name": "GRCh38", "seq_region_name": "7", "start": 100, "end": 100, "allele_string": "T/C", "strand": -1,
            "transcript_consequences": [{"gene_symbol": "KMT2C", "transcript_id": "ENST00000000001", "hgvsc": "ENST00000000001.1:c.10T>C",
                                        "mane_select": "NM_000001.1", "hgvsp": "ENSP00000000001.1:p.Val4Ala", "canonical": 1,
                                        "consequence_terms": ["missense_variant"], "protein_start": 4,
                                        "revel": 0.7, "alphamissense": {"am_pathogenicity": 0.8, "am_class": "likely_pathogenic"},
                                        "sift_score": 0.01, "polyphen_score": 0.9, "cadd_phred": 20}]}


def test_negative_strand(vep):
    assert normalize_vep_record(vep) == "7:g.100A>G"


@pytest.mark.parametrize("update", [{"assembly_name": "GRCh37"}, {"end": 101}, {"strand": 0}, {"allele_string": "T/CA"}])
def test_invalid_resolution(vep, update):
    vep.update(update)
    with pytest.raises(InputError):
        normalize_vep_record(vep)


def test_exact_transcript_version(vep):
    with pytest.raises(InputError, match="version"):
        select_transcript(vep, Variant("x", transcript="NM_000001.2"), load_panel())
    assert select_transcript(vep, Variant("x", transcript="NM_000001.1"), load_panel())[1] == "requested_transcript"


def test_ambiguous_mane(vep):
    vep["transcript_consequences"].append(copy.deepcopy(vep["transcript_consequences"][0]))
    with pytest.raises(InputError, match="unique"):
        select_transcript(vep, Variant("x"), load_panel())


def test_gene_mismatch(vep):
    with pytest.raises(InputError, match="gene"):
        select_transcript(vep, Variant("x", gene="CHD8"), load_panel())


@pytest.mark.parametrize("hgvs", ["KMT2C:p.Val4834Ala", "NM_170606:c.14501T>C", "NM_170606.3:c.14501del", "7:g.100A>A"])
def test_reject_unsupported_hgvs(hgvs):
    with pytest.raises(InputError):
        parse_hgvs(hgvs)


def test_multi_alt_and_invalid_preserved(tmp_path):
    path = tmp_path / "input.vcf"
    path.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\n7\t100\tx\tA\tG,T,AG\n")
    rows = read_variants(path)
    assert len(rows) == 3 and rows[0].hgvs == "7:g.100A>G" and rows[1].hgvs == "7:g.100A>T"
    assert rows[2].error and len({r.input_id for r in rows}) == 3


def test_wrong_assembly(tmp_path):
    path = tmp_path / "input.vcf"
    path.write_text("##reference=GRCh37\n7\t100\t.\tA\tG\n")
    with pytest.raises(InputError, match="GRCh37"):
        read_variants(path)


def test_empty_input(tmp_path):
    path = tmp_path / "input.txt"
    path.write_text("")
    with pytest.raises(InputError, match="no variant"):
        read_variants(path)


def test_single_column_hgvs_header(tmp_path):
    path = tmp_path / "input.csv"
    path.write_text("hgvs\nNM_000001.1:c.10T>C\n")
    rows = read_variants(path)
    assert len(rows) == 1 and not rows[0].error


def test_malformed_provider_schema():
    row = annotate(Variant("test", "7:g.100A>G"), FakeClient("schema_changed"), load_panel())
    assert row["status"] == "resolution_error"


def test_low_confidence_visible(vep):
    vep["transcript_consequences"][0]["sift_prediction"] = "deleterious_low_confidence"
    row = annotate(Variant("test", "7:g.100A>G"), FakeClient(vep), load_panel())
    assert any("low_confidence" in w for w in row["warnings"])


def test_csv_transcript_conflict(tmp_path):
    path = tmp_path / "input.csv"
    path.write_text("hgvs,transcript\nNM_000001.1:c.10T>C,NM_000001.2\n")
    assert "conflicts" in read_variants(path)[0].error


def test_missing_not_zero_frequency():
    data = {"_id": "chr7:g.100A>G"}
    assert extract_myvariant(data, "7:g.100A>G")["gnomad_max_observed_af"] is None
    data["gnomad_exome"] = {"af": {"af": 0}}
    assert extract_myvariant(data, "7:g.100A>G")["gnomad_max_observed_af"] == 0


def test_wrong_allele():
    with pytest.raises(InputError):
        extract_myvariant({"_id": "chr7:g.100A>T"}, "7:g.100A>G")
    with pytest.raises(InputError):
        extract_myvariant({"_id": "chr7:g.100A>G", "clinvar": {"alt": "T"}}, "7:g.100A>G")


def test_conservation_only_no_transcript_array_flattening():
    data = {"_id": "chr7:g.100A>G", "dbnsfp": {"revel": {"score": [0.1, 0.9]}, "phylop": {"100way_vertebrate": {"score": 2.5}}}}
    result = extract_myvariant(data, "7:g.100A>G")
    assert "revel" not in result and result["phylop100way"] == 2.5


def test_clinvar_conflicting_records():
    data = {"_id": "chr7:g.100A>G", "clinvar": {"rcv": [{"clinical_significance": "Benign", "accession": "RCVsynthetic1"}, {"clinical_significance": "Pathogenic", "accession": "RCVsynthetic2"}]}}
    row = extract_myvariant(data, "7:g.100A>G") | {"status": "annotated", "gene": "KMT2C", "input_id": "test"}
    flags, _ = recommend(row, load_panel())
    assert "ClinVar_review_needed" in flags


class FakeClient:
    def __init__(self, vep, fail=False):
        self.vep, self.fail = vep, fail

    def get(self, url):
        if "/vep/" in url:
            return [self.vep]
        if self.fail:
            raise ProviderError("HTTP 404")
        return {"_id": "chr7:g.100A>G"}


def test_annotate_selected_transcript(vep):
    row = annotate(Variant("test", "NM_000001.1:c.10T>C", transcript="NM_000001.1"), FakeClient(vep), load_panel())
    assert row["status"] == "annotated" and row["revel"] == 0.7 and row["sift"] == 0.01
    assert row["genomic_hgvs"] == "7:g.100A>G" and row["gnomad_max_observed_af"] is None


def test_partial_provider_failure_keeps_scores(vep):
    row = annotate(Variant("test", "7:g.100A>G"), FakeClient(vep, fail=True), load_panel())
    assert row["status"] == "annotated" and row["myvariant_status"] == "not_found" and row["revel"] == 0.7
    assert any("does not mean" in w for w in row["warnings"])


def test_out_of_scope_no_assay(vep):
    vep["transcript_consequences"][0]["consequence_terms"] = ["synonymous_variant"]
    row = annotate(Variant("test", "7:g.100A>G"), FakeClient(vep), load_panel())
    assert row["status"] == "out_of_scope" and row["revel"] is None
    assert recommend(row, load_panel())[1][0]["rule"] == "resolve_input"


def test_offline_cache_missing_and_integrity(tmp_path):
    client = Client(tmp_path, offline=True)
    url = "https://example.invalid/synthetic"
    with pytest.raises(ProviderError, match="No valid"):
        client.get(url)
    payload = {"synthetic": True}
    envelope = {"url": url, "retrieved_at": utc_now(), "sha256": hashlib.sha256(json.dumps(payload, sort_keys=True, ensure_ascii=False).encode()).hexdigest(), "payload": payload}
    path = tmp_path / (hashlib.sha256(url.encode()).hexdigest() + ".json")
    path.write_text(json.dumps(envelope))
    assert client.get(url) == payload
    envelope["payload"]["synthetic"] = False
    path.write_text(json.dumps(envelope))
    with pytest.raises(ProviderError, match="integrity"):
        client.get(url)


def test_demo_report_and_safe_html(tmp_path):
    out = tmp_path / "demo"
    manifest = demo(out)
    assert manifest["synthetic_demo"] and manifest["n_annotated"] == 3
    html = (out / "report.html").read_text(encoding="utf-8")
    assert "SYNTHETIC DEMO" in html and "SIFT (lower)" in html
    assert len(list((out / "figures").glob("*.svg"))) == 3
    from chromatin_vus.report import write_html
    row = {"input_id": '<img src=x onerror="alert(1)">', "status": "input_error"}
    write_html([row], [], manifest, [], tmp_path / "xss.html")
    html = (tmp_path / "xss.html").read_text(encoding="utf-8")
    assert '<img src=x' not in html and "&lt;img" in html
    with pytest.raises(InputError, match="not empty"):
        demo(out)


@pytest.mark.parametrize("value", ["=1+1", " +SUM(A1)", "@cmd", "-formula"])
def test_csv_formula_escaped(value):
    assert csv_cell(value).startswith("'")
    assert csv_cell(-1.5) == -1.5


def test_offline_pipeline_reports_failure(tmp_path):
    path = tmp_path / "input.txt"
    path.write_text("7:g.100A>G\n")
    result = run(path, tmp_path / "out", cache=tmp_path / "empty", offline=True)
    assert result["n_annotated"] == 0 and result["status_counts"] == {"provider_error": 1}
    assert (tmp_path / "out" / "report.html").exists()

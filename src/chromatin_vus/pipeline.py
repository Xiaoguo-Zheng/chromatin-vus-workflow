from pathlib import Path
from collections import Counter
import csv
import hashlib
import json
import shutil
import sys
from importlib.metadata import version

from . import __version__
from .annotation import annotate, SCORES, ENSEMBL, MYVARIANT
from .http import Client, ProviderError, utc_now
from .input import read_variants, InputError
from .rules import load_panel, recommend
from .report import plot_results, write_html


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8")


def csv_cell(value):
    if isinstance(value, (list, dict)):
        value = json.dumps(value, ensure_ascii=False)
    if isinstance(value, str) and value.lstrip().startswith(("=", "+", "-", "@")):
        return "'" + value
    return "" if value is None else value


def write_csv(path, rows):
    columns = list(dict.fromkeys(key for row in rows for key in row))
    with Path(path).open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: csv_cell(v) for k, v in row.items()})


def prepare_output(output):
    output = Path(output)
    if output.exists() and any(output.iterdir()):
        raise InputError("Output directory is not empty; choose a new directory to preserve earlier runs.")
    output.mkdir(parents=True, exist_ok=True)
    return output


def finish(rows, panel, output, manifest):
    suggestions = []
    params = manifest["parameters"]
    for row in rows:
        flags, items = recommend(row, panel, params["revel_threshold"], params["af_threshold"], params["context"])
        row["research_flags"] = flags
        suggestions.extend(items)
    manifest["status_counts"] = dict(Counter(r["status"] for r in rows))
    manifest["runtime"] = {"python": sys.version.split()[0], "matplotlib": version("matplotlib"), "numpy": version("numpy")}
    manifest["n_input_alleles"] = len(rows)
    manifest["n_annotated"] = sum(r["status"] == "annotated" for r in rows)
    write_json(output / "annotations.json", rows)
    write_csv(output / "annotations.csv", rows)
    write_json(output / "recommendations.json", suggestions)
    write_csv(output / "recommendations.csv", suggestions)
    write_json(output / "panel_snapshot.json", panel)
    plots = plot_results(rows, output / "figures", manifest.get("synthetic_demo", False))
    write_html(rows, suggestions, manifest, plots, output / "report.html")
    manifest["output_sha256"] = {str(p.relative_to(output)).replace("\\", "/"): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(output.rglob("*")) if p.is_file() and p.name != "provenance.json"}
    write_json(output / "provenance.json", manifest)
    return manifest


def run(input_path, output, cache=".cache/chromatin-vus", panel_path=None, offline=False, refresh=False,
        revel_threshold=0.5, af_threshold=0.001, context="neurodevelopment", limit=200, progress=None):
    if not 0 <= revel_threshold <= 1 or not 0 <= af_threshold <= 1:
        raise InputError("Thresholds must be in [0, 1].")
    if offline and refresh:
        raise InputError("--offline and --refresh cannot be combined.")
    if context not in {"neurodevelopment", "general"}:
        raise InputError("Context must be neurodevelopment or general.")
    variants = read_variants(input_path, limit=limit)
    panel = load_panel(panel_path)
    output = prepare_output(output)
    client = Client(cache, offline=offline, refresh=refresh)
    manifest = {"workflow": "chromatin-vus-workflow", "version": __version__, "created_at": utc_now(), "assembly": "GRCh38",
                "input_sha256": hashlib.sha256(Path(input_path).read_bytes()).hexdigest(),
                "parameters": {"revel_threshold": revel_threshold, "af_threshold": af_threshold, "context": context, "limit": limit, "offline": offline, "refresh": refresh},
                "provider_versions": {}, "version_notes": "VEP plugin dataset versions are not exposed by these responses; do not infer them from the Ensembl release.", "provider_warnings": []}
    for name, url in {"ensembl_software": ENSEMBL + "/info/software?content-type=application/json", "ensembl_data": ENSEMBL + "/info/data?content-type=application/json", "myvariant": MYVARIANT + "/metadata"}.items():
        try:
            data = client.get(url)
            if name == "myvariant":
                data = {"build_version": data.get("build_version"), "sources": {k: v.get("version") for k, v in data.get("src", {}).items() if k in {"dbnsfp", "clinvar", "gnomad_exome", "gnomad_genome"}}}
            manifest["provider_versions"][name] = data
        except (ProviderError, AttributeError, TypeError) as exc:
            manifest["provider_warnings"].append(f"{name}: {exc}")
    rows = []
    for i, variant in enumerate(variants, 1):
        row = annotate(variant, client, panel)
        rows.append(row)
        if progress:
            progress(i, len(variants), row)
    manifest["requests"] = client.events
    raw = output / "raw"
    raw.mkdir()
    for event in client.events:
        if "sha256" not in event:
            continue
        key = hashlib.sha256(event["url"].encode()).hexdigest()
        source = client.cache / (key + ".json")
        if source.exists():
            shutil.copyfile(source, raw / source.name)
    return finish(rows, panel, output, manifest)


def demo(output):
    """Invented data to exercise output rendering without querying or bundling licensed scores."""
    output = prepare_output(output)
    rows = []
    for i, (revel, am, sift, pp, af) in enumerate([(0.82, 0.91, 0.02, 0.94, None), (0.22, 0.81, 0.2, 0.45, 0.0), (None, 0.12, None, 0.08, 0.012)], 1):
        rows.append({"input_id": f"SYNTHETIC_{i}", "gene": "KMT2C", "genomic_hgvs": "synthetic_only", "hgvsc": "synthetic_only", "hgvsp": f"synthetic_position_{100 + i * 10}",
                     "transcript_id": "SYNTHETIC_TRANSCRIPT", "protein_position": 100 + i * 10, "status": "annotated", "error": "",
                     "revel": revel, "alphamissense": am, "sift": sift, "polyphen": pp, "cadd_phred": 10 + i,
                     "phylop100way": None, "gerp_rs": None, "phastcons100way": None, "gnomad_max_observed_af": af,
                     "clinvar_significance": [], "warnings": ["Invented demonstration data; no database query performed."], "score_sources": {"all": "SYNTHETIC"}})
    manifest = {"workflow": "chromatin-vus-workflow", "version": __version__, "created_at": utc_now(), "synthetic_demo": True,
                "parameters": {"revel_threshold": 0.5, "af_threshold": 0.001, "context": "neurodevelopment"}, "requests": []}
    return finish(rows, load_panel(), output, manifest)

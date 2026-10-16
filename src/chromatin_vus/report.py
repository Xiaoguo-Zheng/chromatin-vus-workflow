"""Portable HTML and scientific plots, keeping each score on its own scale."""
from pathlib import Path
from html import escape
import math

from .annotation import SCORES

LABELS = {"revel": "REVEL (higher)", "alphamissense": "AlphaMissense (higher)",
          "sift": "SIFT (lower)", "polyphen": "PolyPhen (higher)", "cadd_phred": "CADD PHRED (higher)",
          "phylop100way": "phyloP 100-way", "gerp_rs": "GERP++ RS", "phastcons100way": "phastCons 100-way"}


def display(value):
    if value is None or value == "":
        return "NA"
    if isinstance(value, float):
        return f"{value:.4g}"
    if isinstance(value, list):
        return "; ".join(map(str, value)) or "NA"
    return str(value)


def plot_results(rows, directory, synthetic=False):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    directory = Path(directory)
    directory.mkdir(exist_ok=True, parents=True)
    paths = []
    title_prefix = "SYNTHETIC DEMO | " if synthetic else ""
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False, "svg.fonttype": "none"})

    def save(fig, name):
        fig.savefig(directory / (name + ".png"), dpi=160, bbox_inches="tight")
        fig.savefig(directory / (name + ".svg"), bbox_inches="tight")
        paths.append("figures/" + name + ".png")
        plt.close(fig)

    valid = [r for r in rows if r["status"] == "annotated"]
    fig, ax = plt.subplots(figsize=(9, 4))
    denominator = len(rows)
    counts = [sum(r.get(k) is not None and r["status"] == "annotated" for r in rows) for k in SCORES]
    bars = ax.barh([LABELS[k] for k in SCORES], counts, color="#157d87")
    ax.bar_label(bars, labels=[f"{n}/{denominator}" for n in counts], padding=5)
    ax.set(xlim=(0, max(1, denominator) * 1.18), xlabel="Input alleles with usable scores (includes failed rows in denominator)",
           title=title_prefix + "Annotation coverage")
    ax.invert_yaxis()
    fig.tight_layout()
    save(fig, "coverage")
    for page in range(math.ceil(len(valid) / 35)):
        chunk = valid[page * 35:(page + 1) * 35]
        keys = ["revel", "alphamissense", "sift", "polyphen"]
        values = np.array([[r.get(k) if r.get(k) is not None else np.nan for k in keys] for r in chunk], dtype=float)
        fig, ax = plt.subplots(figsize=(9, max(3.5, len(chunk) * 0.38 + 1.6)))
        cmap = plt.get_cmap("YlGnBu").with_extremes(bad="#e8e8e8")
        im = ax.imshow(np.ma.masked_invalid(values), cmap=cmap, vmin=0, vmax=1, aspect="auto")
        ax.set_xticks(range(4), [LABELS[k].replace(" (", "\n(") for k in keys], fontsize=9)
        ax.set_yticks(range(len(chunk)), [f"{r['input_id'][:35]} | {r.get('gene', '')} {r.get('hgvsp', '').split(':')[-1]}" for r in chunk])
        for i in range(len(chunk)):
            for j in range(4):
                v = values[i, j]
                ax.text(j, i, "NA" if np.isnan(v) else f"{v:.3f}", ha="center", va="center", color="white" if not np.isnan(v) and v > 0.55 else "#18252f", fontsize=8)
        ax.set_title(title_prefix + "Raw prediction scores | SIFT has the opposite direction", pad=18)
        fig.colorbar(im, ax=ax, label="Raw score", shrink=0.7)
        fig.tight_layout()
        save(fig, f"scores_{page + 1}")
    # Separate tracks per gene AND selected transcript; never mix isoform coordinates.
    groups = sorted({(r["gene"], r.get("transcript_id") or "unknown") for r in valid if r.get("protein_position")})
    for index, (gene, transcript) in enumerate(groups, 1):
        items = [r for r in valid if r["gene"] == gene and (r.get("transcript_id") or "unknown") == transcript and r.get("protein_position")]
        for page in range(math.ceil(len(items) / 30)):
            chunk = items[page * 30:(page + 1) * 30]
            fig, ax = plt.subplots(figsize=(9, max(3.3, len(chunk) * 0.33 + 1.5)))
            positions = [r["protein_position"] for r in chunk]
            ax.scatter(positions, range(len(chunk)), color="#c76934", s=60)
            ax.set_yticks(range(len(chunk)), [r["input_id"][:40] + " | " + r.get("hgvsp", "").split(":")[-1] for r in chunk])
            pad = max(5, (max(positions) - min(positions)) * 0.12)
            ax.set_xlim(max(0, min(positions) - pad), max(positions) + pad)
            ax.set(xlabel="Selected-transcript amino-acid position (zoomed input range; not full protein)",
                   title=title_prefix + f"{gene} / {transcript}")
            ax.grid(axis="x", alpha=0.2)
            fig.tight_layout()
            save(fig, f"positions_{index}_{page + 1}")
    return paths


def write_html(rows, suggestions, manifest, plots, path):
    def e(value):
        return escape(display(value), quote=True)
    columns = ["input_id", "gene", "genomic_hgvs", "hgvsc", "hgvsp", "status"] + SCORES + ["gnomad_max_observed_af", "clinvar_significance", "research_flags", "warnings", "error"]
    headers = "".join(f"<th>{e(LABELS.get(k, k))}</th>" for k in columns)
    body = "".join("<tr>" + "".join(f"<td>{e(row.get(k))}</td>" for k in columns) + "</tr>" for row in rows)
    cards = "".join(f"<article><b>{e(s['input_id'])} · {e(s['rule'])}</b><p>{e(s['rationale'])}</p><p>{e(s['suggestion'])}</p><small>Source: {e(s['source'])}</small></article>" for s in suggestions)
    figures = "".join(f'<figure><img src="{escape(p, quote=True)}" alt="{escape(Path(p).stem)}"></figure>' for p in plots)
    n_ok = sum(r["status"] == "annotated" for r in rows)
    demo = '<p class="notice">SYNTHETIC DEMO — invented scores and variants; not biological observations.</p>' if manifest.get("synthetic_demo") else ""
    html = """<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Chromatin VUS report</title><style>
body{font-family:system-ui,sans-serif;margin:0;color:#20313c;background:#f5f7f8}main{max-width:1300px;margin:auto;padding:36px 24px}h1{font-size:32px;margin-bottom:8px}h2{margin-top:32px}p{line-height:1.6}.meta{color:#536773}.notice{background:#fff0cc;padding:14px;border-left:4px solid #a26300}.tablewrap{overflow:auto;max-height:600px;background:white}table{border-collapse:collapse;font-size:13px;width:100%}td,th{padding:10px;border-bottom:1px solid #dce5e8;text-align:left;min-width:100px;max-width:360px;overflow-wrap:anywhere}th{position:sticky;top:0;background:#163e50;color:white}input{padding:12px;width:min(95%,500px);margin:10px 0;border:1px solid #9cabb2;border-radius:4px}article{background:white;border-left:3px solid #157d87;padding:18px;margin:12px 0}article p{margin:8px 0}small{overflow-wrap:anywhere}figure{margin:20px 0;background:white;padding:12px}img{max-width:100%;height:auto}a{color:#086d86}code{background:#e6edef;padding:3px}footer{margin-top:36px;font-size:13px;color:#536773}
</style><main><div class="meta">CHROMATIN VUS WORKFLOW · RESEARCH REPORT</div><h1>From variant evidence to an assay plan</h1>"""
    html += '<style>td{vertical-align:top}td:nth-child(-n+6){white-space:nowrap}td:nth-last-child(2){min-width:360px}td:nth-last-child(3){min-width:240px}</style>'
    html += demo + f"<p>{n_ok} / {len(rows)} input alleles annotated. Assembly: GRCh38. Run: {e(manifest['created_at'])}.</p>"
    html += "<p>Scores describe model predictions or conservation. They are not measured cellular effects or clinical classifications. NA means unavailable. Research flags use the configured thresholds; this report does not calculate a combined pathogenicity score.</p>"
    html += f"<p>REVEL threshold: {e(manifest['parameters']['revel_threshold'])}; observed AF threshold: {e(manifest['parameters']['af_threshold'])}. AF is the maximum returned dataset-wide exome/genome AF, not ancestry popmax AF.</p>"
    html += '<p>Files: <a href="annotations.csv">annotation CSV</a> · <a href="annotations.json">full JSON</a> · <a href="recommendations.csv">assay suggestions</a> · <a href="provenance.json">provenance</a></p>'
    html += '<h2>Evidence table</h2><label for="filter">Filter rows by gene, variant or status</label><br><input id="filter" type="search" placeholder="e.g. KMT2C"><div class="tablewrap"><table id="variants"><thead><tr>' + headers + "</tr></thead><tbody>" + body + "</tbody></table></div>"
    html += "<h2>Figures</h2><p>Color represents raw values, not a common risk scale. SIFT decreases with predicted impact. SVG versions are available in the figures folder.</p>" + figures
    html += "<h2>Follow-up suggestions</h2><p>These are editable, rule-based research prompts. Select assays based on gene biology, disease mechanism and the available cell model.</p>" + cards
    html += '<footer>Provider responses and hashes are retained in raw/ and provenance.json. Read docs/METHODS.md and docs/DATA_SOURCES.md in the project for assumptions and upstream terms.</footer></main><script>document.getElementById("filter").addEventListener("input",function(){const q=this.value.toLowerCase();document.querySelectorAll("#variants tbody tr").forEach(r=>{r.hidden=!r.textContent.toLowerCase().includes(q)})});</script></html>'
    Path(path).write_text(html, encoding="utf-8")

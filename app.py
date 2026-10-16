"""Optional local upload interface: streamlit run app.py"""
from pathlib import Path
from tempfile import TemporaryDirectory
from io import BytesIO
from zipfile import ZipFile, ZIP_DEFLATED
import json
import streamlit as st
from chromatin_vus.pipeline import run, demo

st.set_page_config(page_title="Chromatin VUS workflow", layout="wide")
st.title("Chromatin VUS workflow")
st.write("Upload a GRCh38 missense SNV list to collect evidence and plan functional assays.")
st.caption("Research workflow · 14-gene starter panel · CSV / TSV / VCF / versioned coding HGVS")
uploaded = st.file_uploader("Variant list", type=["csv", "tsv", "vcf", "txt"])
st.info("Online annotation sends variant alleles to Ensembl and MyVariant.info. Sample genotypes and patient identifiers are not sent by this workflow. Use approved de-identified inputs and follow the data providers' terms.")
context = st.selectbox("Research context", ["neurodevelopment", "general"])
col1, col2 = st.columns(2)
revel = col1.number_input("REVEL research threshold", min_value=0.0, max_value=1.0, value=0.5)
af = col2.number_input("Observed AF research threshold", min_value=0.0, max_value=1.0, value=0.001, format="%.5f")
online = st.button("Annotate uploaded list", disabled=uploaded is None)
synthetic = st.button("Try synthetic demo (no network)")
if online or synthetic:
    try:
        with TemporaryDirectory(prefix="chromatin_vus_") as tmp:
            root = Path(tmp)
            out = root / "report"
            with st.spinner("Preparing evidence and figures..."):
                if synthetic:
                    manifest = demo(out)
                else:
                    input_path = root / ("variants" + Path(uploaded.name).suffix.lower())
                    input_path.write_bytes(uploaded.getvalue())
                    manifest = run(input_path, out, context=context, revel_threshold=revel, af_threshold=af)
                records = json.loads((out / "annotations.json").read_text(encoding="utf-8"))
                archive = BytesIO()
                with ZipFile(archive, "w", ZIP_DEFLATED) as zipped:
                    for path in out.rglob("*"):
                        if path.is_file():
                            zipped.write(path, path.relative_to(out))
                st.session_state["result"] = (manifest, records, archive.getvalue(), [(p.name, p.read_bytes()) for p in sorted((out / "figures").glob("*.png"))])
    except Exception as exc:
        st.error(str(exc))
if "result" in st.session_state:
    manifest, records, archive, figures = st.session_state["result"]
    st.write(f"Annotated {manifest['n_annotated']}/{manifest['n_input_alleles']} alleles")
    st.dataframe([{k: r.get(k) for k in ["input_id", "gene", "hgvsp", "status", "revel", "alphamissense", "sift", "polyphen", "cadd_phred"]} for r in records])
    st.download_button("Download report, tables, figures and provenance", archive, file_name="chromatin_vus_report.zip", mime="application/zip")
    for name, data in figures:
        st.image(data, caption=name)

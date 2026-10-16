# Validation record

Local checks on 2026-10-16, Windows, Python 3.12:

- **33 automated tests passed.** Coverage includes allele/strand conversion, transcript version matching, ambiguous mappings, invalid and multi-ALT input, schema failures, source disagreement, missing versus zero AF, low-confidence labels, cache integrity, offline failure, HTML escaping and CSV formula escaping.
- **Two public KMT2C coding-HGVS examples completed live annotation:** `NM_170606.3:c.14501T>C` and `NM_170606.3:c.14335C>G`. All eight requested score/conservation fields were returned for both. No overall gnomAD AF was returned; missingness remained explicit.
- **Independent genomic representation check:** `7:g.152139219A>G` matched the V4834A coding-HGVS result for genomic identity, coding/protein consequence and five predictors.
- **Offline replay passed** using only the preserved response directory.
- **Optional Streamlit interface:** initialization and its synthetic-demo interaction passed `streamlit.testing.v1.AppTest`; the test produced a downloadable report archive in application state.
- **Distribution build passed:** the wheel contains the bundled gene panel.
- Generated scientific plots and the HTML evidence table were inspected; long heatmap column labels were corrected. Browser row filtering was checked separately.

Live source metadata was Ensembl release **116**, MyVariant build **20250624**, dbNSFP **4.8a** and ClinVar **2025-05**. These describe the returned sources, not the latest available releases of each underlying resource. VEP plugin dataset versions were not exposed by the response. No raw provider response or licensed score table is included in this public repository.

These checks establish that the tested software paths execute and preserve key annotation distinctions. They are not a clinical benchmark, validation across all supported genes, or evidence that predictions accurately estimate the effects of the study's variants. The GitHub Actions matrix is configured for Windows/Linux and Python 3.10/3.12; only completed CI runs should be described as passed.

# Data sources, attribution and reuse

Checked on 2026-10-16. Online services and licenses can change; keep each run's metadata and check the original terms for the intended use. The MIT license covers this repository's code and original documentation. It does not grant rights to third-party annotations.

| Evidence | Access in this workflow | Original resource / reference | Reuse notes |
| --- | --- | --- | --- |
| Consequences, transcripts, domains | [Ensembl VEP REST](https://rest.ensembl.org/documentation/info/vep_hgvs_get) | [McLaren et al., 2016](https://doi.org/10.1186/s13059-016-0974-4) | Follow Ensembl and component terms; plugin data can have separate licenses. |
| REVEL | Selected VEP transcript field | [Ioannidis et al., 2016](https://doi.org/10.1016/j.ajhg.2016.08.016) | VEP documentation specifies noncommercial use. Do not assume MIT permits commercial score reuse. |
| AlphaMissense | Selected VEP transcript field | [Cheng et al., 2023](https://doi.org/10.1126/science.adg7492) | Check the [official repository](https://github.com/google-deepmind/alphamissense) and data terms; model labels are not clinical diagnoses. |
| SIFT | Selected VEP transcript field | [Vaser et al., 2016](https://doi.org/10.1038/nprot.2015.123) | Preserve original attribution and low-confidence annotations. |
| PolyPhen-2 | Selected VEP transcript field | [Adzhubei et al., 2010](https://doi.org/10.1038/nmeth0410-248) | Follow the originating resource's terms. |
| CADD PHRED | Selected VEP transcript field | [CADD](https://cadd.gs.washington.edu/), [Rentzsch et al., 2019](https://doi.org/10.1093/nar/gky1016) | CADD data have their own usage terms; confirm commercial and redistribution permissions. |
| ClinVar, gnomAD and conservation annotations | [MyVariant.info API](https://docs.myvariant.info/en/latest/doc/variant_annotation_service.html) | [MyVariant data sources](https://docs.myvariant.info/en/latest/doc/data.html), [Xin et al., 2016](https://doi.org/10.1186/s13059-016-0953-9) | Aggregator availability does not replace upstream attribution or licenses. |
| ClinVar assertions | MyVariant ClinVar RCV records | [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/), [Landrum et al., 2018](https://doi.org/10.1093/nar/gkx1153) | Inspect submissions, conditions and review status; classifications can change. |
| gnomAD allele frequency | MyVariant overall exome/genome AF | [gnomAD](https://gnomad.broadinstitute.org/), [Karczewski et al., 2020](https://doi.org/10.1038/s41586-020-2308-7) | Report the actual source version; these API fields may lag newer gnomAD releases. |
| phyloP, GERP++, phastCons | MyVariant dbNSFP genomic-site fields | [dbNSFP sources](https://www.dbnsfp.org/download/), [licenses](https://www.dbnsfp.org/license/) | Versions and academic/commercial bundles differ. The smoke query served dbNSFP 4.8a; newer releases do not retroactively identify this response. |
| Gene-to-assay mapping | Manually curated starter panel | [UniProt](https://www.uniprot.org/); individual reviewed entries linked in `panel.csv` | Original short assay suggestions; cite the specific entry and its underlying evidence when used. |

The repository contains public example allele identifiers and synthetic test fixtures. It does not bundle prediction databases, patient records, downloaded annotation caches or unpublished experimental measurements. Local reports include fetched annotations to support reproducibility; check upstream terms before distributing those reports publicly.

## Biological framing and experimental interpretation

- [Rots et al., 2024](https://doi.org/10.1016/j.ajhg.2024.06.009): KMT2C-related neurodevelopmental disorder and variant interpretation. The V4834A example is a reported variant, not a claim of established pathogenicity or a validated positive control.
- [Kubo et al., 2024](https://doi.org/10.1016/j.molcel.2024.02.030): KMT2C/KMT2D catalytic activity, H3K4me1 and neural differentiation provide context for the KMT2C/KMT2D assay suggestions.
- [Brnich et al., 2020](https://doi.org/10.1186/s13073-019-0690-2): evaluating functional evidence for PS3/BS3. This workflow suggests experimental design considerations and does not automatically apply those evidence codes.

## Suggested software citation

Chromatin VUS Workflow contributors. *Chromatin VUS Workflow*, version 0.1.0. https://github.com/Xiaoguo-Zheng/chromatin-vus-workflow . Include the exact commit, access date, parameters and annotation-source versions used in the analysis. No paper acceptance status, journal endorsement or software DOI is implied.

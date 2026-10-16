"""Strict, one-based GRCh38 SNV input handling."""
from dataclasses import dataclass, asdict
from pathlib import Path
import csv
import io
import re


class InputError(ValueError):
    pass


@dataclass
class Variant:
    input_id: str
    hgvs: str = ""
    gene: str = ""
    transcript: str = ""
    original: str = ""
    error: str = ""

    def to_dict(self):
        return asdict(self)


def chromosome(value):
    value = str(value).strip().removeprefix("chr").upper()
    if value == "M":
        value = "MT"
    if value not in {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}:
        raise InputError("Use chromosomes 1–22, X, Y or MT on GRCh38.")
    return value


def genomic_hgvs(chrom, pos, ref, alt):
    chrom = chromosome(chrom)
    try:
        position = int(str(pos))
    except ValueError as exc:
        raise InputError("Position must be a one-based integer.") from exc
    ref, alt = str(ref).upper(), str(alt).upper()
    if position < 1:
        raise InputError("Position must be positive and one-based.")
    if ref not in "ACGT" or alt not in "ACGT" or len(ref) != 1 or len(alt) != 1:
        raise InputError("v0.1 supports SNVs only; indels, symbolic alleles and CNVs are not scored.")
    if ref == alt:
        raise InputError("REF and ALT must differ.")
    return f"{chrom}:g.{position}{ref}>{alt}"


def parse_hgvs(value):
    value = value.strip()
    match = re.fullmatch(r"(?:chr)?([0-9]+|X|Y|M|MT):g\.(\d+)([ACGT])>([ACGT])", value)
    if match:
        return genomic_hgvs(*match.groups()), ""
    match = re.fullmatch(r"((?:NM_\d+|ENST\d+)\.\d+):c\.([0-9]+)([ACGT])>([ACGT])", value)
    if match:
        if match[3] == match[4]:
            raise InputError("REF and ALT must differ.")
        return value, match[1]
    raise InputError("Use GRCh38 chromosome:g.posREF>ALT or versioned NM_/ENST coding SNV HGVS. Protein-only changes are ambiguous and are not accepted in v0.1.")


def read_variants(path, assembly="GRCh38", limit=200):
    if assembly != "GRCh38":
        raise InputError("This release accepts GRCh38/hg38 only; no implicit liftover is performed.")
    path = Path(path)
    if path.stat().st_size > 10_000_000:
        raise InputError("Input exceeds the 10 MB limit; split the variant list.")
    text = path.read_text(encoding="utf-8-sig")
    if not text.strip():
        raise InputError("The input contains no variant records.")
    records = []
    if path.suffix.lower() == ".vcf" or text.startswith("##fileformat=VCF"):
        for line in text.splitlines():
            if line.startswith(("##reference=", "##assembly=")) and re.search(r"hg19|GRCh37|b37", line, re.I):
                raise InputError("The VCF declares a GRCh37 reference; supply a GRCh38 file.")
            if not line or line.startswith("#"):
                continue
            columns = line.split("\t")
            if len(columns) < 5:
                raise InputError("A VCF record needs at least CHROM, POS, ID, REF and ALT.")
            for alt in columns[4].split(","):
                records.append({"chrom": columns[0], "pos": columns[1], "ref": columns[3], "alt": alt,
                                "input_id": columns[2] if columns[2] != "." else ""})
    else:
        delimiter = "\t" if "\t" in text.splitlines()[0] else "," if text.splitlines() and "," in text.splitlines()[0] else None
        if delimiter:
            reader = csv.DictReader(io.StringIO(text), delimiter=delimiter)
            records = [{str(k).strip().lower(): (v or "").strip() for k, v in row.items() if k is not None} for row in reader]
        else:
            records = [{"hgvs": line.strip()} for line in text.splitlines() if line.strip() and not line.startswith("#") and line.strip().lower() != "hgvs"]
    if not records:
        raise InputError("The input contains no variant records.")
    if len(records) > limit:
        raise InputError(f"Input has {len(records)} alleles; the configured limit is {limit}. Split the list or explicitly increase --limit.")
    variants = []
    for number, row in enumerate(records, 1):
        original = row.get("hgvs") or ":".join(row.get(k, "") for k in ("chrom", "pos", "ref", "alt"))
        item = Variant(input_id=f"row{number}" + ("_" + row["input_id"] if row.get("input_id") else ""),
                       gene=row.get("gene", "").upper(), transcript=row.get("transcript", ""), original=original)
        try:
            if row.get("assembly") and row["assembly"] not in {"GRCh38", "hg38"}:
                raise InputError("Per-row assembly is not GRCh38/hg38.")
            if row.get("hgvs"):
                item.hgvs, inferred_transcript = parse_hgvs(row["hgvs"])
                if item.transcript and inferred_transcript and item.transcript != inferred_transcript:
                    raise InputError("The transcript column conflicts with the HGVS transcript.")
                item.transcript = inferred_transcript or item.transcript
            else:
                item.hgvs = genomic_hgvs(row.get("chrom", ""), row.get("pos", ""), row.get("ref", ""), row.get("alt", ""))
            if item.transcript and not re.fullmatch(r"(?:NM_\d+|ENST\d+)\.\d+", item.transcript):
                raise InputError("A transcript must be a versioned NM_ or ENST identifier.")
        except (InputError, ValueError) as exc:
            item.error = str(exc)
        variants.append(item)
    return variants

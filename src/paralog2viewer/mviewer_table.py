# SPDX-License-Identifier: GPL-3.0-or-later
"""Write the 365-column MViewer table from annotation records.

One row is written per allele and per haplotype that carries it (GT 1 in
that haplotype's column of the merged VCF). Allele-level columns come from
the annotation record (:mod:`paralog2viewer.annotate`), haplotype-level
columns (Otherinfo, Genotype) from the haplotype's own call.

Column map (0-based indices of ``resources/mviewer_header.tsv``; "verified"
means checked against the two validation tables read-only, nothing from
them is stored here):

* 0-4 ``Chr``, ``Start``, ``End``, ``Ref``, ``Alt``: ANNOVAR values with
  ``chr`` removed (verified: bare ``1``; indels use ANNOVAR's ``-``
  convention and shifted Start). Without an ANNOVAR row, the VCF values.
* 5-257: ANNOVAR columns copied by case-insensitive name (verified: all 253
  header names exist in the native multianno header, same order).
* 258-270 ``Otherinfo1-13``: zygosity (``hom`` for a haploid GT 1, verified),
  QUAL as a number or ``0`` (verified ``0`` for QUAL ``.``), INFO ``DP`` or
  ``.``, then the VCF columns CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO,
  FORMAT and this haplotype's sample string (verified column order). The
  values come from the merged VCF record. A missing FILTER is written as
  ``PASS``, because only passing calls enter the merged VCF.
* 271-296 legacy names: ``1000G_*`` from ``1000g2015aug_*``, ``ExAC_Freq``
  from ``ExAC_ALL``, ``ESP6500si_*`` from ``esp6500siv2_*``,
  ``gnomAD_exome_*`` and ``gnomAD_genome_*`` from the ``gnomad41_*``
  columns (verified, 0 differences: ``ALL`` from ``AF``, ``OTH`` from
  ``AF_remaining``, the others from ``AF_<population>``; the older
  ``gnomAD_*`` ANNOVAR columns are not the source), ``COSMIC_ID`` and
  ``COSMIC_DIS`` from ``cosmic70`` (``ID=`` and ``OCCURENCE=``, commas
  become ``;`` in both), ``AAChange.refgene`` from
  ``AAChange.refGeneWithVer`` and ``GeneDetail.refGeneWithVer`` (the parts
  that are not ``.``, joined by ``/``; verified), ``dbSNP`` from
  ``avsnp147`` (verified equal to
  ``dbSNP ID``; the VCF ID replaces ``.``, as in the validation tables),
  ``ClinVar_SIG/DIS/DB`` from ``CLNSIG/CLNDN/CLNDISDB`` (verified, 100 %
  equal) and ``MutationAssessor_score_rankscore`` from
  ``MutationAssessor_rankscore`` (verified). The 1000G, ExAC, ESP6500 and
  COSMIC mappings follow the ANNOVAR protocol names (not joined to values,
  the native columns are not in the validation tables). ``CG46``,
  ``ClinVar_ID``, ``ClinVar_DBID``, ``GWAS_*`` and the other legacy
  columns stay empty (verified empty).
* 318 ``Consequence``: VEP consequence of the chosen transcript with ``;``
  between terms; for symbolic alleles the sorted union over all VEP rows
  (the validation tables show that form for the four symbolic rows).
* 319 ``Gene``: distinct VEP symbols in order of appearance, else the chosen
  gene, else the ANNOVAR refGene genes. Verified only in kind: the
  validation tables hold a union of several gene sets, whose rule is not
  recoverable from them.
* 320 ``Variant`` ``Ref>Alt/Alt`` (verified, VCF ref and alt), 321
  ``Coordinate`` = ANNOVAR Start (verified), 322 ``dbSNP ID``, 323 ``Type``
  (verified: snv, deletion, insertion from ANNOVAR Ref/Alt; ``other`` is not
  seen in the validation tables), 324-326 Panel, Phenotypes, Inheritance,
  328 HGMD class, 329-330 ClinVar and ClinVar_VID, 333 Taiwan Biobank, 345
  InterVar text, 350 Interpretation Result (class), all from the record.
* 353 ``Max Allele Freq``: maximum over ``1000G_*``, ``ExAC_*``,
  ``ESP6500si_*``, ``gnomAD_exome_*``, ``gnomAD_genome_*`` and ``Taiwan
  Biobank``, formatted like .NET (``3.086E-05``, ``1``); ``0`` when no
  value. Verified: exact match on all explicit rows of both tables.
* 354-357 Transcript, Exon (``exonN``), Nucleotide (HGVS c.), AAChange
  (HGVS p.) from the record's chosen transcript. Not comparable: the
  validation tables fill them for about 14 % of rows only.
* 359 ``Genotype``: the exact HP value. 360 ``Severity Score``: empty.
* All other columns (legacy, Methylated loci, CentoMD, PharmGKB, MaxEntScan,
  SpliceAI, family genotypes, ...) are empty (verified empty).
* Symbolic and breakend alleles: Chr, Start, End, Ref, Alt, Consequence,
  Gene and Genotype only. The validation tables also fill ``Func.refgene``
  and ``Gene.refgene`` for these rows; the project rule leaves them empty.

Coordinates: ``Start`` and ``End`` are 1-based, closed, as ANNOVAR writes
them.
"""

from __future__ import annotations

import csv
import gzip
import io
import re
from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from urllib.parse import unquote

import pysam

from paralog2viewer.annotate import read_records
from paralog2viewer.datalayout import allele_key

HEADER_RESOURCE = Path(__file__).parent / "resources" / "mviewer_header.tsv"
N_COLUMNS = 365
GENOTYPE_INDEX = 359

# gnomAD population label -> suffix of the gnomad41 ``AF`` column.
GNOMAD41_POPULATIONS: dict[str, str] = {
    "ALL": "",
    "AFR": "_afr",
    "AMR": "_amr",
    "ASJ": "_asj",
    "EAS": "_eas",
    "FIN": "_fin",
    "NFE": "_nfe",
    "OTH": "_remaining",
    "SAS": "_sas",
}

# Legacy-named column -> native ANNOVAR column (case-insensitive).
LEGACY_FROM_ANNOVAR: dict[str, str] = {
    "1000G_ALL": "1000g2015aug_all",
    "1000G_AFR": "1000g2015aug_afr",
    "1000G_AMR": "1000g2015aug_amr",
    "1000G_EAS": "1000g2015aug_eas",
    "1000G_SAS": "1000g2015aug_sas",
    "1000G_EUR": "1000g2015aug_eur",
    "ExAC_Freq": "ExAC_ALL",
    "ESP6500si_ALL": "esp6500siv2_all",
    "ESP6500si_AA": "esp6500siv2_aa",
    "ESP6500si_EA": "esp6500siv2_ea",
    **{
        f"gnomAD_{kind}_{pop}": f"gnomad41_{kind}_AF{suffix}"
        for kind in ("exome", "genome")
        for pop, suffix in GNOMAD41_POPULATIONS.items()
    },
    "ClinVar_SIG": "CLNSIG",
    "ClinVar_DIS": "CLNDN",
    "ClinVar_DB": "CLNDISDB",
    "MutationAssessor_score_rankscore": "MutationAssessor_rankscore",
}
FREQ_PREFIXES = (
    "1000G_",
    "ExAC_",
    "ESP6500si_",
    "gnomAD_exome_",
    "gnomAD_genome_",
)

# Columns the writer rebuilds by rule or leaves empty; compare_tables
# ignores them by default.
DEFAULT_IGNORED_COLUMNS: tuple[str, ...] = (
    "Consequence",
    "Gene",
    "Variant",
    "Type",
    "Max Allele Freq",
    "Transcript",
    "Exon",
    "Nucleotide",
    "AAChange",
    "Severity Score",
    # Documented deviations: richer INFO/FORMAT/sample from the merged VCF,
    # ID '.', and the caller IDs the old pipeline wrote into dbSNP.
    "Otherinfo6",
    "Otherinfo11",
    "Otherinfo12",
    "Otherinfo13",
    "dbSNP",
    "dbSNP ID",
)
# Compared on small variants only: for symbolic alleles the baseline fills
# them from ANNOVAR, the project rule leaves them empty.
SYMBOLIC_IGNORED_COLUMNS: tuple[str, ...] = ("Func.refgene", "Gene.refgene")
DEFAULT_SORT_COLUMNS: tuple[str, ...] = ("rmsk",)

_SYMBOLIC_ALT = re.compile(r"[<>\[\]]")
_DIGITS = re.compile(r"(\d+)")
_CHROM_RANK = {"X": 23, "Y": 24, "M": 25, "MT": 25}


def load_header(path: Path | None = None) -> list[str]:
    """Read and check the 365 column names.

    Args:
        path: Header file (one tab-separated line); the packaged copy when
            None.

    Raises:
        ValueError: If the header is not 365 unique names with ``Genotype``
            at index 359.
    """
    text = Path(path or HEADER_RESOURCE).read_text(encoding="utf-8")
    header = text.rstrip("\n").split("\t")
    if (
        len(header) != N_COLUMNS
        or header[GENOTYPE_INDEX] != "Genotype"
        or len(set(header)) != N_COLUMNS
    ):
        raise ValueError("invalid MViewer header")
    return header


def natural_key(text: str) -> tuple:
    """Return a sort key that orders digit runs numerically."""
    return tuple(
        (0, int(p), "") if p.isdigit() else (1, 0, p.lower())
        for p in _DIGITS.split(text)
        if p
    )


def chrom_key(chrom: str) -> tuple:
    """Return the sort key chr1..chr22, chrX, chrY, chrM, others after."""
    name = chrom.removeprefix("chr")
    if name.isdigit() and 1 <= int(name) <= 22:
        return (int(name), ())
    if name.upper() in _CHROM_RANK:
        return (_CHROM_RANK[name.upper()], ())
    return (26, natural_key(name))


def net_float(value: float) -> str:
    """Format a number the way MViewer (.NET) does.

    ``0.0075`` -> ``0.0075``, ``3.086e-05`` -> ``3.086E-05``, ``1.0`` ->
    ``1``.
    """
    return repr(value).upper().removesuffix(".0")


def _number(text: str) -> float | None:
    """Return the float of a table cell, or None for blanks and ``.``."""
    try:
        return float(text)
    except ValueError:
        return None


def max_allele_freq(cells: Mapping[str, str]) -> str:
    """Return ``Max Allele Freq`` from the frequency cells of a row.

    Args:
        cells: Column name -> value of the row.

    Returns:
        The maximum over the columns named by ``FREQ_PREFIXES`` and
        ``Taiwan Biobank`` as .NET-style text; ``0`` if none is numeric.
    """
    values = [
        v
        for name, text in cells.items()
        if name.startswith(FREQ_PREFIXES) or name == "Taiwan Biobank"
        if (v := _number(text)) is not None
    ]
    top = max(values, default=0.0)
    return "0" if top == 0 else net_float(top)


def variant_type(ref: str, alt: str) -> str:
    """Classify an allele as snv, insertion, deletion or other.

    ANNOVAR style (``-`` for the empty allele) and VCF style are both
    accepted.
    """
    if ref == "-":
        return "insertion"
    if alt == "-":
        return "deletion"
    if len(ref) == 1 and len(alt) == 1:
        return "snv"
    if len(alt) > len(ref) and alt.startswith(ref):
        return "insertion"
    if len(ref) > len(alt) and ref.startswith(alt):
        return "deletion"
    return "other"


def _after_colon(text: str) -> str:
    """Return the part of a HGVS string after the transcript prefix."""
    return unquote(text.split(":", 1)[1]) if ":" in text else unquote(text)


def _exon_label(exon: str) -> str:
    """Convert VEP ``4/10`` to ``exon4``; empty when absent."""
    first = exon.split("/")[0].split("-")[0]
    return f"exon{first}" if first.isdigit() else ""


def _vep_symbols(record: Mapping) -> list[str]:
    """Return distinct VEP gene symbols in order of appearance."""
    symbols = (
        t["symbol"]
        for t in record["vep_transcripts"]
        if t["feature_type"] == "Transcript" and t["symbol"]
    )
    return list(dict.fromkeys(symbols))


def _gene(record: Mapping) -> str:
    """Return the ``Gene`` cell for a record."""
    symbols = _vep_symbols(record)
    if not symbols:
        symbols = [record["gene"]] if record["gene"] else list(record["genes"])
    return ";".join(symbols)


def _symbolic_consequence(record: Mapping) -> str:
    """Return the sorted union of VEP consequence terms of all rows."""
    terms = {
        term
        for t in record["vep_transcripts"]
        for term in t["consequence"].split(",")
        if term
    }
    return ";".join(sorted(terms))


def _cosmic(text: str) -> tuple[str, str]:
    """Split ANNOVAR ``cosmic70`` into the COSMIC_ID and COSMIC_DIS cells."""
    if text in ("", "."):
        return "", ""
    fields = dict(p.split("=", 1) for p in text.split(";") if "=" in p)
    return (
        fields.get("ID", "").replace(",", ";"),
        fields.get("OCCURENCE", "").replace(",", ";"),
    )


def _aachange_refgene(annovar: Mapping[str, str]) -> str:
    """Return the ``AAChange.refgene`` cell from lower-case ANNOVAR keys.

    The baseline joins the versioned ``AAChange`` and ``GeneDetail`` texts
    (those that are not ``.``) with ``/``; ``.`` when both are missing. The
    plain ``AAChange.refGene`` text is kept when the versioned columns are
    absent.
    """
    keys = ("aachange.refgenewithver", "genedetail.refgenewithver")
    if not any(k in annovar for k in keys):
        return annovar.get("aachange.refgene", ".")
    parts = [annovar[k] for k in keys if annovar.get(k, ".") not in ("", ".")]
    return "/".join(parts) or "."


class _Writer:
    """Holds the header index and builds rows."""

    def __init__(self, header: Sequence[str]) -> None:
        self.header = list(header)
        self.index = {name: i for i, name in enumerate(self.header)}
        self.annovar_names = self.header[5:258]
        self.freq_names = [
            n for n in self.header if n.startswith(FREQ_PREFIXES)
        ]

    def _put(self, row: list[str], name: str, value: str) -> None:
        row[self.index[name]] = value

    def row(
        self,
        record: Mapping,
        vcf: Mapping,
        hap: str,
    ) -> list[str]:
        """Build one 365-cell row.

        Args:
            record: Annotation record of the allele.
            vcf: Haplotype-level VCF values: ``chrom``, ``pos``, ``end``,
                ``ref``, ``alt``, ``id``, ``qual``, ``filter``, ``info``,
                ``format``, ``sample`` and ``gt`` (tuple).
            hap: Exact HP value.
        """
        row = [""] * N_COLUMNS
        put = self._put
        annovar = {k.lower(): v for k, v in record["annovar"].items()}
        symbolic = bool(record["symbolic"])
        if annovar:
            put(row, "Chr", annovar["chr"].removeprefix("chr"))
            for name in ("Start", "End", "Ref", "Alt"):
                put(row, name, annovar[name.lower()])
            for name in self.annovar_names:
                put(row, name, annovar.get(name.lower(), "."))
        else:
            put(row, "Chr", vcf["chrom"].removeprefix("chr"))
            put(row, "Start", str(vcf["pos"]))
            put(row, "End", str(vcf["end"]))
            put(row, "Ref", vcf["ref"])
            put(row, "Alt", vcf["alt"])
        put(row, "Genotype", hap)
        if symbolic:
            put(row, "Consequence", _symbolic_consequence(record))
            put(row, "Gene", _gene(record))
            return row

        self._otherinfo(row, vcf)
        self._legacy(row, annovar, vcf)
        ref, alt = vcf["ref"], vcf["alt"]
        put(row, "Gene", _gene(record))
        put(row, "Consequence", record["consequence"].replace(",", ";"))
        put(row, "Variant", f"{ref}>{alt}/{alt}")
        put(row, "Coordinate", row[self.index["Start"]])
        put(row, "dbSNP ID", row[self.index["dbSNP"]])
        if annovar:
            kind = variant_type(annovar["ref"], annovar["alt"])
        else:
            kind = variant_type(ref, alt)
        put(row, "Type", kind)
        put(row, "Panel", record["panel"])
        put(row, "Phenotypes", record["phenotypes"])
        put(row, "Inheritance", record["inheritance"])
        put(row, "HGMD Variant Class", record["hgmd_class"])
        put(row, "ClinVar", record["clinvar"])
        put(row, "ClinVar_VID", record["clinvar_vid"])
        put(row, "Taiwan Biobank", record["taiwan_biobank"])
        put(row, "InterVar", record["intervar"])
        put(row, "Interpretation Result", record["intervar_class"])
        put(
            row,
            "Max Allele Freq",
            max_allele_freq(
                {n: row[self.index[n]] for n in self.freq_names}
                | {"Taiwan Biobank": record["taiwan_biobank"]}
            ),
        )
        if record["transcript"]:
            put(row, "Transcript", record["transcript"])
            put(row, "Exon", _exon_label(record["exon"]))
            put(row, "Nucleotide", _after_colon(record["hgvsc"]))
            put(row, "AAChange", _after_colon(record["hgvsp"]))
        return row

    def _otherinfo(self, row: list[str], vcf: Mapping) -> None:
        """Fill Otherinfo1-13 from the haplotype's own VCF call."""
        gt = vcf["gt"]
        qual = _number(vcf["qual"])
        depth = re.search(r"(?:^|;)DP=([^;]+)", vcf["info"])
        values = [
            "hom" if len(set(gt)) == 1 else "het",
            "0" if qual is None else vcf["qual"],
            depth.group(1) if depth else ".",
            vcf["chrom"],
            str(vcf["pos"]),
            vcf["id"],
            vcf["ref"],
            vcf["alt"],
            vcf["qual"],
            vcf["filter"],
            vcf["info"],
            vcf["format"],
            vcf["sample"],
        ]
        for number, value in enumerate(values, start=1):
            self._put(row, f"Otherinfo{number}", value)

    def _legacy(
        self, row: list[str], annovar: Mapping[str, str], vcf: Mapping
    ) -> None:
        """Fill the legacy-named frequency and ClinVar columns."""
        if not annovar:
            self._put(row, "dbSNP", vcf["id"])
            return
        for name, native in LEGACY_FROM_ANNOVAR.items():
            if name not in self.index:  # no gnomAD_genome_SAS column
                continue
            self._put(row, name, annovar.get(native.lower(), "."))
        cosmic_id, cosmic_dis = _cosmic(annovar.get("cosmic70", "."))
        self._put(row, "COSMIC_ID", cosmic_id)
        self._put(row, "COSMIC_DIS", cosmic_dis)
        self._put(row, "AAChange.refgene", _aachange_refgene(annovar))
        rsid = annovar.get("avsnp147", ".")
        self._put(row, "dbSNP", vcf["id"] if rsid in ("", ".") else rsid)


def _haplotype_calls(
    merged_vcf: Path,
) -> Iterable[tuple[dict, str]]:
    """Yield (VCF values, HP name) for every allele a haplotype carries."""
    with pysam.VariantFile(str(merged_vcf)) as vcf:
        samples = list(vcf.header.samples)
        for rec in vcf:
            alts = rec.alts or ()
            cols = None
            for sample_index, name in enumerate(samples):
                gt = rec.samples[name]["GT"]
                if not gt or None in gt or not any(g > 0 for g in gt):
                    continue
                if cols is None:
                    cols = str(rec).rstrip("\n").split("\t")
                filter_ = cols[6] if cols[6] not in ("", ".") else "PASS"
                for alt_index in sorted({g for g in gt if g > 0}):
                    yield (
                        {
                            "chrom": rec.chrom,
                            "pos": rec.pos,
                            "end": rec.stop,
                            "ref": rec.ref,
                            "alt": alts[alt_index - 1],
                            "id": rec.id or ".",
                            "qual": cols[5],
                            "filter": filter_,
                            "info": cols[7],
                            "format": cols[8],
                            "sample": cols[9 + sample_index],
                            "gt": gt,
                        },
                        name,
                    )


def write_mviewer_table(
    records_path: Path,
    merged_vcf: Path,
    sample_id: str,
    sex: str,
    out_path: Path,
    header_path: Path | None = None,
) -> int:
    """Write the MViewer table.

    Args:
        records_path: ``records.jsonl.gz`` from the annotation step.
        merged_vcf: Merged VCF; haplotypes are sample columns named by HP.
        sample_id: Sample name (kept for the caller's file naming and
            messages; not written into the table).
        sex: ``Female``, ``Male`` or ``Unknown`` (case-insensitive).
        out_path: Output path (gzip text), for example
            ``{SAMPLE_ID}.paralog.tsv.gz``.
        header_path: Header file; the packaged 365-name header when None.

    Returns:
        Number of data rows written.

    Raises:
        ValueError: On an unknown sex or an allele without an annotation
            record.
    """
    if sex.capitalize() not in ("Female", "Male", "Unknown"):
        raise ValueError(
            f"sex must be Female, Male or Unknown, got {sex!r} ({sample_id})"
        )
    header = load_header(header_path)
    records = {r["allele_key"]: r for r in read_records(Path(records_path))}
    writer = _Writer(header)
    rows: list[tuple[tuple, list[str]]] = []
    for vcf, hap in _haplotype_calls(Path(merged_vcf)):
        key = allele_key(vcf["chrom"], vcf["pos"], vcf["ref"], vcf["alt"])
        if key not in records:
            raise ValueError(f"no annotation record for {key}")
        sort_key = (
            chrom_key(vcf["chrom"]),
            vcf["pos"],
            natural_key(hap),
            vcf["ref"],
            vcf["alt"],
        )
        rows.append((sort_key, writer.row(records[key], vcf, hap)))
    rows.sort(key=lambda item: item[0])
    out_path = Path(out_path)
    with out_path.open("wb") as raw:
        # mtime 0 and no name: identical input gives identical bytes.
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
            with io.TextIOWrapper(gz, encoding="utf-8", newline="") as text:
                text.write(f"#keywords=\n#sex={sex.capitalize()}\n")
                out = csv.writer(text, delimiter="\t", lineterminator="\n")
                out.writerow(header)
                out.writerows(row for _, row in rows)
    return len(rows)


def _read_table(path: Path) -> tuple[list[str], list[list[str]]]:
    """Read a table with two comment lines; return header and rows."""
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        lines = handle.read().split("\n")
    header = lines[2].split("\t")
    return header, [x.split("\t") for x in lines[3:] if x]


def _normalize(name: str, value: str, sort_columns: Sequence[str]) -> str:
    """Return a cell value, with repeat names sorted for sort columns."""
    if name in sort_columns and value.startswith("Name="):
        return "Name=" + ",".join(sorted(value[5:].split(",")))
    return value


def compare_tables(
    baseline_gz: Path,
    new_gz: Path,
    ignore_columns: Sequence[str] = DEFAULT_IGNORED_COLUMNS,
    sort_columns: Sequence[str] = DEFAULT_SORT_COLUMNS,
    symbolic_ignore_columns: Sequence[str] = SYMBOLIC_IGNORED_COLUMNS,
) -> dict[str, int]:
    """Count per-column differences between two MViewer tables.

    Rows are matched by (Chr, Start, Ref, Alt, Genotype); duplicate keys
    are paired in file order. Cells of ``sort_columns`` that start with
    ``Name=`` have their comma-separated names sorted before comparing
    (their order varies between runs of the same pipeline).

    Args:
        baseline_gz: Reference table (gzip).
        new_gz: Table to check (gzip).
        ignore_columns: Columns not compared. The default is the nine
            rebuilt columns, ``Severity Score`` and the documented
            deviations (``Otherinfo6/11/12/13``, ``dbSNP``, ``dbSNP ID``).
        sort_columns: Columns whose repeat names are sorted first.
        symbolic_ignore_columns: Columns not compared on rows whose Alt is
            symbolic (``<DEL>``) or a breakend (brackets).

    Returns:
        Column name -> number of matched rows that differ (columns without
        a difference are present with 0), plus ``__rows_baseline__``,
        ``__rows_new__``, ``__unmatched_baseline__`` and
        ``__unmatched_new__``.

    Raises:
        ValueError: If the two headers differ.
    """
    head_a, rows_a = _read_table(Path(baseline_gz))
    head_b, rows_b = _read_table(Path(new_gz))
    if head_a != head_b:
        raise ValueError("table headers differ")
    pos = {n: i for i, n in enumerate(head_a)}
    key_cols = [pos[n] for n in ("Chr", "Start", "Ref", "Alt", "Genotype")]

    def keyed(rows: list[list[str]]) -> dict[tuple, list[list[str]]]:
        table: dict[tuple, list[list[str]]] = defaultdict(list)
        for r in rows:
            table[tuple(r[i] for i in key_cols)].append(r)
        return table

    table_a, table_b = keyed(rows_a), keyed(rows_b)
    ignored = set(ignore_columns)
    compared = [n for n in head_a if n not in ignored]
    symbolic_skip = {pos[n] for n in symbolic_ignore_columns if n in pos}
    alt_col = pos["Alt"]
    diff: Counter[str] = Counter({n: 0 for n in compared})
    unmatched_a = unmatched_b = 0
    for key in table_a.keys() | table_b.keys():
        left, right = table_a.get(key, []), table_b.get(key, [])
        unmatched_a += max(0, len(left) - len(right))
        unmatched_b += max(0, len(right) - len(left))
        for ra, rb in zip(left, right, strict=False):
            symbolic = bool(_SYMBOLIC_ALT.search(ra[alt_col]))
            for name in compared:
                i = pos[name]
                if symbolic and i in symbolic_skip:
                    continue
                if _normalize(name, ra[i], sort_columns) != _normalize(
                    name, rb[i], sort_columns
                ):
                    diff[name] += 1
    result = dict(diff)
    result["__rows_baseline__"] = len(rows_a)
    result["__rows_new__"] = len(rows_b)
    result["__unmatched_baseline__"] = unmatched_a
    result["__unmatched_new__"] = unmatched_b
    return result

# SPDX-License-Identifier: GPL-3.0-or-later
"""Write the merged, annotated and indexed VCF.

The output is a copy of ``merged.vcf.gz`` with these INFO fields added:

* ``CSQ``: the VEP consequence field, one entry per VEP row, subfields as
  listed in the header description (``Allele|Consequence|...``), including
  the MaxEntScan and AlphaMissense plugin values. Terms inside a field are
  separated by ``&``; characters that VCF forbids are percent-encoded.
* ``INTERVAR``: the InterVar class (spaces written as ``_``).
* ``CLNSIG``: the ClinVar significance (spaces written as ``_``).
* mode ``full`` also adds ``ANN_*`` fields with selected ANNOVAR columns.

HGMD and Taiwan Biobank values are licensed or restricted and are never
written, in any mode.
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path

import pysam

from paralog2viewer.annotate import read_records
from paralog2viewer.datalayout import allele_key

MODES = ("default", "full")
CSQ_FIELDS: tuple[str, ...] = (
    "Allele",
    "Consequence",
    "IMPACT",
    "SYMBOL",
    "Gene",
    "Feature_type",
    "Feature",
    "BIOTYPE",
    "EXON",
    "INTRON",
    "HGVSc",
    "HGVSp",
    "cDNA_position",
    "CDS_position",
    "Protein_position",
    "Amino_acids",
    "Codons",
    "CANONICAL",
    "MANE_SELECT",
    "MANE_PLUS_CLINICAL",
)
# Plugin fields always present in the CSQ header; further keys seen in the
# records are appended in sorted order.
PLUGIN_FIELDS: tuple[str, ...] = (
    "MaxEntScan_alt",
    "MaxEntScan_diff",
    "MaxEntScan_ref",
    "am_class",
    "am_pathogenicity",
)
# Record key of each CSQ subfield that is not a plugin value.
CSQ_SOURCE = {
    "Consequence": "consequence",
    "IMPACT": "impact",
    "SYMBOL": "symbol",
    "Gene": "gene_id",
    "Feature_type": "feature_type",
    "Feature": "feature",
    "BIOTYPE": "biotype",
    "EXON": "exon",
    "INTRON": "intron",
    "HGVSc": "hgvsc",
    "HGVSp": "hgvsp",
    "cDNA_position": "cdna_position",
    "CDS_position": "cds_position",
    "Protein_position": "protein_position",
    "Amino_acids": "amino_acids",
    "Codons": "codons",
    "CANONICAL": "canonical",
    "MANE_SELECT": "mane_select",
    "MANE_PLUS_CLINICAL": "mane_plus_clinical",
}
# ANNOVAR columns added in mode "full" (public-source columns only).
FULL_ANNOVAR_FIELDS: tuple[str, ...] = (
    "Func.refGene",
    "Gene.refGene",
    "ExonicFunc.refGene",
    "AAChange.refGene",
    "avsnp147",
    "gnomad41_exome_AF",
    "gnomad41_genome_AF",
    "REVEL_score",
    "CADD_phred",
    "AlphaMissense_score",
    "AlphaMissense_pred",
)
# Never written to the VCF.
WITHHELD_RECORD_KEYS = ("hgmd_class", "taiwan_biobank")

_BAD_CHARS = re.compile(r"[;=,|\s]|%(?![0-9A-Fa-f]{2})")


def _escape(value: str) -> str:
    """Percent-encode characters that VCF INFO or CSQ cannot hold."""
    return _BAD_CHARS.sub(lambda m: f"%{ord(m.group()):02X}", value)


def _ann_tag(column: str) -> str:
    """Return the INFO ID for an ANNOVAR column (``ANN_`` + safe name)."""
    return "ANN_" + re.sub(r"[^0-9A-Za-z_]", "_", column)


def plugin_fields(records: Iterable[Mapping]) -> list[str]:
    """Return the plugin subfield names: defaults first, others sorted."""
    seen = {
        key
        for rec in records
        for tx in rec["vep_transcripts"]
        for key in tx["plugins"]
    }
    return list(PLUGIN_FIELDS) + sorted(seen - set(PLUGIN_FIELDS))


def csq_entries(
    record: Mapping, alt: str, fields: Sequence[str]
) -> tuple[str, ...]:
    """Return the CSQ entries of one allele.

    Args:
        record: Annotation record.
        alt: ALT allele written in the Allele subfield.
        fields: All CSQ subfield names, in header order.

    Returns:
        One pipe-separated entry per VEP row.
    """
    entries = []
    for tx in record["vep_transcripts"]:
        cells = []
        for name in fields:
            if name == "Allele":
                value = alt
            elif name in CSQ_SOURCE:
                value = tx[CSQ_SOURCE[name]]
                if name == "Consequence":
                    value = value.replace(",", "&")
            else:
                value = tx["plugins"].get(name, "")
            cells.append(_escape(value))
        entries.append("|".join(cells))
    return tuple(entries)


def _clnsig(record: Mapping) -> str:
    """Return the ClinVar significance of a record, or ``""``."""
    return record["clinvar"] or _annovar_value(record, "CLNSIG")


def _tag_value(text: str) -> str:
    """Return a string INFO value: spaces become ``_``, specials encoded."""
    return _escape(text.replace(" ", "_"))


def _header_lines(fields: Sequence[str], mode: str) -> list[str]:
    """Return the INFO header lines this module adds."""
    lines = [
        '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequences '
        "annotated by Ensembl VEP. Format: " + "|".join(fields) + '">',
        '##INFO=<ID=INTERVAR,Number=A,Type=String,Description="InterVar '
        'class of the ALT allele">',
        '##INFO=<ID=CLNSIG,Number=A,Type=String,Description="ClinVar '
        'clinical significance of the ALT allele">',
    ]
    if mode == "full":
        for column in FULL_ANNOVAR_FIELDS:
            lines.append(
                f"##INFO=<ID={_ann_tag(column)},Number=A,Type=String,"
                f'Description="ANNOVAR column {column}">'
            )
    return lines


def write_annotated_vcf(
    merged_vcf: Path,
    records_path: Path,
    out_path: Path,
    mode: str = "default",
) -> Path:
    """Write the annotated, bgzip-compressed and tabix-indexed VCF.

    Args:
        merged_vcf: ``merged.vcf.gz`` (haplotypes as sample columns).
        records_path: ``records.jsonl.gz`` from the annotation step.
        out_path: Output path, for example ``{SAMPLE_ID}.paralog.vcf.gz``.
        mode: ``default`` (CSQ, INTERVAR, CLNSIG) or ``full`` (also ANNOVAR
            fields). HGMD and Taiwan Biobank are never written.

    Returns:
        ``out_path``; the index is ``out_path`` + ``.tbi``.

    Raises:
        ValueError: On an unknown mode or an allele without a record.
    """
    if mode not in MODES:
        raise ValueError(f"mode must be one of {MODES}, got {mode!r}")
    records = {r["allele_key"]: r for r in read_records(Path(records_path))}
    fields = [*CSQ_FIELDS, *plugin_fields(records.values())]
    out_path = Path(out_path)
    with pysam.VariantFile(str(merged_vcf)) as src:
        header = src.header.copy()
        for line in _header_lines(fields, mode):
            header.add_line(line)
        with pysam.VariantFile(str(out_path), "wz", header=header) as out:
            for rec in src:
                alts = rec.alts or ()
                found = []
                for alt in alts:
                    key = allele_key(rec.chrom, rec.pos, rec.ref, alt)
                    if key not in records:
                        raise ValueError(f"no annotation record for {key}")
                    found.append(records[key])
                new = rec.copy()
                new.translate(header)
                csq = tuple(
                    e
                    for r, alt in zip(found, alts, strict=True)
                    for e in csq_entries(r, alt, fields)
                )
                if csq:
                    new.info["CSQ"] = csq
                _set_per_alt(
                    new, "INTERVAR", [r["intervar_class"] for r in found]
                )
                _set_per_alt(new, "CLNSIG", [_clnsig(r) for r in found])
                if mode == "full":
                    for column in FULL_ANNOVAR_FIELDS:
                        _set_per_alt(
                            new,
                            _ann_tag(column),
                            [_annovar_value(r, column) for r in found],
                        )
                out.write(new)
    pysam.tabix_index(str(out_path), preset="vcf", force=True)
    return out_path


def _annovar_value(record: Mapping, column: str) -> str:
    """Return an ANNOVAR column of a record; ``.`` counts as missing."""
    value = next(
        (
            v
            for k, v in record["annovar"].items()
            if k.lower() == column.lower()
        ),
        "",
    )
    return "" if value == "." else value


def _set_per_alt(
    rec: pysam.VariantRecord, tag: str, values: Sequence[str]
) -> None:
    """Set a Number=A string INFO field; skip when every value is empty."""
    if any(values):
        rec.info[tag] = tuple(_tag_value(v) if v else "." for v in values)

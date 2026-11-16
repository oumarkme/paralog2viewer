# SPDX-License-Identifier: GPL-3.0-or-later
"""Annotation chain: one annotation record per unique allele.

Flow of :func:`annotate_alleles`:

1. :func:`sites_only_vcf` writes sites-only VCFs from the merged VCF, one
   line per unique allele with ``ID`` = ``datalayout.allele_key``. Explicit
   sequence alleles go to ``small.vcf``; symbolic and breakend alleles go to
   ``symbolic.vcf``.
2. A SQLite cache (:class:`RecordCache`) returns finished records keyed by
   (allele key, SHA-256 of the data manifest, record schema version and
   paralog2viewer version). Only missing alleles continue.
3. VEP (small variants with ``--hgvs``) and ANNOVAR ``-vcfinput`` run at the
   same time in two threads; InterVar runs after ANNOVAR; symbolic alleles
   go through VEP only, after the small-variant VEP run.
4. Records are joined and stored; ``W/annotation/records.jsonl.gz`` has one
   JSON object per allele in the order of first appearance in the merged VCF.

Record schema (``RECORD_FIELDS``; every key is always present; absent values
are ``""``, ``[]`` or ``{}``):

* ``allele_key``, ``chrom``, ``pos`` (1-based), ``ref``, ``alt``,
  ``symbolic`` (bool).
* ``annovar``: dict of ANNOVAR native column name -> string, in file order:
  Chr, Start, End, Ref, Alt, then the protocol columns, without the
  ``Otherinfo*`` columns (those are per haplotype, not per allele). Native
  names differ in case from ``mviewer_header.tsv`` (for example
  ``Func.refGene`` versus ``Func.refgene``); match case-insensitively.
  Empty ``{}`` for symbolic alleles and for alleles ANNOVAR dropped.
* ``intervar``: the InterVar summary string as written by InterVar (leading
  space kept), ``intervar_class``: Pathogenic, Likely pathogenic, Uncertain
  significance, Likely benign, Benign, or ``""``.
* ``vep_variant``: dict of allele-level VEP fields (``Existing_variation``,
  ``VARIANT_CLASS``, ``MAX_AF``, ``*_AF``, ``CLIN_SIG``, ``PUBMED`` ...).
* ``vep_transcripts``: list of dicts, one per VEP output row, with
  ``feature``, ``feature_type``, ``gene_id``, ``symbol``, ``consequence``
  (VEP comma-separated terms; the baseline table uses ``;``), ``impact``,
  ``biotype``, ``canonical``, ``mane_select``, ``mane_plus_clinical``,
  ``exon`` (VEP ``n/total``), ``intron``, ``hgvsc``, ``hgvsp`` (VEP
  strings with the transcript prefix), ``cdna_position``, ``cds_position``,
  ``protein_position``, ``amino_acids``, ``codons``, and ``plugins`` (the
  MaxEntScan_* and AlphaMissense ``am_class``/``am_pathogenicity`` keys).
* ``clinvar``, ``clinvar_vid``, ``taiwan_biobank``, ``hgmd_class``:
  position-and-allele lookups (:mod:`paralog2viewer.lookups`).
* ``genes``: gene symbols used for the gene lookups (ANNOVAR refGene and
  knownGene genes, else the chosen VEP symbol); ``panel``,
  ``phenotypes``, ``inheritance``: PanelList/genemap2 strings.
* Chosen transcript: ``transcript_source`` (``default_transcript``,
  ``mane_select``, ``canonical_refseq``, ``canonical``, ``first_transcript``,
  ``first_row`` or ``none``), ``consequence``, ``gene`` (symbol),
  ``transcript``, ``exon``, ``hgvsc``, ``hgvsp``, ``am_class``,
  ``am_pathogenicity``, ``maxentscan`` (dict).

Coordinates are 1-based VCF coordinates of the normalized merged VCF.

Not verified (see the final report of the task): whether VEP echoes the
VCF ``ID`` in ``#Uploaded_variation`` (a Location+Allele fallback exists),
VEP ``--mane``, ANNOVAR ``-thread``/``--maxgenethread`` (ported from
the Nextflow version), the formats of the profile tables, and the
DefaultTranscript.tsv layout.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import logging
import re
import sqlite3
from collections.abc import Iterator, Mapping, Sequence
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path

import pysam

from paralog2viewer import __version__, proc
from paralog2viewer.datalayout import DataLayout, allele_key, bare_chrom
from paralog2viewer.lookups import GeneTables, LookupTables

LOGGER = logging.getLogger("paralog2viewer")

RECORD_SCHEMA_VERSION = "2"
# Cache key part: a code upgrade never reuses records of another version.
CACHE_SCHEMA = f"{RECORD_SCHEMA_VERSION}+{__version__}"
RECORD_FIELDS: tuple[str, ...] = (
    "allele_key",
    "chrom",
    "pos",
    "ref",
    "alt",
    "symbolic",
    "annovar",
    "intervar",
    "intervar_class",
    "vep_variant",
    "vep_transcripts",
    "clinvar",
    "clinvar_vid",
    "taiwan_biobank",
    "hgmd_class",
    "genes",
    "panel",
    "phenotypes",
    "inheritance",
    "transcript_source",
    "consequence",
    "gene",
    "transcript",
    "exon",
    "hgvsc",
    "hgvsp",
    "am_class",
    "am_pathogenicity",
    "maxentscan",
)

# Column of the VEP tab output that holds the VCF ID; a constant so the
# join can be adapted if VEP names it differently.
VEP_KEY_COLUMN = "#Uploaded_variation"
VEP_ASSEMBLY = "GRCh38"
VEP_CACHE_VERSION = "116"
# Added to the earlier VEP options. --mane is needed for MANE_SELECT in the
# output; its name in VEP 116 is not verified.
VEP_ADDED_OPTIONS = ("--hgvs", "--mane")
VEP_SYMBOLIC_OPTIONS = ("--canonical", "--numbers", "--mane")
VEP_TRANSCRIPT_KEYS = (
    "IMPACT",
    "SYMBOL",
    "BIOTYPE",
    "CANONICAL",
    "MANE_SELECT",
    "MANE_PLUS_CLINICAL",
    "EXON",
    "INTRON",
    "HGVSc",
    "HGVSp",
)
VEP_VARIANT_KEYS = (
    "VARIANT_CLASS",
    "MAX_AF",
    "MAX_AF_POPS",
    "AF",
    "CLIN_SIG",
    "SOMATIC",
    "PHENO",
    "PUBMED",
    "VAR_SYNONYMS",
)
VEP_PLUGIN_PREFIXES = ("MaxEntScan", "am_", "AlphaMissense")

ANNOVAR_PROTOCOLS = (
    "refGene,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_amr,"
    "1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas,exac03,"
    "esp6500siv2_all,esp6500siv2_aa,esp6500siv2_ea,nci60,avsnp147,cosmic70,"
    "clinvar_20240611,dbnsfp47a,gnomad_exome,gnomad_genome,dbscsnv11,rmsk,"
    "ensGene,knownGene,kaviar_20150923,hrcr1,regsnpintron,"
    "gene4denovo201907,icgc28,refGeneWithVer,gnomad41_exome,gnomad41_genome"
)
ANNOVAR_OPERATIONS = (
    "g,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,r,g,g,f,f,f,f,f,g,f,f"
)
# With -vcfinput -otherinfo, ANNOVAR writes zygosity, quality and depth as
# Otherinfo1-3, then the original VCF columns CHROM, POS, ID, REF, ALT, ...
OTHERINFO_VCF_START = 3
OTHERINFO_CHROM, OTHERINFO_POS, OTHERINFO_REF, OTHERINFO_ALT = 0, 1, 3, 4
INTERVAR_COLUMN_PREFIX = "InterVar and Evidence"
ENSEMBL_ID_RE = re.compile(r"ENS[A-Z]*[GT]\d")
INTERVAR_CLASS_RE = re.compile(r"InterVar:\s*(.*?)\s+PVS1=")
INTERVAR_CLASSES = (
    "Pathogenic",
    "Likely pathogenic",
    "Uncertain significance",
    "Likely benign",
    "Benign",
)
PLP_CLINVAR = frozenset(
    {"Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"}
)
PLP_INTERVAR = frozenset({"Pathogenic", "Likely pathogenic"})
SMALL_ALLELE_RE = re.compile(r"^[ACGTNacgtn]+$")
SYMBOLIC_INFO_KEYS = ("END", "SVTYPE", "SVLEN", "CIPOS", "CIEND", "MATEID")
REFSEQ_RE = re.compile(r"^(NM|NR|XM|XR)_")


# ANNOVAR ``-thread``: measured on one sample (4,111 alleles, one run each,
# whole step): 1 thread 49.9 s, 4 threads 52.0 s, 8 threads 69.9 s, 16
# threads 109.3 s. More threads were slower, and the multianno table is the
# same except the order of repeat names in the rmsk column, which varies from
# run to run.
ANNOVAR_THREADS = 1
# VEP ``--fork``. Larger values (16: 37 s, 32: 35 s against 41 s for 8) do not
# shorten the run because the ANNOVAR -> InterVar chain (about 65 s) is the
# longer branch, and they change a few VEP rows (HGNC ids differ between fork
# counts), so the value stays fixed whatever the budget.
VEP_FORK = 8


@dataclass(frozen=True)
class ToolPaths:
    """Locations of the annotators inside the container.

    Attributes:
        vep: VEP executable.
        perl: Perl executable used for ANNOVAR.
        annovar_dir: Directory with ``table_annovar.pl`` and friends.
        python: Python executable used for InterVar.
        intervar: ``Intervar.py`` path.
    """

    vep: str = "vep"
    perl: str = "perl"
    annovar_dir: str = "/opt/annovar"
    python: str = "python3"
    intervar: str = "/opt/intervar/Intervar.py"


@dataclass(frozen=True)
class Allele:
    """One unique allele of the merged VCF (1-based POS)."""

    chrom: str
    pos: int
    ref: str
    alt: str
    symbolic: bool
    info: str = "."

    @property
    def key(self) -> str:
        """The allele key (``chrom:pos:ref:alt``)."""
        return allele_key(self.chrom, self.pos, self.ref, self.alt)


@dataclass
class SitesVcfs:
    """Result of :func:`sites_only_vcf`.

    Attributes:
        small: Sites-only VCF of explicit sequence alleles, or None if empty.
        symbolic: Sites-only VCF of symbolic/breakend alleles, or None.
        alleles: All written alleles in output order.
    """

    small: Path | None
    symbolic: Path | None
    alleles: list[Allele] = field(default_factory=list)


def _info_string(rec: pysam.VariantRecord) -> str:
    """Return the symbolic-allele INFO fields kept for VEP."""
    parts = []
    for key in SYMBOLIC_INFO_KEYS:
        if key not in rec.info or key == "END":
            continue
        value = rec.info[key]
        if isinstance(value, tuple):
            value = ",".join(str(v) for v in value)
        parts.append(f"{key}={value}")
    if rec.alts and any(a.startswith("<") for a in rec.alts):
        parts.insert(0, f"END={rec.stop}")
    return ";".join(parts) or "."


def read_alleles(merged_vcf: Path) -> tuple[list[str], list[Allele]]:
    """Read unique alleles from a merged VCF.

    Multi-allelic records give one allele per ALT. Alleles repeated in
    several records, or present in several haplotype columns, are listed
    once (first occurrence). ``*`` and ``.`` ALT values are skipped.

    Args:
        merged_vcf: VCF or VCF.gz, haplotypes as sample columns.

    Returns:
        Header lines without FORMAT lines and sample names, and the alleles.
    """
    alleles: dict[str, Allele] = {}
    with pysam.VariantFile(str(merged_vcf)) as vcf:
        header = [
            line
            for line in str(vcf.header).splitlines()
            if line.startswith("##") and not line.startswith("##FORMAT")
        ]
        for rec in vcf:
            for alt in rec.alts or ():
                if alt in ("*", "."):
                    continue
                small = bool(
                    SMALL_ALLELE_RE.match(alt)
                    and SMALL_ALLELE_RE.match(rec.ref)
                )
                allele = Allele(
                    rec.chrom,
                    rec.pos,
                    rec.ref,
                    alt,
                    not small,
                    "." if small else _info_string(rec),
                )
                previous = alleles.setdefault(allele.key, allele)
                if previous.info != allele.info:
                    LOGGER.warning(
                        "allele key %s: differing INFO (%s vs %s); first kept",
                        allele.key,
                        previous.info,
                        allele.info,
                    )
    return header, list(alleles.values())


def _write_sites(
    path: Path, header: Sequence[str], rows: Sequence[Allele]
) -> None:
    """Write a sites-only VCF."""
    lines = list(header)
    have = "\n".join(lines)
    for key, number, kind, text in (
        ("END", "1", "Integer", "End position"),
        ("SVTYPE", "1", "String", "SV type"),
        ("SVLEN", ".", "Integer", "SV length"),
    ):
        if f"##INFO=<ID={key}," not in have:
            lines.append(
                f"##INFO=<ID={key},Number={number},Type={kind},"
                f'Description="{text}">'
            )
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    for allele in rows:
        lines.append(
            "\t".join(
                [
                    allele.chrom,
                    str(allele.pos),
                    allele.key,
                    allele.ref,
                    allele.alt,
                    ".",
                    ".",
                    allele.info,
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def sites_only_vcf(
    merged_vcf: Path, out_dir: Path, only: set[str] | None = None
) -> SitesVcfs:
    """Write sites-only VCFs, one line per unique allele.

    Args:
        merged_vcf: The merged VCF (``W/merged.vcf.gz``).
        out_dir: Output directory; ``small.vcf`` and ``symbolic.vcf`` are
            written there when they have alleles.
        only: If given, write only alleles with these keys.

    Returns:
        The paths and the alleles written.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    header, alleles = read_alleles(merged_vcf)
    if only is not None:
        alleles = [a for a in alleles if a.key in only]
    small = [a for a in alleles if not a.symbolic]
    symbolic = [a for a in alleles if a.symbolic]
    result = SitesVcfs(None, None, alleles)
    if small:
        result.small = out_dir / "small.vcf"
        _write_sites(result.small, header, small)
    if symbolic:
        result.symbolic = out_dir / "symbolic.vcf"
        _write_sites(result.symbolic, header, symbolic)
    return result


# --- cache -----------------------------------------------------------------


def manifest_checksum(layout: DataLayout) -> str:
    """Return the SHA-256 of ``manifest.json`` (the data version)."""
    return hashlib.sha256(layout.manifest.read_bytes()).hexdigest()


class RecordCache:
    """SQLite cache of finished records.

    Keyed by (allele key, data manifest checksum, schema version and tool
    version), so a data, schema or code change never returns stale records.
    """

    def __init__(self, path: Path, data_checksum: str) -> None:
        """Open or create the cache file.

        Args:
            path: SQLite file.
            data_checksum: SHA-256 of the data manifest.
        """
        path.parent.mkdir(parents=True, exist_ok=True)
        self._checksum = data_checksum
        self._db = sqlite3.connect(path)
        self._db.execute(
            "CREATE TABLE IF NOT EXISTS records ("
            "allele_key TEXT NOT NULL, data_checksum TEXT NOT NULL, "
            "schema_version TEXT NOT NULL, record TEXT NOT NULL, "
            "PRIMARY KEY (allele_key, data_checksum, schema_version))"
        )
        self._db.commit()

    def close(self) -> None:
        """Close the database."""
        self._db.close()

    def get_many(self, keys: Sequence[str]) -> dict[str, dict]:
        """Return the cached records among ``keys``."""
        found: dict[str, dict] = {}
        for start in range(0, len(keys), 500):
            chunk = keys[start : start + 500]
            marks = ",".join("?" * len(chunk))
            rows = self._db.execute(
                "SELECT allele_key, record FROM records WHERE "
                f"data_checksum=? AND schema_version=? "
                f"AND allele_key IN ({marks})",
                [self._checksum, CACHE_SCHEMA, *chunk],
            )
            found.update({k: json.loads(r) for k, r in rows})
        return found

    def put_many(self, records: Sequence[Mapping]) -> None:
        """Store records (replacing existing ones)."""
        self._db.executemany(
            "INSERT OR REPLACE INTO records VALUES (?, ?, ?, ?)",
            [
                (
                    r["allele_key"],
                    self._checksum,
                    CACHE_SCHEMA,
                    json.dumps(r),
                )
                for r in records
            ],
        )
        self._db.commit()


# --- annotator commands ----------------------------------------------------


def vep_command(
    layout: DataLayout,
    reference: Path,
    vcf: Path,
    output: Path,
    threads: int,
    symbolic: bool,
    tools: ToolPaths = ToolPaths(),
) -> list[str]:
    """Return the VEP command line (ported, plus --hgvs).

    Args:
        layout: Data directory (cache, plugins, plugin data).
        reference: Reference FASTA.
        vcf: Sites-only VCF.
        output: Output tab file.
        threads: VEP ``--fork`` (omitted when 1).
        symbolic: Use the short symbolic/breakend option set.
        tools: Executable locations.
    """
    cmd = [
        tools.vep,
        "--input_file", str(vcf),
        "--output_file", str(output),
        "--format", "vcf",
    ]  # fmt: skip
    if symbolic:
        cmd += ["--symbol", *VEP_SYMBOLIC_OPTIONS]
    else:
        cmd += [
            "--sift", "b", "--polyphen", "b", "--ccds", "--symbol",
            "--numbers", "--domains", "--regulatory", "--canonical",
            "--protein", "--biotype", "--uniprot", "--tsl", "--appris",
            "--gene_phenotype", "--af", "--af_1kg", "--af_gnomad",
            "--max_af", "--pubmed", "--var_synonyms", "--variant_class",
            *VEP_ADDED_OPTIONS,
        ]  # fmt: skip
    cmd += ["--assembly", VEP_ASSEMBLY, "--force_overwrite"]
    if threads > 1:
        cmd += ["--fork", str(threads)]
    if not symbolic:
        cmd += ["--use_given_ref"]
    cmd += [
        "--dir", str(layout.vep_cache),
        "--offline",
        "--fasta", str(reference),
        "--merged",
    ]  # fmt: skip
    if not symbolic:
        cmd += [
            "--buffer_size", "262144",
            "--dir_plugins", str(layout.vep_plugins),
            "--plugin", f"MaxEntScan,{layout.maxentscan}",
            "--plugin", f"AlphaMissense,file={layout.alphamissense}",
        ]  # fmt: skip
    cmd += ["--cache_version", VEP_CACHE_VERSION, "--tab"]
    return cmd


def annovar_command(
    layout: DataLayout,
    vcf: Path,
    out_prefix: Path,
    threads: int,
    tools: ToolPaths = ToolPaths(),
) -> list[str]:
    """Return the ANNOVAR ``table_annovar.pl`` command.

    ``-thread``/``--maxgenethread`` are kept from the Nextflow version and are
    not independently verified against the ANNOVAR documentation.
    """
    return [
        tools.perl, f"{tools.annovar_dir}/table_annovar.pl", str(vcf),
        str(layout.annovar_humandb),
        "-buildver", "hg38", "-out", str(out_prefix),
        "-protocol", ANNOVAR_PROTOCOLS,
        "-operation", ANNOVAR_OPERATIONS,
        "-vcfinput", "-otherinfo",
        "-thread", str(threads), "--maxgenethread", str(threads),
        "-mingenelinecount", "1",
    ]  # fmt: skip


def intervar_command(
    layout: DataLayout,
    avinput: Path,
    out_prefix: Path,
    tools: ToolPaths = ToolPaths(),
) -> list[str]:
    """Return the InterVar command (``--skip_annovar``)."""
    d = tools.annovar_dir
    return [
        tools.python, tools.intervar,
        f"--input={avinput}", "--input_type=AVinput",
        f"--output={out_prefix}", "--buildver=hg38",
        f"--database_intervar={layout.intervardb}",
        f"--table_annovar={d}/table_annovar.pl",
        f"--convert2annovar={d}/convert2annovar.pl",
        f"--annotate_variation={d}/annotate_variation.pl",
        f"--database_locat={layout.annovar_humandb}",
        "--skip_annovar",
    ]  # fmt: skip


# --- parsing and joining ---------------------------------------------------


def read_vep_tab(path: Path) -> list[dict[str, str]]:
    """Read a VEP ``--tab`` file into row dicts.

    Lines starting ``##`` are skipped; the line starting with ``#`` is the
    header. The ``Extra`` column is split into ``KEY=VALUE`` pairs stored
    under the key ``"_extra"`` (a dict).
    """
    rows: list[dict[str, str]] = []
    header: list[str] | None = None
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("##") or not line:
                continue
            if line.startswith("#"):
                header = line.split("\t")
                continue
            if header is None:
                raise ValueError(f"{path}: data before VEP header")
            row = dict(zip(header, line.split("\t"), strict=False))
            extra = {}
            for item in row.get("Extra", "").split(";"):
                name, _, value = item.partition("=")
                if name:
                    extra[name] = value
            row["_extra"] = extra  # type: ignore[assignment]
            rows.append(row)
    return rows


def vep_location_allele(allele: Allele) -> tuple[str, int, int, str]:
    """Return how VEP reports an allele: (bare chrom, lo, hi, allele).

    VEP drops the first base when REF and ALT share it, writes ``-`` for an
    empty allele, and reports an insertion with start after end; ``lo`` and
    ``hi`` are sorted so the order does not matter.
    """
    ref, alt, start = allele.ref, allele.alt, allele.pos
    if ref[:1] == alt[:1]:
        ref, alt, start = ref[1:], alt[1:], start + 1
    end = start + len(ref) - 1
    lo, hi = sorted((start, end))
    return bare_chrom(allele.chrom), lo, hi, alt or "-"


def _row_location_allele(row: Mapping) -> tuple[str, int, int, str] | None:
    """Return (bare chrom, lo, hi, allele) of a VEP row, or None."""
    chrom, _, span = row.get("Location", "").partition(":")
    numbers = [int(n) for n in re.findall(r"\d+", span)]
    if not chrom or not numbers:
        return None
    return bare_chrom(chrom), min(numbers), max(numbers), row.get("Allele", "")


def group_vep_rows(
    rows: Sequence[dict], alleles: Sequence[Allele]
) -> dict[str, list[dict]]:
    """Assign VEP rows to allele keys.

    Primary join: ``#Uploaded_variation`` equals the allele key. Fallback
    for rows whose ID is not a known key: Location and Allele with VEP's
    indel trimming.
    """
    known = {a.key for a in alleles}
    by_loc: dict[tuple[str, int, int, str], list[str]] = {}
    for a in alleles:
        by_loc.setdefault(vep_location_allele(a), []).append(a.key)
    grouped: dict[str, list[dict]] = {}
    for row in rows:
        ident = row.get(VEP_KEY_COLUMN, "")
        if ident in known:
            grouped.setdefault(ident, []).append(row)
            continue
        loc = _row_location_allele(row)
        keys = by_loc.get(loc, []) if loc else []
        if len(keys) == 1:
            grouped.setdefault(keys[0], []).append(row)
        else:
            LOGGER.warning("VEP row not joined to an allele: %s", ident)
    return grouped


def _dash(value: str | None) -> str:
    """Map VEP's ``-`` placeholder to an empty string."""
    return "" if value in (None, "-") else str(value)


def vep_transcripts(rows: Sequence[Mapping]) -> list[dict]:
    """Convert VEP rows of one allele to transcript dicts (see module doc)."""
    out = []
    for row in rows:
        extra = row.get("_extra", {})
        out.append(
            {
                "feature": _dash(row.get("Feature")),
                "feature_type": _dash(row.get("Feature_type")),
                "gene_id": _dash(row.get("Gene")),
                "symbol": _dash(extra.get("SYMBOL")),
                "consequence": _dash(row.get("Consequence")),
                "impact": _dash(extra.get("IMPACT")),
                "biotype": _dash(extra.get("BIOTYPE")),
                "canonical": _dash(extra.get("CANONICAL")),
                "mane_select": _dash(extra.get("MANE_SELECT")),
                "mane_plus_clinical": _dash(extra.get("MANE_PLUS_CLINICAL")),
                "exon": _dash(extra.get("EXON")),
                "intron": _dash(extra.get("INTRON")),
                "hgvsc": _dash(extra.get("HGVSc")),
                "hgvsp": _dash(extra.get("HGVSp")),
                "cdna_position": _dash(row.get("cDNA_position")),
                "cds_position": _dash(row.get("CDS_position")),
                "protein_position": _dash(row.get("Protein_position")),
                "amino_acids": _dash(row.get("Amino_acids")),
                "codons": _dash(row.get("Codons")),
                "plugins": {
                    k: v
                    for k, v in extra.items()
                    if k.startswith(VEP_PLUGIN_PREFIXES) and v not in ("", "-")
                },
            }
        )
    return out


def vep_variant_fields(rows: Sequence[Mapping]) -> dict[str, str]:
    """Return allele-level VEP fields (first non-empty value over rows)."""
    out: dict[str, str] = {}
    for row in rows:
        existing = _dash(row.get("Existing_variation"))
        if existing and not out.get("Existing_variation"):
            out["Existing_variation"] = existing
        for key, value in row.get("_extra", {}).items():
            wanted = key in VEP_VARIANT_KEYS or key.endswith("_AF")
            if wanted and value and key not in out:
                out[key] = value
    return out


def _unversioned(accession: str) -> str:
    """Strip the version suffix from a transcript accession."""
    return accession.split(".")[0]


def choose_transcript(
    transcripts: Sequence[dict], default_accessions: set[str]
) -> tuple[dict | None, str]:
    """Pick the transcript whose fields fill the table.

    Order: a transcript listed in DefaultTranscript.tsv (version ignored);
    the MANE Select transcript (the row whose accession equals a row's
    ``MANE_SELECT`` value, else any row carrying that value); the canonical
    RefSeq transcript; any canonical transcript; the first transcript row;
    the first row of any kind.

    Returns:
        The chosen dict (or None without rows) and the source label.
    """
    rows = [t for t in transcripts if t["feature_type"] == "Transcript"]
    for t in rows:
        if _unversioned(t["feature"]) in default_accessions:
            return t, "default_transcript"
    mane = {_unversioned(t["mane_select"]) for t in rows if t["mane_select"]}
    for t in rows:
        if _unversioned(t["feature"]) in mane:
            return t, "mane_select"
    for t in rows:
        if t["mane_select"]:
            return t, "mane_select"
    for t in rows:
        if t["canonical"] == "YES" and REFSEQ_RE.match(t["feature"]):
            return t, "canonical_refseq"
    for t in rows:
        if t["canonical"] == "YES":
            return t, "canonical"
    if rows:
        return rows[0], "first_transcript"
    if transcripts:
        return transcripts[0], "first_row"
    return None, "none"


def read_multianno(path: Path) -> dict[str, dict[str, str]]:
    """Read ANNOVAR's multianno table keyed by the original VCF allele.

    The key is built from the VCF chrom, pos, ref and alt that
    ``-vcfinput -otherinfo`` keeps after the zygosity/quality/depth columns,
    because ANNOVAR rewrites indel coordinates in its own columns.

    Returns:
        Allele key -> column dict without ``Otherinfo*`` columns.
    """
    out: dict[str, dict[str, str]] = {}
    with open(path, encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        first_other = next(
            (i for i, n in enumerate(header) if n.startswith("Otherinfo")),
            None,
        )
        if first_other is None:
            raise ValueError(f"{path}: no Otherinfo columns (-otherinfo)")
        start = first_other + OTHERINFO_VCF_START
        for line in handle:
            cells = line.rstrip("\n").split("\t")
            if len(cells) <= start + OTHERINFO_ALT:
                raise ValueError(f"{path}: row lacks VCF columns: {line[:80]}")
            key = allele_key(
                cells[start + OTHERINFO_CHROM],
                int(cells[start + OTHERINFO_POS]),
                cells[start + OTHERINFO_REF],
                cells[start + OTHERINFO_ALT],
            )
            out[key] = dict(
                zip(header[:first_other], cells[:first_other], strict=True)
            )
    return out


def _annovar_id(cols: Mapping[str, str]) -> tuple[str, ...]:
    """Return the (Chr, Start, End, Ref, Alt) join key of a row."""
    return (
        bare_chrom(cols["Chr"]),
        cols["Start"],
        cols["End"],
        cols["Ref"],
        cols["Alt"],
    )


def read_intervar(path: Path) -> dict[tuple[str, ...], str]:
    """Read InterVar's ``.intervar`` file.

    InterVar writes one line per gene for an intergenic variant; the first
    line is kept, as in the baseline (verified on all 40 rows that differ).

    Returns:
        (bare Chr, Start, End, Ref, Alt) -> the raw ``InterVar and
        Evidence`` text.
    """
    out: dict[tuple[str, ...], str] = {}
    with open(path, encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").lstrip("#").split("\t")
        column = next(
            (i for i, n in enumerate(header) if INTERVAR_COLUMN_PREFIX in n),
            None,
        )
        if column is None:
            raise ValueError(f"{path}: no {INTERVAR_COLUMN_PREFIX!r} column")
        for line in handle:
            cells = line.rstrip("\n").split("\t")
            if len(cells) <= column:
                continue
            out.setdefault((bare_chrom(cells[0]), *cells[1:5]), cells[column])
    return out


def intervar_class(text: str) -> str:
    """Return the class in an InterVar summary string, or ``""``."""
    match = INTERVAR_CLASS_RE.search(text)
    if match and match.group(1) in INTERVAR_CLASSES:
        return match.group(1)
    return ""


def _annovar_genes(annovar: Mapping[str, str]) -> list[str]:
    """Return gene symbols from ANNOVAR's refGene and knownGene columns.

    This is the gene set the baseline Panel, Phenotypes and Inheritance
    lookups use (verified on all baseline rows): ``Gene.refGene`` then
    ``Gene.knownGene``, distance annotations removed, Ensembl transcript or
    gene IDs dropped, duplicates removed. ``Gene.ensGene`` is not used.
    """
    lower = {name.lower(): value for name, value in annovar.items()}
    symbols: list[str] = []
    for column in ("gene.refgene", "gene.knowngene"):
        text = re.sub(r"\([^)]*\)", "", lower.get(column, "."))
        symbols += [
            s
            for s in re.split(r"[;,]", text)
            if s and s != "." and not ENSEMBL_ID_RE.match(s)
        ]
    return list(dict.fromkeys(symbols))


def build_record(
    allele: Allele,
    vep_rows: Sequence[Mapping],
    annovar: Mapping[str, str] | None,
    intervar: str,
    tables: LookupTables,
) -> dict:
    """Join all sources for one allele into a record.

    Args:
        allele: The allele.
        vep_rows: VEP rows of this allele (from :func:`read_vep_tab`).
        annovar: ANNOVAR columns, or None (symbolic or dropped allele).
        intervar: Raw InterVar text, or ``""``.
        tables: Open lookup tables.
    """
    transcripts = vep_transcripts(vep_rows)
    gene_tables: GeneTables = tables.genes
    chosen, source = choose_transcript(
        transcripts, gene_tables.all_default_accessions
    )
    chosen = chosen or {}
    plugins = chosen.get("plugins", {})
    symbols = _annovar_genes(annovar or {})
    if not symbols and chosen.get("symbol"):
        symbols = [chosen["symbol"]]
    lookups = tables.variant_fields(
        allele.chrom, allele.pos, allele.ref, allele.alt
    )
    gene_info = gene_tables.gene_info(symbols)
    record = {
        "allele_key": allele.key,
        "chrom": allele.chrom,
        "pos": allele.pos,
        "ref": allele.ref,
        "alt": allele.alt,
        "symbolic": allele.symbolic,
        "annovar": dict(annovar or {}),
        "intervar": intervar,
        "intervar_class": intervar_class(intervar),
        "vep_variant": vep_variant_fields(vep_rows),
        "vep_transcripts": transcripts,
        **lookups,
        "genes": symbols,
        **gene_info,
        "transcript_source": source,
        "consequence": chosen.get("consequence", ""),
        "gene": chosen.get("symbol", ""),
        "transcript": chosen.get("feature", ""),
        "exon": chosen.get("exon", ""),
        "hgvsc": chosen.get("hgvsc", ""),
        "hgvsp": chosen.get("hgvsp", ""),
        "am_class": plugins.get("am_class", ""),
        "am_pathogenicity": plugins.get("am_pathogenicity", ""),
        "maxentscan": {
            k: v for k, v in plugins.items() if k.startswith("MaxEntScan")
        },
    }
    return {name: record[name] for name in RECORD_FIELDS}


# --- orchestration ---------------------------------------------------------


def _run_env(directory: Path) -> dict[str, str]:
    """Return HOME/cache variables that keep tool state in ``directory``."""
    return {"HOME": str(directory), "XDG_CACHE_HOME": str(directory / ".cache")}


def _vep_branch(
    layout: DataLayout,
    reference: Path,
    sites: SitesVcfs,
    ann_dir: Path,
    threads: int,
    tools: ToolPaths,
) -> None:
    """Run VEP on the small and then the symbolic sites-only VCF."""
    env = _run_env(ann_dir)
    for vcf, name, symbolic in (
        (sites.small, "small.vep.txt", False),
        (sites.symbolic, "symbolic.vep.txt", True),
    ):
        if vcf is not None:
            cmd = vep_command(
                layout, reference, vcf, ann_dir / name, threads, symbolic, tools
            )
            proc.run(cmd, cwd=ann_dir, env=env)


def _annovar_branch(
    layout: DataLayout,
    sites: SitesVcfs,
    ann_dir: Path,
    tools: ToolPaths,
) -> None:
    """Run ANNOVAR then InterVar on the small-variant VCF."""
    if sites.small is None:
        return
    env = _run_env(ann_dir)
    prefix = ann_dir / "small"
    proc.run(
        annovar_command(layout, sites.small, prefix, ANNOVAR_THREADS, tools),
        cwd=ann_dir,
        env=env,
    )
    proc.run(
        intervar_command(layout, ann_dir / "small.avinput", prefix, tools),
        cwd=ann_dir,
        env=env,
    )
    if not (ann_dir / "small.hg38_multianno.txt.intervar").stat().st_size:
        raise proc.ToolError("InterVar wrote an empty output")


def run_annotators(
    layout: DataLayout,
    reference: Path,
    sites: SitesVcfs,
    ann_dir: Path,
    threads: int,
    tools: ToolPaths = ToolPaths(),
) -> None:
    """Run VEP and ANNOVAR->InterVar concurrently in two threads.

    Raises:
        ToolError: If any command fails.
    """
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures = [
            pool.submit(
                _vep_branch, layout, reference, sites, ann_dir, threads, tools
            ),
            pool.submit(_annovar_branch, layout, sites, ann_dir, tools),
        ]
        for future in futures:
            future.result()


def annotate_new_alleles(
    layout: DataLayout,
    reference: Path,
    merged_vcf: Path,
    keys: set[str],
    ann_dir: Path,
    threads: int,
    tools: ToolPaths = ToolPaths(),
) -> list[dict]:
    """Annotate the alleles in ``keys`` and return their records."""
    sites = sites_only_vcf(merged_vcf, ann_dir, only=keys)
    run_annotators(layout, reference, sites, ann_dir, threads, tools)
    vep_rows: list[dict] = []
    for name in ("small.vep.txt", "symbolic.vep.txt"):
        if (ann_dir / name).exists():
            vep_rows += read_vep_tab(ann_dir / name)
    grouped = group_vep_rows(vep_rows, sites.alleles)
    annovar: dict[str, dict[str, str]] = {}
    intervar: dict[tuple[str, ...], str] = {}
    if sites.small is not None:
        annovar = read_multianno(ann_dir / "small.hg38_multianno.txt")
        intervar = read_intervar(ann_dir / "small.hg38_multianno.txt.intervar")
    records = []
    with LookupTables(layout) as tables:
        for allele in sites.alleles:
            cols = annovar.get(allele.key)
            text = intervar.get(_annovar_id(cols), "") if cols else ""
            records.append(
                build_record(
                    allele, grouped.get(allele.key, []), cols, text, tables
                )
            )
    return records


def annotate_alleles(
    merged_vcf: Path,
    reference: Path,
    data_layout: DataLayout,
    work_dir: Path,
    cache_path: Path,
    threads: int,
    tools: ToolPaths = ToolPaths(),
) -> Path:
    """Annotate every unique allele of the merged VCF.

    Args:
        merged_vcf: ``W/merged.vcf.gz``.
        reference: Reference FASTA (matches the data directory build).
        data_layout: Data directory layout; its manifest keys the cache.
        work_dir: Run work directory ``W``; files go to ``W/annotation``.
        cache_path: SQLite cache file (may be shared across runs).
        threads: VEP ``--fork`` (ANNOVAR uses ``ANNOVAR_THREADS``).
        tools: Annotator locations.

    Returns:
        ``W/annotation/records.jsonl.gz``, one JSON record per allele.
    """
    ann_dir = work_dir / "annotation"
    ann_dir.mkdir(parents=True, exist_ok=True)
    _, alleles = read_alleles(merged_vcf)
    order = [a.key for a in alleles]
    cache = RecordCache(cache_path, manifest_checksum(data_layout))
    try:
        records = cache.get_many(order)
        missing = {k for k in order if k not in records}
        LOGGER.info(
            "annotation: %d alleles, %d cached, %d to annotate",
            len(order),
            len(records),
            len(missing),
        )
        if missing:
            fresh = annotate_new_alleles(
                data_layout,
                reference,
                merged_vcf,
                missing,
                ann_dir,
                threads,
                tools,
            )
            cache.put_many(fresh)
            records.update({r["allele_key"]: r for r in fresh})
    finally:
        cache.close()
    out = ann_dir / "records.jsonl.gz"
    with open(out, "wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as gz:
            for key in order:
                line = json.dumps(records[key]) + "\n"
                gz.write(line.encode("utf-8"))
    return out


def read_records(path: Path) -> Iterator[dict]:
    """Yield records from ``records.jsonl.gz``."""
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            yield json.loads(line)


def is_pathogenic(record: Mapping) -> bool:
    """Return True if ClinVar or InterVar says Pathogenic/Likely pathogenic.

    Exact match of the classification: ANNOVAR ``ClinVar_SIG`` against
    {Pathogenic, Likely_pathogenic, Pathogenic/Likely_pathogenic}, the
    ClinVar table value (spaces read as underscores) against the same set,
    and the InterVar class against {Pathogenic, Likely pathogenic}.
    """
    values = [record.get("clinvar", "").replace(" ", "_")]
    values += [
        v
        for k, v in record.get("annovar", {}).items()
        if k.lower() == "clinvar_sig"
    ]
    if any(v in PLP_CLINVAR for v in values):
        return True
    return record.get("intervar_class", "") in PLP_INTERVAR


def pathogenic_allele_keys(records_path: Path) -> set[str]:
    """Return keys of alleles classified Pathogenic or Likely pathogenic.

    Used for the white-diamond markers in the figures.
    """
    return {
        r["allele_key"] for r in read_records(records_path) if is_pathogenic(r)
    }

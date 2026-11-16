# SPDX-License-Identifier: GPL-3.0-or-later
"""Table lookups that replace the MViewer.CLI profile tables.

Variant lookups (ClinVar, Taiwan Biobank, HGMD) use tabix on the tables in
the data directory (see :class:`paralog2viewer.datalayout.DataLayout`) and
match on position and alleles (Taiwan Biobank: exact position and REF).
Gene lookups read ``PanelList.tsv``,
``genemap2.txt`` (OMIM) and ``DefaultTranscript.tsv`` into memory once.

Coordinates: all inputs are 1-based VCF-style (POS, REF, ALT). Table rows
written in ANNOVAR style (``-`` for an empty allele, start of the base before
an insertion) are converted to the same minimal form before comparison, so
VCF-style and ANNOVAR-style rows both match (:func:`minimal_allele`).

Validated against the baseline MViewer tables (two samples, 9,901 rows); the
rules are in :meth:`AlleleFreqTable.query` and :meth:`GeneTables.gene_info`.
Column names are found by candidate lists; a table whose header matches none
raises ``LookupSchemaError`` rather than returning silent blanks. Callers
can pass their own candidate lists. The baseline ``Gene`` column is a
multi-source union of symbols; its sources are not resolved here, the caller
supplies the symbols.
"""

from __future__ import annotations

import gzip
import re
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path

import pysam

from paralog2viewer.datalayout import DataLayout, bare_chrom

# Candidate header names (compared case-insensitively, first hit wins).
CLINVAR_COLUMNS: Mapping[str, Sequence[str]] = {
    "ref": ("ref", "reference", "reference_allele", "referenceallelevcf"),
    "alt": ("alt", "alternate", "alternate_allele", "alternateallelevcf"),
    "assembly": ("assembly", "genome_build", "build"),
    "clinvar": (
        "clinvar",
        "clnsig",
        "clinical_significance",
        "clinicalsignificance",
        "significance",
        "clinvar_sig",
    ),
    "clinvar_vid": (
        "clinvar_vid",
        "vid",
        "variationid",
        "variation_id",
        "variation id",
    ),
}
# The real table has one row per site, ``Chr Pos Ref Freq.``; ``Freq.`` is a
# space-separated list of ``allele:frequency`` pairs that includes the
# reference allele (alleles are VCF-style anchored; ``*`` may occur).
TAIWAN_BIOBANK_COLUMNS: Mapping[str, Sequence[str]] = {
    "ref": ("ref", "reference", "reference_allele"),
    "freqs": ("freq.", "freq", "allele_freqs"),
}
# Values of an ``assembly`` column that denote GRCh38 (lower case).
GRCH38_NAMES = frozenset({"grch38", "hg38"})
PANEL_COLUMNS: Mapping[str, Sequence[str]] = {
    "gene": ("gene", "genes", "gene symbol", "gene_symbol", "symbol"),
    "panel": ("panel", "panels", "panel name", "panel_name", "list"),
}
DEFAULT_TRANSCRIPT_RE = re.compile(r"\b((?:NM|NR|XM|XR)_\d+|ENST\d+)(?:\.\d+)?")
HGMD_CLASS_KEY = "CLASS"
OMIM_SYMBOL_COLUMN = "Approved Gene Symbol"
OMIM_PHENOTYPE_COLUMN = "Phenotypes"
OMIM_ALIASES_COLUMN = "Gene Symbols"
# OMIM inheritance term -> short code (assumed to follow MViewer's display;
# unknown terms are kept unchanged).
INHERITANCE_CODES: Mapping[str, str] = {
    "autosomal dominant": "AD",
    "autosomal recessive": "AR",
    "x-linked": "XL",
    "x-linked dominant": "XLD",
    "x-linked recessive": "XLR",
    "y-linked": "YL",
    "mitochondrial": "Mi",
    "multifactorial": "Multi",
    "isolated cases": "IC",
    "digenic dominant": "DD",
    "digenic recessive": "DR",
    "somatic mutation": "SMu",
    "somatic mosaicism": "SMo",
    "pseudoautosomal dominant": "PAD",
    "pseudoautosomal recessive": "PAR",
}
# Display order of inheritance codes in the baseline. Only these orders are
# seen there (AD, AR, Multi and XLR, XLD, XL); the position of Multi relative
# to the X-linked codes and the order of all other codes are not constrained
# by the data, so other codes follow alphabetically.
INHERITANCE_ORDER = ("AD", "AR", "Multi", "XLR", "XLD", "XL")
_PHENOTYPE_RE = re.compile(
    r"^(?P<name>.*?),?\s*(?P<mim>\d{6})?\s*\((?P<key>[1-4])\)"
    r"(?:,\s*(?P<inh>.*))?$"
)


class LookupSchemaError(ValueError):
    """A lookup table header does not contain the expected columns."""


def minimal_allele(pos: int, ref: str, alt: str) -> tuple[int, str, str]:
    """Return the minimal representation of an allele.

    Common trailing bases are removed first, then common leading bases (the
    position moves right by the number of leading bases removed). No left
    alignment is done. ANNOVAR-style ``-`` alleles are first converted: an
    empty REF (insertion) moves the position one base right, because ANNOVAR
    stores the base before the insertion.

    Args:
        pos: 1-based position of the first REF base.
        ref: Reference allele, or ``-``.
        alt: Alternate allele, or ``-``.

    Returns:
        Tuple of 1-based position, REF and ALT; REF or ALT may be empty.
    """
    if ref == "-":
        # ANNOVAR insertion: POS is the base before the inserted bases.
        return pos + 1, "", alt.upper()
    ref = ref.upper()
    alt = "" if alt == "-" else alt.upper()
    while ref and alt and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    shift = 0
    while ref and alt and ref[0] == alt[0]:
        ref, alt = ref[1:], alt[1:]
        shift += 1
    return pos + shift, ref, alt


def _candidate_contigs(chrom: str) -> list[str]:
    """Return contig spellings to try for ``chrom`` (with/without chr)."""
    bare = bare_chrom(chrom)
    names = [chrom, bare, f"chr{bare}"]
    return list(dict.fromkeys(names))


def _find_columns(
    names: Sequence[str],
    wanted: Mapping[str, Sequence[str]],
    path: Path,
    optional: Iterable[str] = (),
) -> dict[str, int]:
    """Map logical names to column indexes by candidate header names."""
    lowered = [n.strip().lower() for n in names]
    found: dict[str, int] = {}
    for logical, candidates in wanted.items():
        for candidate in candidates:
            if candidate in lowered:
                found[logical] = lowered.index(candidate)
                break
        else:
            if logical not in optional:
                raise LookupSchemaError(
                    f"{path}: no column for {logical!r}; header is {names}"
                )
    return found


def _read_header(path: Path) -> list[str]:
    """Return the first line of a gzip text file split on tabs."""
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        line = handle.readline().rstrip("\n")
    return line.lstrip("#").split("\t")


class VariantTable:
    """A tabix-indexed TSV (header in line 1) queried by position/alleles.

    Attributes:
        path: The bgzip file; ``path + '.tbi'`` must exist.
    """

    def __init__(
        self, path: Path, columns: Mapping[str, Sequence[str]]
    ) -> None:
        """Open the table.

        Args:
            path: Bgzip TSV with chromosome in column 1, 1-based position in
                column 2 and a header line.
            columns: Logical name -> candidate header names; must include
                ``ref`` and ``alt``. Other names are returned as values,
                except the optional ``assembly``: if the table has that
                column, rows whose value is not GRCh38 are ignored.

        Raises:
            LookupSchemaError: If a column cannot be found.
        """
        self.path = path
        names = _read_header(path)
        self._index = _find_columns(names, columns, path, ("assembly",))
        self._tabix = pysam.TabixFile(str(path))

    def close(self) -> None:
        """Close the file."""
        self._tabix.close()

    def query(
        self, chrom: str, pos: int, ref: str, alt: str
    ) -> list[dict[str, str]]:
        """Return rows whose position and alleles equal the query.

        Args:
            chrom: Contig name; ``chr`` prefix is tried both ways.
            pos: 1-based VCF position.
            ref: VCF REF.
            alt: VCF ALT.

        Returns:
            One dict per matching row with the value columns as strings.
        """
        target = minimal_allele(pos, ref, alt)
        window = (max(0, pos - 2), pos + len(ref) + 2)
        rows: list[dict[str, str]] = []
        for contig in _candidate_contigs(chrom):
            try:
                lines = list(self._tabix.fetch(contig, *window))
            except ValueError:
                continue
            for line in lines:
                cells = line.split("\t")
                if "assembly" in self._index and (
                    cells[self._index["assembly"]].strip().lower()
                    not in GRCH38_NAMES
                ):
                    continue
                row_pos = int(cells[1])
                row_ref = cells[self._index["ref"]]
                row_alt = cells[self._index["alt"]]
                if minimal_allele(row_pos, row_ref, row_alt) != target:
                    continue
                rows.append(
                    {
                        name: cells[idx]
                        for name, idx in self._index.items()
                        if name not in ("ref", "alt", "assembly")
                    }
                )
            break
        return rows


class AlleleFreqTable:
    """A tabix-indexed per-site allele-frequency TSV (Taiwan Biobank).

    Each row is one site (``Chr Pos Ref Freq.``) and ``Freq.`` lists
    ``allele:frequency`` pairs including the reference allele.
    """

    def __init__(
        self, path: Path, columns: Mapping[str, Sequence[str]]
    ) -> None:
        """Open the table.

        Args:
            path: Bgzip TSV with chromosome in column 1, 1-based position in
                column 2 and a header line.
            columns: Logical names ``ref`` and ``freqs`` -> candidate header
                names.

        Raises:
            LookupSchemaError: If a column cannot be found.
        """
        self.path = path
        self._index = _find_columns(_read_header(path), columns, path)
        self._tabix = pysam.TabixFile(str(path))

    def close(self) -> None:
        """Close the file."""
        self._tabix.close()

    def query(self, chrom: str, pos: int, ref: str, alt: str) -> str:
        """Return the Taiwan Biobank value the way MViewer reports it.

        The row must have the same position and the same REF as the query
        (no normalisation: against the baseline, a differently anchored row
        never matched). If the ALT is one of the row's alleles, the value is
        its frequency with five decimals (``0.0010`` becomes ``0.00100``).
        If the row exists but lacks the ALT, MViewer reports the row's whole
        ``Freq.`` text; this is kept so that the table equals the baseline.

        Args:
            chrom: Contig name; ``chr`` prefix is tried both ways.
            pos: 1-based VCF position.
            ref: VCF REF.
            alt: VCF ALT.

        Returns:
            The frequency, the whole ``Freq.`` text, or ``""`` without a row.
        """
        ref, alt = ref.upper(), alt.upper()
        for contig in _candidate_contigs(chrom):
            try:
                lines = list(self._tabix.fetch(contig, pos - 1, pos))
            except ValueError:
                continue
            for line in lines:
                cells = line.split("\t")
                if (
                    int(cells[1]) != pos
                    or cells[self._index["ref"]].upper() != ref
                ):
                    continue
                text = cells[self._index["freqs"]]
                for pair in text.split():
                    allele, _, freq = pair.rpartition(":")
                    if allele.upper() == alt:
                        return f"{float(freq):.5f}"
                return text
            break
        return ""


def _inheritance_sort_key(code: str) -> tuple[int, str]:
    """Sort key placing the codes of INHERITANCE_ORDER first, in order."""
    if code in INHERITANCE_ORDER:
        return INHERITANCE_ORDER.index(code), code
    return len(INHERITANCE_ORDER), code


def _inheritance_codes(text: str) -> list[str]:
    """Convert an OMIM inheritance string to short codes."""
    codes = []
    for part in text.split(","):
        term = part.strip()
        if term:
            codes.append(INHERITANCE_CODES.get(term.lower(), term))
    return codes


def parse_omim_phenotypes(text: str) -> list[dict[str, str]]:
    """Parse the ``Phenotypes`` cell of genemap2.

    Args:
        text: For example ``Name, 123456 (3), Autosomal dominant; ...``.

    Returns:
        One dict per phenotype with ``text`` (the entry as written),
        ``name``, ``mim``, ``mapping_key`` and ``inheritance`` (short codes
        joined by ``, ``); entries that do not parse are kept with only
        ``text`` and ``name``.
    """
    out = []
    for item in text.split(";"):
        item = item.strip()
        if not item:
            continue
        match = _PHENOTYPE_RE.match(item)
        if match is None:
            out.append(
                {
                    "text": item,
                    "name": item,
                    "mim": "",
                    "mapping_key": "",
                    "inheritance": "",
                }
            )
            continue
        out.append(
            {
                "text": item,
                "name": match["name"].strip().rstrip(","),
                "mim": match["mim"] or "",
                "mapping_key": match["key"],
                "inheritance": ", ".join(
                    _inheritance_codes(match["inh"] or "")
                ),
            }
        )
    return out


@dataclass
class GeneTables:
    """Gene-level tables held in memory.

    Attributes:
        panels: Gene symbol -> panel names in file order.
        omim: Gene symbol (approved, else alias) -> phenotype dicts
            (see parse_omim_phenotypes).
        default_transcripts: Gene symbol -> unversioned accessions, or the
            key ``""`` when a line has no recognisable gene column.
    """

    panels: dict[str, list[str]] = field(default_factory=dict)
    omim: dict[str, list[dict[str, str]]] = field(default_factory=dict)
    default_transcripts: dict[str, set[str]] = field(default_factory=dict)

    @property
    def all_default_accessions(self) -> set[str]:
        """Every default transcript accession, without version."""
        return set().union(*self.default_transcripts.values())

    def gene_info(self, symbols: Iterable[str]) -> dict[str, str]:
        """Return Panel, Phenotypes and Inheritance strings for genes.

        The rules reproduce the baseline MViewer tables: each gene
        contributes its panels (sorted, unique within the gene; not
        deduplicated across genes) and its OMIM phenotype entries as written
        in genemap2 (``; `` joined); genes are joined by ``,`` (panels) or
        ``, `` (phenotypes), repeated texts removed; inheritance is the
        union of codes in INHERITANCE_ORDER.

        Args:
            symbols: Gene symbols (HGNC); unknown symbols are ignored.

        Returns:
            Dict with ``panel``, ``phenotypes`` and ``inheritance``; empty
            strings when nothing matches.
        """
        panels: list[str] = []
        texts: list[str] = []
        codes: set[str] = set()
        for symbol in dict.fromkeys(symbols):
            panels += sorted(set(self.panels.get(symbol, [])))
            entries = self.omim.get(symbol, [])
            texts.append("; ".join(dict.fromkeys(e["text"] for e in entries)))
            for entry in entries:
                codes.update(c for c in entry["inheritance"].split(", ") if c)
        if codes & {"XLR", "XLD"}:
            # Baseline: MViewer matches the generic term "X-linked" inside
            # "X-linked recessive/dominant" too, so XL always accompanies.
            codes.add("XL")
        return {
            "panel": ",".join(panels),
            "phenotypes": ", ".join(dict.fromkeys(t for t in texts if t)),
            "inheritance": ", ".join(sorted(codes, key=_inheritance_sort_key)),
        }


def load_panel_list(
    path: Path, columns: Mapping[str, Sequence[str]] = PANEL_COLUMNS
) -> dict[str, list[str]]:
    """Read PanelList.tsv into gene -> panel names.

    Args:
        path: Tab-separated file with a header line.
        columns: Candidate header names for ``gene`` and ``panel``. A gene
            cell may hold several symbols separated by comma, semicolon or
            space.

    Raises:
        LookupSchemaError: If the header lacks a gene or panel column.
    """
    panels: dict[str, list[str]] = {}
    with open(path, encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").lstrip("#").split("\t")
        index = _find_columns(header, columns, path)
        for line in handle:
            cells = line.rstrip("\n").split("\t")
            if len(cells) <= max(index.values()):
                continue
            panel = cells[index["panel"]].strip()
            for symbol in re.split(r"[,;\s]+", cells[index["gene"]]):
                if symbol and panel:
                    panels.setdefault(symbol, []).append(panel)
    return panels


def load_genemap2(path: Path) -> dict[str, list[dict[str, str]]]:
    """Read OMIM genemap2.txt into gene -> phenotype dicts.

    Args:
        path: genemap2.txt; the column header is the last ``#`` line that
            contains ``Approved Gene Symbol``.

    Raises:
        LookupSchemaError: If that header line is missing.
    """
    header: list[str] | None = None
    omim: dict[str, list[dict[str, str]]] = {}
    aliases: dict[str, list[dict[str, str]]] = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("#"):
                if OMIM_SYMBOL_COLUMN in line:
                    header = [c.strip() for c in line.lstrip("# ").split("\t")]
                continue
            if header is None:
                raise LookupSchemaError(f"{path}: no genemap2 header line")
            cells = line.split("\t")
            row = dict(zip(header, cells, strict=False))
            symbol = row.get(OMIM_SYMBOL_COLUMN, "").strip()
            entries = parse_omim_phenotypes(row.get(OMIM_PHENOTYPE_COLUMN, ""))
            if symbol:
                omim.setdefault(symbol, []).extend(entries)
            for alias in row.get(OMIM_ALIASES_COLUMN, "").split(","):
                if alias.strip():
                    aliases.setdefault(alias.strip(), []).extend(entries)
    # Baseline: a symbol that is not an approved symbol in genemap2 (for
    # example GBA, now GBA1) falls back to the rows listing it as an alias.
    for alias, entries in aliases.items():
        omim.setdefault(alias, entries)
    return omim


def load_default_transcripts(path: Path) -> dict[str, set[str]]:
    """Read DefaultTranscript.tsv into gene -> unversioned accessions.

    The layout is not verified: the first tab-separated cell that is not an
    accession is taken as the gene symbol, and every NM_/NR_/XM_/XR_/ENST
    accession on the line is a default transcript for it. Lines without an
    accession (headers) are skipped.

    Args:
        path: The DefaultTranscript.tsv file.
    """
    out: dict[str, set[str]] = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            cells = line.rstrip("\n").split("\t")
            accessions = {
                m.group(1) for m in DEFAULT_TRANSCRIPT_RE.finditer(line)
            }
            if not accessions:
                continue
            gene = next(
                (
                    c.strip()
                    for c in cells
                    if c.strip() and not DEFAULT_TRANSCRIPT_RE.search(c)
                ),
                "",
            )
            out.setdefault(gene, set()).update(accessions)
    return out


class LookupTables:
    """All lookups for one data directory.

    Use as a context manager so the tabix files are closed.
    """

    def __init__(
        self,
        layout: DataLayout,
        clinvar_columns: Mapping[str, Sequence[str]] = CLINVAR_COLUMNS,
        taiwan_columns: Mapping[str, Sequence[str]] = TAIWAN_BIOBANK_COLUMNS,
    ) -> None:
        """Open the tables named by ``layout``.

        Args:
            layout: Data directory layout.
            clinvar_columns: Candidate header names for ClinVar.tsv.gz.
            taiwan_columns: Candidate header names for the Taiwan Biobank
                table.
        """
        self.clinvar = VariantTable(layout.clinvar, clinvar_columns)
        self.taiwan = AlleleFreqTable(layout.taiwan_biobank, taiwan_columns)
        self.hgmd = pysam.VariantFile(str(layout.hgmd))
        self.genes = GeneTables(
            panels=load_panel_list(layout.panel_list),
            omim=load_genemap2(layout.genemap2),
            default_transcripts=load_default_transcripts(
                layout.default_transcript
            ),
        )

    def __enter__(self) -> LookupTables:
        """Return self."""
        return self

    def __exit__(self, *exc: object) -> None:
        """Close the open files."""
        self.close()

    def close(self) -> None:
        """Close the tabix files."""
        self.clinvar.close()
        self.taiwan.close()
        self.hgmd.close()

    def hgmd_class(self, chrom: str, pos: int, ref: str, alt: str) -> str:
        """Return the HGMD variant class (``DM``, ``DM?``, ...) or ``""``.

        Matching is on position and alleles after minimal normalisation.
        """
        target = minimal_allele(pos, ref, alt)
        for contig in _candidate_contigs(chrom):
            try:
                records = list(
                    self.hgmd.fetch(contig, max(0, pos - 2), pos + len(ref) + 2)
                )
            except ValueError:
                continue
            for rec in records:
                for rec_alt in rec.alts or ():
                    if minimal_allele(rec.pos, rec.ref, rec_alt) == target:
                        value = rec.info.get(HGMD_CLASS_KEY, "")
                        if isinstance(value, tuple):
                            value = ",".join(str(v) for v in value)
                        return str(value)
            break
        return ""

    def variant_fields(
        self, chrom: str, pos: int, ref: str, alt: str
    ) -> dict[str, str]:
        """Return ClinVar, ClinVar_VID, Taiwan Biobank and HGMD values.

        The first matching row of each table is used; missing values are
        empty strings.

        Returns:
            Dict with ``clinvar``, ``clinvar_vid``, ``taiwan_biobank`` and
            ``hgmd_class``.
        """
        clinvar = self.clinvar.query(chrom, pos, ref, alt)
        return {
            "clinvar": clinvar[0]["clinvar"] if clinvar else "",
            "clinvar_vid": clinvar[0]["clinvar_vid"] if clinvar else "",
            "taiwan_biobank": self.taiwan.query(chrom, pos, ref, alt),
            "hgmd_class": self.hgmd_class(chrom, pos, ref, alt),
        }

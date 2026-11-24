# SPDX-License-Identifier: GPL-3.0-or-later
"""Draw read-pileup figures per region or locus.

One JPEG per region (or per locus when the haplotypes of a region are
projected to different loci). Top to bottom: gene track, allele track,
and one band per haplotype (coverage plus packed reads) ordered gene
first, then paralog, pseudogene and fusion haplotypes. Haplotype roles
come from the Paraphase HP names (``resources/haplotype_roles.tsv``).
Coordinates: VCF positions are 1-based; read spans from pysam are 0-based
half-open; the x axis is in GRCh38 bp.

Port of the former Nextflow process ``PLOT_REGION``.
"""

from __future__ import annotations

import bisect
import collections
import csv
import functools
import json
import math
import os
import random
import re
import string
import zlib
from collections.abc import Mapping, Sequence
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field, replace
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np  # noqa: E402
import pysam  # noqa: E402
from matplotlib import pyplot as plt  # noqa: E402
from matplotlib.axes import Axes  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402
from matplotlib.ticker import FuncFormatter, MaxNLocator  # noqa: E402

from paralog2viewer.datalayout import allele_key, bare_chrom  # noqa: E402
from paralog2viewer.figstyle import (  # noqa: E402
    FIGURE_DPI,
    FONT_PT,
    MIXED_COLOR,
    ROLE_COLOR,
    base_rc,
    tint,
)

RESOURCES = Path(__file__).parent / "resources"
READ_COLOR = "#7F7F7F"
BAND_COLOR = "#F2F2F2"
LINE_COLOR = "#4D4D4D"
GUIDE_COLOR = "#BDBDBD"
GRAY_TICK = "#6E6E6E"
NEUTRAL_GENE = "#7F7F7F"
MM = 25.4
PT_MM = MM / 72
WIDTH_MM = 180.0
MAX_HEIGHT_MM = 247.0
MIN_SV_BP = 30
# Reads drawn per haplotype, and the packing rule for their rows.
MAX_READS = 40
STUB_FRACTION = 0.1
MIN_GENE_FRACTION = 0.03
PACK_GAP_BP = 300
MAX_PITCH_MM = 1.0
MIN_PITCH_MM = 0.3
# Layout in mm.
LEFT_MM = 50.0
TOP_MM = 2.0
SECTION_GAP_MM = 1.5
BAND_GAP_MM = 2.0
MIN_BAND_MM = 5.0
MAX_BAND_MM = 26.0
GENE_ROW_MM = 10.0
GENE_BOX_Y_MM = 7.6
MAX_GENE_ROWS = 4
SUBTITLE_MM = 3.8
# Reads area of a band: rows x this pitch, clamped.
TARGET_PITCH_MM = 0.8
MIN_READS_MM = 4.0
MAX_READS_MM = 21.0
# Okabe-Ito then grays, one per haplotype when all share one role.
HAP_COLORS = (
    "#009E73",
    "#CC79A7",
    "#56B4E9",
    "#E69F00",
    "#D55E00",
    "#0072B2",
    "#F0E442",
    "#595959",
    "#9A9A9A",
)
VAR_MM = 10.5
AXIS_MM = 9.0
LEGEND_ROW_MM = 3.6
KEY_ROW_MM = 3.2
LINE_MM = 3.0
CHAR_MM = 1.3
EXON_LABEL_GAP_MM = 2.9
LETTER_GAP_MM = 3.0
MAX_LETTERS = 18


@dataclass(frozen=True, eq=False)
class Allele:
    """One called allele (identity-hashed, so it can sit in id sets).

    Attributes:
        chrom: Contig name.
        pos: 1-based VCF POS.
        end: 1-based inclusive end (equals ``pos + len(ref) - 1`` unless
            the record has an END or SVLEN).
        ref: REF allele.
        alt: ALT allele (symbolic or breakend text kept verbatim).
        symbolic: True for ``<..>`` and breakend alleles.
        haps: Haplotypes whose GT carries this allele.
        plp: True if the allele key is in the pathogenic set.
        label: Text in the key for P/LP and known alleles.
    """

    chrom: str
    pos: int
    end: int
    ref: str
    alt: str
    symbolic: bool
    haps: frozenset[str]
    plp: bool = False
    label: str = ""


@dataclass
class HapReads:
    """Reads of one haplotype, grouped by read name.

    Attributes:
        rows: One list of alignments (primary and supplementary) per read,
            sorted by start then end.
        chrom: Contig carrying most reads.
    """

    rows: list[list[pysam.AlignedSegment]] = field(default_factory=list)
    chrom: str = ""


# ---------------------------------------------------------------- inputs


def _read_tsv(path: Path | str) -> list[dict[str, str]]:
    """Read a tab-separated file with a header, skipping ``#`` lines."""
    with open(path, newline="", encoding="utf-8") as handle:
        lines = (ln for ln in handle if not ln.startswith("#"))
        return list(csv.DictReader(lines, delimiter="\t"))


KNOWN_COLUMNS = (
    "functional_gene",
    "hgvs",
    "label",
    "chrom",
    "pos",
    "ref",
    "alt",
    "reference",
)


def read_known(path: Path | str) -> list[dict[str, str]]:
    """Read known_gene_pseudogene_variants.tsv.

    Args:
        path: Tab-separated file; ``#`` lines are comments; needs the
            columns in ``KNOWN_COLUMNS`` (pos is 1-based).

    Returns:
        One dict per row with all columns of the file (strings).

    Raises:
        ValueError: If a required column is missing or a pos is not an
            integer; the message names the file, column or row.
    """
    with open(path, newline="", encoding="utf-8") as handle:
        lines = (ln for ln in handle if not ln.startswith("#"))
        reader = csv.DictReader(lines, delimiter="\t")
        header = reader.fieldnames or []
        rows = list(reader)
    missing = [c for c in KNOWN_COLUMNS if c not in header]
    if missing:
        raise ValueError(f"{path}: missing column(s) {', '.join(missing)}")
    for number, row in enumerate(rows, start=1):
        if not row["pos"].isdigit():
            raise ValueError(f"{path}: row {number}: pos {row['pos']!r}")
    return rows


def known_hits(
    known_rows: Sequence[Mapping[str, str]],
    alleles: Sequence[Allele],
    region_of: Mapping[str, str] | None = None,
    roles: Mapping[tuple[str, str], HapRole] | None = None,
    hap_symbols: Mapping[str, str] | None = None,
) -> list[tuple[str, str, int, str, str, str, str, str]]:
    """List known gene/pseudogene variants present in the merged VCF.

    A row is a hit when an allele has the identical chrom, pos, ref and
    alt and at least one haplotype carrying it belongs to the row's
    ``functional_gene``: role ``gene`` and the same symbol in the
    haplotype role table, or the symbol given in ``hap_symbols``. A
    haplotype with no resolved role (for example NCF1 copies) never
    qualifies, so a pseudogene copy is not mistaken for the gene. The
    allele is the pseudogene-type allele inside the functional gene.

    Args:
        known_rows: Output of ``read_known``.
        alleles: Output of ``read_alleles`` (haplotypes that carry them).
        region_of: Optional HP -> region name for roles of HP names
            without a region part.
        roles: Role table (``load_roles``); default the packaged one.
        hap_symbols: Optional HP -> gene symbol that replaces the role
            table (used for RCCX haplotypes, see ``rccx.cyp21_symbols``).

    Returns:
        (gene, chrom, pos, ref, alt, label, hgvs, reference) per hit, in
        file order.
    """
    return [
        (
            row["functional_gene"],
            row["chrom"],
            int(row["pos"]),
            row["ref"],
            row["alt"],
            row["label"],
            row["hgvs"],
            row["reference"],
        )
        for row, _allele, _haps in _known_carriers(
            known_rows, alleles, region_of, roles, hap_symbols
        )
    ]


def _known_carriers(
    known_rows: Sequence[Mapping[str, str]],
    alleles: Sequence[Allele],
    region_of: Mapping[str, str] | None,
    roles: Mapping[tuple[str, str], HapRole] | None,
    hap_symbols: Mapping[str, str] | None = None,
) -> list[tuple[Mapping[str, str], Allele, frozenset[str]]]:
    """Return (row, allele, carrier haplotypes of the functional gene)."""
    table = roles if roles is not None else _default_roles()
    region_of = region_of or {}
    hap_symbols = hap_symbols or {}

    def symbol_of(hap: str) -> str:
        if hap in hap_symbols:
            return hap_symbols[hap].upper()
        role = resolve_role(hap, table, region_of.get(hap, ""))
        return role.symbol.upper() if role.role == "gene" else ""

    by_site: dict[tuple[str, int, str, str], list[Allele]] = (
        collections.defaultdict(list)
    )
    for a in alleles:
        by_site[(a.chrom, a.pos, a.ref, a.alt)].append(a)
    found = []
    for row in known_rows:
        gene = row["functional_gene"].upper()
        site = (row["chrom"], int(row["pos"]), row["ref"], row["alt"])
        for a in by_site.get(site, ()):
            carriers = frozenset(h for h in a.haps if symbol_of(h) == gene)
            if carriers:
                found.append((row, a, carriers))
    return found


def read_alleles(
    merged_vcf: Path | str,
    pathogenic_keys: set[str],
    labels: Mapping[str, str] | None = None,
) -> list[Allele]:
    """Read alleles carried (GT 1) by at least one haplotype.

    Args:
        merged_vcf: Merged VCF; sample columns are HP values, haploid GT.
        pathogenic_keys: Allele keys classified P/LP.
        labels: Optional allele key -> label for P/LP markers.

    Returns:
        Alleles sorted by (chrom, pos). Multi-allelic records are split by
        the GT allele index.
    """
    labels = labels or {}
    alleles = []
    with pysam.VariantFile(str(merged_vcf)) as vcf:
        samples = list(vcf.header.samples)
        for rec in vcf:
            carriers: dict[int, list[str]] = collections.defaultdict(list)
            for name in samples:
                for index in rec.samples[name]["GT"] or ():
                    if index:
                        carriers[index].append(name)
            for index, haps in carriers.items():
                alt = rec.alts[index - 1]
                symbolic = (
                    alt.startswith("<")
                    or "[" in alt
                    or "]" in alt
                    or alt == "*"
                )
                key = allele_key(rec.chrom, rec.pos, rec.ref, alt)
                alleles.append(
                    Allele(
                        rec.chrom,
                        rec.pos,
                        rec.stop,
                        rec.ref,
                        alt,
                        symbolic,
                        frozenset(haps),
                        key in pathogenic_keys,
                        labels.get(key, ""),
                    )
                )
    return sorted(alleles, key=lambda a: (a.chrom, a.pos))


def load_reads(bam_path: Path | str) -> HapReads:
    """Read mapped, non-secondary alignments of one haplotype BAM once.

    Args:
        bam_path: reads.bam of the haplotype.

    Returns:
        Reads on the contig that carries most of them, row-sorted.
    """
    per_name: dict[str, list[pysam.AlignedSegment]] = {}
    with pysam.AlignmentFile(str(bam_path)) as bam:
        for read in bam.fetch(until_eof=True):
            if not (read.is_unmapped or read.is_secondary):
                per_name.setdefault(read.query_name, []).append(read)
    counts = collections.Counter(
        a.reference_name for group in per_name.values() for a in group
    )
    if not counts:
        return HapReads()
    chrom = counts.most_common(1)[0][0]
    rows = []
    for group in per_name.values():
        group = [a for a in group if a.reference_name == chrom]
        if group:
            rows.append(group)
    rows.sort(key=_bounds)
    return HapReads(rows, chrom)


def _bounds(group: Sequence[pysam.AlignedSegment]) -> tuple[int, int]:
    """Return the 0-based half-open span of all alignments of one read."""
    return (
        min(a.reference_start for a in group),
        max(a.reference_end for a in group),
    )


def span(reads: HapReads) -> tuple[int, int]:
    """Return the 0-based half-open span of the reads, minus stray stubs.

    Reads shorter than ``STUB_FRACTION`` of the median read span do not
    widen the span (they stay drawn where they fall inside it).
    """
    bounds = [_bounds(g) for g in reads.rows]
    lengths = sorted(e - s for s, e in bounds)
    floor = STUB_FRACTION * lengths[len(lengths) // 2]
    kept = [b for b in bounds if b[1] - b[0] >= floor] or bounds
    return min(s for s, _ in kept), max(e for _, e in kept)


def cluster_loci(spans: Mapping[str, tuple[str, int, int]]) -> list[dict]:
    """Group haplotypes whose callable spans overlap into loci.

    Haplotypes projected to different loci (native gene2 coordinates)
    cannot share an x axis.

    Args:
        spans: HP -> (chrom, start0, end0).

    Returns:
        Loci sorted by position, each a dict with chrom, start0, end0, haps.
    """
    loci: list[dict] = []
    for hap, (chrom, start0, end0) in sorted(
        spans.items(), key=lambda kv: (kv[1][0], kv[1][1], kv[0])
    ):
        last = loci[-1] if loci else None
        if last and last["chrom"] == chrom and start0 < last["end0"]:
            last["haps"].append(hap)
            last["end0"] = max(last["end0"], end0)
        else:
            loci.append(
                {"chrom": chrom, "start0": start0, "end0": end0, "haps": [hap]}
            )
    return loci


# ------------------------------------------------------- read-vs-allele


def alt_state(
    read: pysam.AlignedSegment,
    ref_to_query: Mapping[int, int],
    s0: int,
    ref: str,
    alt: str,
) -> bool | None:
    """Tell whether a read carries a sequence allele.

    Args:
        read: Alignment with a stored query sequence.
        ref_to_query: 0-based reference position -> query index (M only).
        s0: 0-based anchor base (VCF POS - 1).
        ref: REF allele.
        alt: ALT allele.

    Returns:
        True if the read carries ALT; False if it covers the site with REF
        or another allele; None if it does not cover the site.
    """
    if not read.reference_start <= s0 < read.reference_end:
        return None
    if read.query_sequence is None:
        return False
    seq = read.query_sequence.upper()
    anchor = ref_to_query.get(s0)
    if len(ref) == len(alt):
        if any(s0 + i not in ref_to_query for i in range(len(ref))):
            return False
        return all(
            seq[ref_to_query[s0 + i]] == alt[i].upper() for i in range(len(ref))
        )
    right = ref_to_query.get(s0 + len(ref))
    if anchor is None or right is None:
        return False
    # Bases the read has between the anchor and the first base after REF,
    # whatever the CIGAR looks like, must be the ALT bases after the anchor.
    observed = seq[anchor + 1 : right]
    expected = alt[1:].upper()
    if len(observed) != len(expected):
        return False
    return sum(
        a == b for a, b in zip(observed, expected, strict=True)
    ) >= 0.8 * len(expected)


def sv_state(
    alignments: Sequence[pysam.AlignedSegment], s0: int, e0: int
) -> bool | None:
    """Tell whether a read supports a symbolic or breakend allele.

    A long CIGAR indel overlapping the call, or a split-alignment boundary
    next to it, counts as support.

    Args:
        alignments: All alignments of one read.
        s0: 0-based start of the call.
        e0: 0-based end of the call.

    Returns:
        True (supports), False (covers, no support), None (no coverage).
    """
    covered = any(a.reference_start <= s0 < a.reference_end for a in alignments)
    if not covered and len(alignments) < 2:
        return None
    flank = 100
    for a in alignments:
        position = a.reference_start
        for op, length in a.cigartuples or ():
            if op in (0, 7, 8):
                position += length
            elif op == 2:
                near = position < e0 + flank and position + length > s0 - flank
                if length >= MIN_SV_BP and near:
                    return True
                position += length
            elif op == 1 and length >= MIN_SV_BP:
                if s0 - flank <= position <= e0 + flank:
                    return True
            elif op == 3:
                position += length
    if len(alignments) > 1:
        for a in alignments:
            gap = min(
                abs(a.reference_start - s0),
                abs(a.reference_end - s0),
                abs(a.reference_start - e0),
                abs(a.reference_end - e0),
            )
            if gap <= 200:
                return True
    return False


# ------------------------------------------------- haplotype roles (names)

ROLE_ORDER = (
    "gene",
    "paralog_gene",
    "pseudogene",
    "fusion",
    "other",
    "unresolved",
)
ROLE_TEXT = {
    "gene": "gene",
    "paralog_gene": "paralog",
    "pseudogene": "pseudogene",
    "fusion": "fusion",
    "other": "other",
    "unresolved": "unresolved",
}
ALLELE_TEXT = {
    "gene": "gene-only allele",
    "paralog_gene": "paralog-only allele",
    "pseudogene": "pseudogene-only allele",
    "fusion": "fusion-only allele",
}
_HAP_RE = re.compile(
    r"^(?:(?P<region>.+)_)?(?P<prefix>[A-Za-z0-9]*hap)(?P<number>\d+)$"
)


@dataclass(frozen=True)
class HapRole:
    """Display symbol and role of one haplotype.

    Attributes:
        symbol: Gene symbol shown in the label (the region name if unknown).
        role: One of ``ROLE_ORDER``.
    """

    symbol: str
    role: str


def parse_hap_name(hap: str) -> tuple[str, str, int]:
    """Split a Paraphase HP name such as ``pms2_pms2clhap2``.

    Returns:
        (region, prefix, number); region is ``""`` without an underscore
        and prefix is ``""`` (number 0) if the name has no ``hap`` part.
    """
    match = _HAP_RE.match(hap)
    if not match:
        return "", "", 0
    return (
        match["region"] or "",
        match["prefix"],
        int(match["number"]),
    )


def load_roles(
    path: Path | str | None = None,
) -> dict[tuple[str, str], HapRole]:
    """Read haplotype_roles.tsv.

    Args:
        path: TSV with columns region, prefix, symbol, role, source;
            ``#`` lines are comments. Default: the packaged table.

    Returns:
        (region lower case or ``*``, prefix lower case) -> HapRole.
    """
    path = Path(path) if path else RESOURCES / "haplotype_roles.tsv"
    table = {}
    for row in _read_tsv(path):
        if row["role"] not in ROLE_ORDER:
            raise ValueError(f"unknown role {row['role']!r} in {path}")
        key = (row["region"].lower(), row["prefix"].lower())
        table[key] = HapRole(row["symbol"], row["role"])
    return table


def resolve_role(
    hap: str,
    table: Mapping[tuple[str, str], HapRole],
    region_hint: str = "",
) -> HapRole:
    """Return the role of a haplotype from its HP name.

    Args:
        hap: HP value, ``<region>_<prefix><number>``.
        table: Output of ``load_roles``.
        region_hint: Region used when the HP name carries none.

    Returns:
        The table row, with the region name as symbol when the table
        leaves it empty; ``unresolved`` when no row matches.
    """
    region, prefix, _ = parse_hap_name(hap)
    region = region or region_hint or hap
    for key in ((region.lower(), prefix.lower()), ("*", prefix.lower())):
        row = table.get(key)
        if row:
            return HapRole(row.symbol or region, row.role)
    return HapRole(region, "unresolved")


def order_haps(haps: Sequence[str], roles: Mapping[str, HapRole]) -> list[str]:
    """Sort haplotypes: gene first, then paralog, pseudogene, fusion, rest."""

    def key(hap: str) -> tuple:
        _, prefix, number = parse_hap_name(hap)
        return (ROLE_ORDER.index(roles[hap].role), prefix, number, hap)

    return sorted(haps, key=key)


def load_hap_info(json_path: Path | str) -> dict[str, dict]:
    """Collect per-haplotype facts from a Paraphase JSON.

    Args:
        json_path: ``<sample>.paraphase.json``.

    Returns:
        HP -> dict with region, boundary, boundary_gene2, is_truncated
        (list of ``'5p'``/``'3p'``), haplogroup, and two_copy (True if
        Paraphase lists the haplotype in ``two_copy_haplotypes``).
    """
    with open(json_path, encoding="utf-8") as handle:
        data = json.load(handle)
    info: dict[str, dict] = {}
    for region, rec in data.items():
        details = rec.get("haplotype_details") if isinstance(rec, dict) else 0
        if not isinstance(details, dict):
            continue
        final = rec.get("final_haplotypes")
        final = final if isinstance(final, dict) else {}
        two = {final.get(k, k) for k in rec.get("two_copy_haplotypes") or ()}
        for hap, det in details.items():
            info[hap] = {
                "region": region,
                "boundary": det.get("boundary"),
                "boundary_gene2": det.get("boundary_gene2"),
                "is_truncated": list(det.get("is_truncated") or ()),
                "haplogroup": det.get("haplogroup"),
                "two_copy": hap in two,
            }
    return info


# ----------------------------------------------------------- gene models


@dataclass(frozen=True)
class Transcript:
    """One transcript of a refGene table (0-based half-open coordinates)."""

    symbol: str
    accession: str
    chrom: str
    strand: str
    start0: int
    end0: int
    cds_start0: int
    cds_end0: int
    exons: tuple[tuple[int, int], ...]

    @property
    def exon_bp(self) -> int:
        """Total exon length."""
        return sum(e - s for s, e in self.exons)


def read_default_transcripts(path: Path | str) -> dict[str, str]:
    """Read DefaultTranscript.tsv (Name, Chr, DefaultTranscript).

    Returns:
        Gene symbol -> accession without the version suffix.
    """
    out = {}
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            out[row["Name"]] = row["DefaultTranscript"].split(".")[0]
    return out


def read_gene_models(
    refgene: Path | str,
    default_transcripts: Mapping[str, str] | None,
    chrom: str,
    start0: int,
    end0: int,
    prefer: Sequence[str] = (),
    max_genes: int = MAX_GENE_ROWS,
) -> list[Transcript]:
    """Pick one transcript per gene that overlaps a window.

    Args:
        refgene: UCSC refGene table (``hg38_refGene.txt``).
        default_transcripts: Gene symbol -> preferred accession; a gene
            without one uses its transcript with the most exon bases.
        chrom: Contig of the window.
        start0: 0-based window start.
        end0: 0-based exclusive window end.
        prefer: Gene symbols kept first when more than ``max_genes``
            genes overlap (those of the haplotypes on display).
        max_genes: Upper bound on genes returned.

    Returns:
        Transcripts sorted by start.
    """
    want = bare_chrom(chrom)
    by_gene: dict[str, list[Transcript]] = collections.defaultdict(list)
    with open(refgene, encoding="utf-8") as handle:
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if len(f) < 13 or bare_chrom(f[2]) != want:
                continue
            if int(f[5]) <= start0 or int(f[4]) >= end0:
                continue
            starts = [int(x) for x in f[9].split(",") if x]
            ends = [int(x) for x in f[10].split(",") if x]
            by_gene[f[12]].append(
                Transcript(
                    f[12],
                    f[1],
                    f[2],
                    f[3],
                    int(f[4]),
                    int(f[5]),
                    int(f[6]),
                    int(f[7]),
                    tuple(zip(starts, ends, strict=True)),
                )
            )
    chosen = []
    for symbol, items in by_gene.items():
        wanted = (default_transcripts or {}).get(symbol)
        best = [t for t in items if t.accession.split(".")[0] == wanted]
        chosen.append(best[0] if best else max(items, key=lambda t: t.exon_bp))
    favored = {s.upper() for s in prefer}
    window = max(end0 - start0, 1)
    chosen = [
        t
        for t in chosen
        if t.symbol.upper() in favored
        or min(t.end0, end0) - max(t.start0, start0)
        >= MIN_GENE_FRACTION * window
    ]

    def rank(t: Transcript) -> tuple:
        inside = min(t.end0, end0) - max(t.start0, start0)
        return (t.symbol.upper() not in favored, -inside, t.symbol)

    return sorted(sorted(chosen, key=rank)[:max_genes], key=lambda t: t.start0)


# ------------------------------------------------------- read selection


def sample_reads(
    rows: Sequence[list[pysam.AlignedSegment]], cap: int, key: str
) -> list[list[pysam.AlignedSegment]]:
    """Choose at most ``cap`` reads, the same ones on every run.

    Args:
        rows: Alignment groups, one per read.
        cap: Maximum number of reads kept.
        key: Seed text (the haplotype name).

    Returns:
        All rows if ``len(rows) <= cap``; otherwise a random sample drawn
        from the rows sorted by read name with a seed derived from ``key``.
    """
    if len(rows) <= cap:
        return list(rows)
    ordered = sorted(rows, key=lambda g: g[0].query_name)
    rng = random.Random(zlib.crc32(key.encode()))
    return [ordered[i] for i in sorted(rng.sample(range(len(ordered)), cap))]


def pack_rows(
    spans: Sequence[tuple[int, int]], gap: int = PACK_GAP_BP
) -> list[int]:
    """Assign reads to rows so that reads in one row do not overlap.

    Args:
        spans: (start, end) per read, same coordinate system.
        gap: Minimum free bases between two reads of one row.

    Returns:
        Row index per read; row 0 holds the leftmost read.
    """
    ends: list[int] = []
    rows = [0] * len(spans)
    for i in sorted(range(len(spans)), key=lambda i: spans[i]):
        start, end = spans[i]
        for r, last in enumerate(ends):
            if last + gap <= start:
                ends[r], rows[i] = end, r
                break
        else:
            ends.append(end)
            rows[i] = len(ends) - 1
    return rows


def coverage(
    reads: HapReads, x_lo: int, x_hi: int, bins: int = 400
) -> tuple[np.ndarray, np.ndarray]:
    """Return alignment depth at ``bins`` evenly spaced positions.

    Counts every alignment of every read, not only the drawn sample.
    """
    grid = np.linspace(x_lo, x_hi, bins)
    starts = np.sort([a.reference_start for g in reads.rows for a in g])
    ends = np.sort([a.reference_end for g in reads.rows for a in g])
    depth = np.searchsorted(starts, grid, "right") - np.searchsorted(
        ends, grid, "right"
    )
    return grid, depth


def read_hits(
    alignments: Sequence[pysam.AlignedSegment],
    carriers: Sequence[Allele],
    positions: Sequence[int],
) -> list[Allele]:
    """List the carried alleles one read supports.

    Args:
        alignments: All alignments of the read.
        carriers: Alleles called in the read's haplotype, sorted by pos.
        positions: ``[v.pos for v in carriers]``.

    Returns:
        The alleles whose ALT the read carries (SV: supports).
    """
    lo, hi = _bounds(alignments)
    maps: dict[int, dict[int, int]] = {}
    hits = []
    first = bisect.bisect_right(positions, lo)
    for v in carriers[first : bisect.bisect_right(positions, hi)]:
        s0 = v.pos - 1
        state = None
        if v.symbolic:
            state = sv_state(alignments, s0, v.end - 1)
        else:
            for a in alignments:
                if a.reference_start <= s0 < a.reference_end:
                    if id(a) not in maps:
                        maps[id(a)] = {
                            r: q
                            for q, r in a.get_aligned_pairs(matches_only=True)
                        }
                    state = alt_state(a, maps[id(a)], s0, v.ref, v.alt)
                    break
        if state:
            hits.append(v)
    return hits


def allele_color(roles: set[str]) -> str:
    """Tick color of an allele from the roles of the haplotypes carrying it."""
    if len(roles) > 1:
        return MIXED_COLOR
    (role,) = roles
    return "#000000" if role in ("other", "unresolved") else ROLE_COLOR[role]


def allele_text(roles: set[str]) -> str:
    """Legend text for ``allele_color(roles)``."""
    if len(roles) > 1:
        return "allele in several roles"
    (role,) = roles
    return ALLELE_TEXT.get(role, "allele in some haplotypes")


# --------------------------------------------------------------- drawing


def _set_style() -> None:
    plt.rcParams.update(
        {
            **base_rc(),
            "axes.labelsize": FONT_PT,
            "xtick.labelsize": FONT_PT,
            "ytick.labelsize": FONT_PT,
            "legend.fontsize": FONT_PT,
            "xtick.major.width": 0.5,
            "xtick.major.size": 2.0,
            "ytick.major.width": 0.5,
            "ytick.major.size": 1.5,
            "figure.dpi": FIGURE_DPI,
        }
    )


def _marker(
    shape: str, size: float, face: str, edge: str, width: float, label: str
) -> Line2D:
    return Line2D(
        [],
        [],
        marker=shape,
        ls="",
        markersize=size,
        markerfacecolor=face,
        markeredgecolor=edge,
        markeredgewidth=width,
        label=label,
    )


def _truncation_text(flags: Sequence[str]) -> str:
    ends = [{"5p": "5'", "3p": "3'"}.get(f, f) for f in flags]
    return f"{' and '.join(ends)} end truncated" if ends else ""


def _draw_transcript(
    ax: Axes,
    tx: Transcript,
    y0: float,
    color: str,
    bp_per_mm: float,
    x_lo: int,
    x_hi: int,
) -> None:
    """Draw exons (boxes, thin UTR), intron line, chevrons, exon numbers.

    ``y0`` is the row top in mm; the axis y runs downward.
    """
    mid = y0 + GENE_BOX_Y_MM
    ax.plot([tx.start0, tx.end0], [mid, mid], lw=0.5, color=color, zorder=2)
    coding = tx.cds_end0 > tx.cds_start0
    for s, e in tx.exons:
        parts = []
        if coding:
            c0, c1 = max(s, tx.cds_start0), min(e, tx.cds_end0)
            if c1 > c0:
                parts.append((c0, c1, 2.8))
            parts += [(s, min(e, c0), 1.4), (max(s, c1), e, 1.4)]
            parts = [p for p in parts if p[1] > p[0]]
        else:
            parts.append((s, e, 2.0))
        for a, b, height in parts:
            ax.add_patch(
                Rectangle(
                    (a, mid - height / 2),
                    b - a,
                    height,
                    facecolor=color,
                    edgecolor="none",
                    zorder=3,
                )
            )
    glyph = ">" if tx.strand == "+" else "<"
    for (_, a), (b, _) in zip(tx.exons[:-1], tx.exons[1:], strict=True):
        centre = (a + b) / 2
        if (b - a) / bp_per_mm >= 3.0 and x_lo < centre < x_hi:
            ax.plot(
                [centre],
                [mid],
                ls="",
                marker=glyph,
                markersize=2.2,
                color=color,
                zorder=3,
            )
    last = [-99.0, -99.0]
    count = len(tx.exons)
    width_mm = (x_hi - x_lo) / bp_per_mm
    for i, (s, e) in enumerate(tx.exons):
        centre = (s + e) / 2
        at = (centre - x_lo) / bp_per_mm
        if not 1.5 < at < width_mm - 1.5:
            continue
        # First row, else a second row above it; skip if both are taken.
        row = next(
            (r for r in (0, 1) if at - last[r] >= EXON_LABEL_GAP_MM), None
        )
        if row is None:
            continue
        last[row] = at
        number = i + 1 if tx.strand == "+" else count - i
        ax.text(
            centre,
            mid - 1.9 - 2.6 * row,
            str(number),
            ha="center",
            va="bottom",
            fontsize=FONT_PT,
            color="black",
            clip_on=True,
        )


def _legend_handles(
    text_of: Mapping[str, str],
    discr: Sequence[Allele],
    plp: Sequence[Allele],
    known_here: Sequence[Allele],
    n_plain: int,
) -> list:
    """Return the legend entries of a figure, in drawing order."""
    handles = [
        Line2D(
            [],
            [],
            ls="",
            marker="|",
            markersize=7,
            markeredgewidth=1.2,
            color=color,
            label=text,
        )
        for color, text in text_of.items()
    ]
    if any(v.symbolic for v in discr):
        handles.append(
            _marker("v", 3.5, "white", "black", 0.4, "structural allele")
        )
    if plp:
        handles.append(
            _marker(
                "D", 3.2, "white", "black", 0.6, "P/LP (ClinVar or InterVar)"
            )
        )
    if known_here:
        handles.append(
            _marker("s", 3.0, "black", "white", 0.3, "pseudogene-type site")
        )
    handles.append(
        Patch(
            facecolor=tint(ROLE_COLOR["unresolved"], 0.6),
            label="read depth (in the band color)",
        )
    )
    if n_plain:
        handles.append(
            Line2D(
                [], [], ls="", label=f"{n_plain} more named sites not lettered"
            )
        )
    return handles


def _at_x(pos: int) -> float:
    """Plot x of a 1-based VCF position (centre of the base)."""
    return pos - 0.5


def _draw_allele_track(
    ax_v: Axes,
    x_lo: int,
    x_hi: int,
    relevant: Sequence[Allele],
    known_ids: set[int],
    color_of: Mapping[int, str],
    plp: Sequence[Allele],
    known_here: Sequence[Allele],
    lettered: Sequence[Allele],
    quiet: bool,
) -> None:
    """Draw ticks, SV bars, P/LP and known markers, and site letters.

    ``quiet`` adds the note that no allele separates the haplotypes.
    """
    ax_v.set_xlim(x_lo, x_hi)
    ax_v.set_ylim(0, VAR_MM)
    ax_v.axis("off")
    ax_v.axhline(0.15, lw=0.4, color=LINE_COLOR)
    ticks = collections.defaultdict(list)
    bars = collections.defaultdict(list)
    for v in relevant:
        if id(v) in known_ids:
            continue
        color = color_of.get(id(v), GRAY_TICK)
        ticks[color].append([(_at_x(v.pos), 0.4), (_at_x(v.pos), 4.4)])
        if v.symbolic and v.end > v.pos:
            bars[color].append([(_at_x(v.pos), 2.4), (_at_x(v.end), 2.4)])
    for color, segs in ticks.items():
        ax_v.add_collection(LineCollection(segs, colors=color, linewidths=0.4))
    for color, segs in bars.items():
        ax_v.add_collection(LineCollection(segs, colors=color, linewidths=1.6))
    plp_x = [_at_x(v.pos) for v in plp]
    known_x = [_at_x(v.pos) for v in known_here]
    marker_y = 5.7
    if plp_x:
        ax_v.scatter(
            plp_x,
            [marker_y] * len(plp_x),
            s=3.2**2,
            marker="D",
            c="white",
            edgecolors="black",
            linewidths=0.5,
            zorder=5,
        )
    if known_x:
        ax_v.scatter(
            known_x,
            [marker_y] * len(known_x),
            s=3.0**2,
            marker="s",
            c="black",
            edgecolors="white",
            linewidths=0.3,
            zorder=5,
        )
    for letter, v in zip(string.ascii_uppercase, lettered, strict=False):
        ax_v.text(
            _at_x(v.pos),
            marker_y + 1.5,
            letter,
            ha="center",
            va="bottom",
            fontsize=FONT_PT,
        )
    if quiet:
        ax_v.text(
            0.5,
            0.2,
            "no allele differs between the haplotypes shown",
            transform=ax_v.transAxes,
            ha="center",
            va="center",
            fontsize=FONT_PT,
            color=LINE_COLOR,
        )


def _draw_reads(
    ax_r: Axes,
    chosen: Sequence[Sequence[pysam.AlignedSegment]],
    row_of: Sequence[int],
    carriers: Sequence[Allele],
    positions: Sequence[int],
    sep_ids: set[int],
    known_ids: set[int],
    color_of: Mapping[int, str],
    pitch_pt: float,
    guides: Sequence[float],
) -> None:
    """Draw read lines, allele ticks, SV, known and P/LP marks, guides."""
    segs, tick_at = [], collections.defaultdict(list)
    plp_xy: tuple[list, list] = ([], [])
    known_xy: tuple[list, list] = ([], [])
    sv_xy = collections.defaultdict(lambda: ([], []))
    for g, row in zip(chosen, row_of, strict=True):
        segs += [[(a.reference_start, row), (a.reference_end, row)] for a in g]
        for v in read_hits(g, carriers, positions):
            x = _at_x(v.pos)
            if id(v) in known_ids:
                known_xy[0].append(x)
                known_xy[1].append(row)
                continue
            if id(v) in sep_ids:
                if v.symbolic:
                    sv_xy[color_of[id(v)]][0].append(x)
                    sv_xy[color_of[id(v)]][1].append(row)
                else:
                    tick_at[color_of[id(v)]].append(
                        [(x, row - 0.45), (x, row + 0.45)]
                    )
            if v.plp:
                plp_xy[0].append(x)
                plp_xy[1].append(row)
    ax_r.add_collection(
        LineCollection(
            segs,
            colors=READ_COLOR,
            linewidths=min(1.8, max(0.45, 0.5 * pitch_pt)),
            zorder=2,
            capstyle="butt",
        )
    )
    for tick_color, tick_segs in tick_at.items():
        ax_r.add_collection(
            LineCollection(
                tick_segs, colors=tick_color, linewidths=0.6, zorder=3
            )
        )
    mark_pt = min(3.6, max(2.4, 1.7 * pitch_pt))
    for tick_color, (xs, ys) in sv_xy.items():
        ax_r.scatter(
            xs,
            ys,
            s=mark_pt**2,
            marker="v",
            c=tick_color,
            edgecolors="black",
            linewidths=0.25,
            zorder=4,
        )
    if known_xy[0]:
        ax_r.scatter(
            *known_xy,
            s=mark_pt**2,
            marker="s",
            c="black",
            edgecolors="white",
            linewidths=0.3,
            zorder=6,
        )
    if plp_xy[0]:
        ax_r.scatter(
            *plp_xy,
            s=(1.1 * mark_pt) ** 2,
            marker="D",
            c="white",
            edgecolors="black",
            linewidths=0.45,
            zorder=5,
        )
    for g in guides:
        ax_r.axvline(g, lw=0.3, ls=(0, (1.5, 1.5)), color=GUIDE_COLOR, zorder=1)


def _short_name(hap: str) -> str:
    """HP name without its region prefix (``pms2_pms2hap1`` -> ``pms2hap1``)."""
    return hap.rsplit("_", 1)[-1]


def draw_locus(
    locus: dict,
    reads: Mapping[str, HapReads],
    alleles: Sequence[Allele],
    known: Sequence[Allele],
    target: Path,
    *,
    hap_roles: Mapping[str, HapRole] | None = None,
    hap_info: Mapping[str, dict] | None = None,
    transcripts: Sequence[Transcript] = (),
    max_reads: int = MAX_READS,
    subtitle: str = "",
) -> Path:
    """Render one locus to a JPEG (150 dpi, 180 mm wide).

    Layout, top to bottom, one shared x axis (GRCh38 bp): gene track,
    allele track, one band per haplotype (coverage plus packed reads),
    axis, legend row, key of lettered sites. Height: fixed parts plus
    ``n_haplotypes x (band + 2 mm)`` with the band clamped to 5-26 mm, so
    the figure is at most 247 mm high for up to about 30 haplotypes.

    Args:
        locus: Dict from ``cluster_loci``.
        reads: HP -> reads.
        alleles: Called alleles (any locus; filtered here).
        known: Known-variant hits as alleles (haps = carrier haplotypes,
            end = pos + len(ref) - 1); see ``plot_regions``.
        target: Output JPEG path.
        hap_roles: HP -> HapRole; default: parsed from the HP names.
        hap_info: HP -> dict from ``load_hap_info`` (truncation, copies).
        transcripts: Gene models to draw (``read_gene_models``).
        max_reads: Reads drawn per haplotype (seeded sample).
        subtitle: Line above the gene track (the paired locus figure).

    Returns:
        ``target``.
    """
    _set_style()
    chrom, x_lo, x_hi = locus["chrom"], locus["start0"], locus["end0"]
    hap_info = hap_info or {}
    table = _default_roles()
    roles = {
        h: (hap_roles or {}).get(h) or resolve_role(h, table)
        for h in locus["haps"]
    }
    names = order_haps(list(locus["haps"]), roles)
    nameset = set(names)
    here = [
        v
        for v in alleles
        if v.chrom == chrom and x_lo < v.pos <= x_hi and v.haps & nameset
    ]
    # Only alleles called in some, not all, haplotypes separate them.
    discr = [v for v in here if not nameset <= v.haps] if len(names) > 1 else []
    sep_ids = {id(v) for v in discr}
    plp = [v for v in here if v.plp]
    known_here = [
        replace(k, haps=k.haps & nameset)
        for k in known
        if k.chrom == chrom
        and k.pos <= x_hi
        and k.end > x_lo
        and k.haps & nameset
    ]
    known_ids = {id(k) for k in known_here}
    relevant = sorted(
        discr + [v for v in plp if id(v) not in sep_ids] + known_here,
        key=lambda v: v.pos,
    )
    # Same role on every haplotype shown: role colors would not compare
    # anything, so color by carrier haplotype instead.
    per_hap = len(names) > 1 and len({roles[h].role for h in names}) == 1
    hap_color = {
        h: HAP_COLORS[i % len(HAP_COLORS)]
        if per_hap
        else ROLE_COLOR[roles[h].role]
        for i, h in enumerate(names)
    }
    color_of: dict[int, str] = {}
    text_of: dict[str, str] = {}
    for v in discr:
        carriers = v.haps & nameset
        if per_hap:
            if len(carriers) == 1:
                (h,) = carriers
                color_of[id(v)] = hap_color[h]
                text_of[hap_color[h]] = f"only in {_short_name(h)}"
            else:
                color_of[id(v)] = "#000000"
                text_of["#000000"] = "allele in several haplotypes"
        else:
            kinds = {roles[h].role for h in carriers}
            color_of[id(v)] = allele_color(kinds)
            text_of[color_of[id(v)]] = allele_text(kinds)

    # Geometry (mm). All x positions share one mapping.
    left_mm, right_mm = LEFT_MM, 3.0
    ax_w = WIDTH_MM - left_mm - right_mm
    bp_per_mm = (x_hi - x_lo) / ax_w

    # Named sites get a letter when it fits; the key lists them.
    named = sorted(
        {id(v): v for v in plp + known_here}.values(), key=lambda v: v.pos
    )
    lettered: list[Allele] = []
    last_mm = -99.0
    for v in named:
        at = (v.pos - 0.5 - x_lo) / bp_per_mm
        if len(lettered) < MAX_LETTERS and at - last_mm >= LETTER_GAP_MM:
            lettered.append(v)
            last_mm = at
    n_plain = len(named) - len(lettered)
    key_mm = KEY_ROW_MM * math.ceil(len(lettered) / 3) + 1 if lettered else 0.0

    handles = _legend_handles(text_of, discr, plp, known_here, n_plain)
    chars = sum(len(h.get_label()) + 5 for h in handles)
    legend_rows = 1 if chars * CHAR_MM <= WIDTH_MM - 6 else 2
    legend_mm = LEGEND_ROW_MM * legend_rows + 1

    # A lone haplotype with no named allele has nothing to put in the track.
    show_vars = len(names) > 1 or bool(plp or known_here)
    var_mm = VAR_MM + SECTION_GAP_MM if show_vars else 0.0
    n_gene = len(transcripts)
    gene_mm = n_gene * GENE_ROW_MM + (SECTION_GAP_MM if n_gene else 0)
    fixed = (
        TOP_MM
        + (SUBTITLE_MM if subtitle else 0.0)
        + gene_mm
        + var_mm
        + AXIS_MM
        + legend_mm
        + key_mm
    )
    n_hap = len(names)
    cov_mm = 4.2 if n_hap <= 8 else 3.0
    plan = {}
    for h in names:
        chosen = sample_reads(reads[h].rows, max_reads, h)
        row_of = pack_rows([_bounds(g) for g in chosen])
        plan[h] = (chosen, row_of, max(row_of, default=0) + 1)
    # A band is as tall as its packed rows need, then all shrink together
    # if the figure would pass the height limit.
    want = {
        h: min(MAX_READS_MM, max(MIN_READS_MM, plan[h][2] * TARGET_PITCH_MM))
        for h in names
    }
    overhead = n_hap * (cov_mm + 0.8 + BAND_GAP_MM)
    room = MAX_HEIGHT_MM - fixed - overhead
    scale = min(1.0, max(room, 2.0 * n_hap) / sum(want.values()))
    reads_of = {h: max(1.5, want[h] * scale) for h in names}
    height_mm = fixed + overhead + sum(reads_of.values())

    fig = plt.figure(figsize=(WIDTH_MM / MM, height_mm / MM))

    def axes_at(top: float, height: float) -> Axes:
        return fig.add_axes(
            [
                left_mm / WIDTH_MM,
                1 - (top + height) / height_mm,
                ax_w / WIDTH_MM,
                height / height_mm,
            ]
        )

    def text_at(x: float, y: float, s: str, **kw: object) -> None:
        fig.text(x / WIDTH_MM, 1 - y / height_mm, s, fontsize=FONT_PT, **kw)

    y = TOP_MM
    if subtitle:
        text_at(2.0, y, subtitle, ha="left", va="top")
        y += SUBTITLE_MM
    if n_gene:
        ax_g = axes_at(y, n_gene * GENE_ROW_MM)
        ax_g.set_xlim(x_lo, x_hi)
        ax_g.set_ylim(n_gene * GENE_ROW_MM, 0)
        ax_g.axis("off")
        by_symbol = {r.symbol.upper(): r.role for r in roles.values()}
        for row, tx in enumerate(transcripts):
            role = by_symbol.get(tx.symbol.upper())
            _draw_transcript(
                ax_g,
                tx,
                row * GENE_ROW_MM,
                ROLE_COLOR[role] if role else NEUTRAL_GENE,
                bp_per_mm,
                x_lo,
                x_hi,
            )
            arrow = "→" if tx.strand == "+" else "←"
            tail = f" ({ROLE_TEXT[role]})" if role else ""
            text_at(
                left_mm - 1.5,
                y + row * GENE_ROW_MM + GENE_BOX_Y_MM,
                f"{tx.symbol} {arrow}{tail}",
                ha="right",
                va="center",
            )
        y += gene_mm

    # Allele track: one tick per allele that separates haplotypes.
    if show_vars:
        ax_v = axes_at(y, VAR_MM)
        quiet = n_hap > 1 and not discr
        _draw_allele_track(
            ax_v,
            x_lo,
            x_hi,
            relevant,
            known_ids,
            color_of,
            plp,
            known_here,
            lettered,
            quiet,
        )
        if not quiet:
            text_at(
                left_mm - 1.5,
                y + VAR_MM / 2,
                "Alleles differing\nbetween haplotypes"
                if n_hap > 1
                else "Named alleles",
                ha="right",
                va="center",
                linespacing=1.1,
            )
    y += var_mm

    # One band per haplotype.
    guides = [_at_x(v.pos) for v in lettered]
    for index, hap in enumerate(names):
        hr, role = reads[hap], roles[hap]
        color = hap_color[hap]
        info = hap_info.get(hap, {})
        chosen, row_of, n_rows = plan[hap]
        reads_mm = reads_of[hap]
        band_mm = cov_mm + 0.8 + reads_mm
        pitch_mm = min(MAX_PITCH_MM, reads_mm / n_rows)
        pitch_pt = pitch_mm / PT_MM

        ax_c = axes_at(y, cov_mm)
        ax_r = axes_at(y + cov_mm + 0.8, reads_mm)
        for ax in (ax_c, ax_r):
            ax.set_xlim(x_lo, x_hi)
            for side in ("top", "right", "left", "bottom"):
                ax.spines[side].set_visible(False)
            ax.set_xticks([])
        ax_r.set_ylim(reads_mm / pitch_mm - 0.5, -0.5)
        ax_r.set_yticks([])
        ax_r.tick_params(left=False)
        if index % 2 == 0:
            ax_r.set_facecolor(BAND_COLOR)
        if hr.rows:
            grid, depth = coverage(hr, x_lo, x_hi)
            peak = max(int(depth.max()), 1)
            ax_c.fill_between(
                grid, 0, depth, color=tint(color, 0.55), lw=0, step="mid"
            )
            ax_c.set_ylim(0, peak * 1.05)
            if cov_mm >= 3.8:
                ax_c.set_yticks([peak])
                ax_c.tick_params(pad=1, length=1.5)
            else:
                ax_c.set_yticks([])
        else:
            ax_c.set_yticks([])

        carriers = [v for v in relevant if hap in v.haps]
        positions = [v.pos for v in carriers]
        _draw_reads(
            ax_r,
            chosen,
            row_of,
            carriers,
            positions,
            sep_ids,
            known_ids,
            color_of,
            pitch_pt,
            guides,
        )

        # Left label with the role strip.
        fig.add_artist(
            Rectangle(
                (0.8 / WIDTH_MM, 1 - (y + band_mm) / height_mm),
                1.3 / WIDTH_MM,
                band_mm / height_mm,
                transform=fig.transFigure,
                facecolor=color,
                edgecolor="none",
            )
        )
        n_total = len(hr.rows)
        shown = f" ({len(chosen)} shown)" if len(chosen) < n_total else ""
        lines = [
            (f"{role.symbol} ({ROLE_TEXT[role.role]})", "bold"),
            (hap, "normal"),
            (f"n = {n_total:,} reads{shown}", "normal"),
        ]
        note = _truncation_text(info.get("is_truncated") or ())
        if note:
            lines.append((note, "normal"))
        if info.get("two_copy"):
            lines.append(("2 identical copies (depth)", "normal"))
        for k, (line, weight) in enumerate(lines[: int(band_mm / LINE_MM)]):
            text_at(
                3.4,
                y + 1.0 + k * LINE_MM,
                line,
                ha="left",
                va="top",
                fontweight=weight,
            )
        y += band_mm + BAND_GAP_MM

    ax_x = axes_at(y - BAND_GAP_MM + 0.5, 0.1)
    ax_x.set_xlim(x_lo, x_hi)
    ax_x.set_yticks([])
    for side in ("top", "right", "left"):
        ax_x.spines[side].set_visible(False)
    # Drop a tick whose label would run past the right figure edge.
    ticks_bp = [
        t
        for t in MaxNLocator(6).tick_values(x_lo, x_hi)
        if x_lo <= t <= x_hi and (t - x_lo) / bp_per_mm <= ax_w - 5.0
    ]
    ax_x.set_xticks(ticks_bp)
    ax_x.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax_x.set_xlabel(f"{chrom} position (GRCh38, bp)", labelpad=2)
    y += AXIS_MM - BAND_GAP_MM + 0.5

    fig.legend(
        handles=handles,
        loc="upper left",
        ncol=len(handles) if legend_rows == 1 else math.ceil(len(handles) / 2),
        frameon=False,
        bbox_to_anchor=(2.0 / WIDTH_MM, 1 - y / height_mm),
        handletextpad=0.2,
        columnspacing=1.0,
        borderaxespad=0,
        borderpad=0,
    )
    y += legend_mm
    col_mm = (WIDTH_MM - 4.0) / 3
    for k, (letter, v) in enumerate(
        zip(string.ascii_uppercase, lettered, strict=False)
    ):
        label = v.label or f"{v.pos:,}"
        entry = f"{letter}  {label[:40]} ({v.pos:,})"
        text_at(
            2.0 + (k % 3) * col_mm,
            y + (k // 3) * KEY_ROW_MM,
            entry,
            ha="left",
            va="top",
        )
    fig.savefig(
        target,
        format="jpeg",
        dpi=FIGURE_DPI,
        pil_kwargs={"quality": 95, "dpi": (FIGURE_DPI, FIGURE_DPI)},
    )
    plt.close(fig)
    return target


@functools.cache
def _default_roles() -> dict[tuple[str, str], HapRole]:
    return load_roles()


# ------------------------------------------------------------ the driver


def _locus_symbols(
    locus: Mapping[str, object], roles: Mapping[str, HapRole]
) -> list[str]:
    """Gene symbols of the resolved haplotypes of a locus, sorted."""
    return sorted(
        {
            roles[h].symbol
            for h in locus["haps"]
            if roles[h].role != "unresolved"
        }
    )


def locus_file_names(
    region: str,
    loci: Sequence[Mapping[str, object]],
    roles: Mapping[str, HapRole],
) -> list[str]:
    """Name the figure of each locus of a region.

    One locus: ``<region>.jpeg``. Several loci: ``<region>_<SYMBOLS>.jpeg``
    from the role-table symbols of the haplotypes on that locus; if any
    locus has no symbol or two loci share one, ``<region>_locus<N>.jpeg``.
    """
    if len(loci) == 1:
        return [f"{region}.jpeg"]
    labels = [
        "-".join(_locus_symbols(lc, roles)).replace("/", "-") for lc in loci
    ]
    if all(labels) and len(set(labels)) == len(labels):
        return [f"{region}_{label}.jpeg" for label in labels]
    return [f"{region}_locus{n}.jpeg" for n in range(1, len(loci) + 1)]


def _region_task(
    region: str,
    hap_bams: Mapping[str, str],
    alleles: Sequence[Allele],
    known: Sequence[Allele],
    out_dir: str,
    options: Mapping[str, object],
) -> list[str]:
    """Draw every locus of one region (runs inside a worker process)."""
    reads = {h: load_reads(b) for h, b in hap_bams.items()}
    spans = {h: (r.chrom, *span(r)) for h, r in reads.items() if r.rows}
    loci = cluster_loci(spans)
    table = load_roles(options.get("roles_tsv"))
    roles = {h: resolve_role(h, table, region) for h in spans}
    refgene = options.get("gene_models")
    defaults = options.get("default_transcripts")
    default_map = read_default_transcripts(defaults) if defaults else None
    names = locus_file_names(region, loci, roles)
    symbols = [_locus_symbols(lc, roles) for lc in loci]
    written = []
    for number, locus in enumerate(loci):
        name = names[number]
        subtitle = ""
        if len(loci) > 1:
            others = ", ".join(n for k, n in enumerate(names) if k != number)
            subtitle = (
                f"{region}: {', '.join(symbols[number]) or 'unresolved'}"
                f" \u00b7 paired figure: {others}"
            )
        transcripts = (
            read_gene_models(
                refgene,
                default_map,
                locus["chrom"],
                locus["start0"],
                locus["end0"],
                prefer=[roles[h].symbol for h in locus["haps"]],
            )
            if refgene
            else []
        )
        draw_locus(
            locus,
            reads,
            alleles,
            known,
            Path(out_dir) / name,
            hap_roles=roles,
            hap_info=options.get("hap_info"),
            transcripts=transcripts,
            max_reads=options.get("max_reads", MAX_READS),
            subtitle=subtitle,
        )
        written.append(str(Path(out_dir) / name))
    return written


def _init_worker() -> None:
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["MPLBACKEND"] = "Agg"


def default_region(hap: str) -> str:
    """Return the region name of an HP value (text before the last ``_``)."""
    return hap.rsplit("_", 1)[0] if "_" in hap else hap


def plot_regions(
    merged_vcf: Path | str,
    hap_dirs: Mapping[str, Path | str],
    pathogenic_keys: set[str],
    pseudogene_tsv: Path | str,
    reference: Path | str,
    out_dir: Path | str,
    workers: int = 1,
    region_of: Mapping[str, str] | None = None,
    allele_labels: Mapping[str, str] | None = None,
    *,
    gene_models: Path | str | None = None,
    default_transcripts: Path | str | None = None,
    hap_info: Mapping[str, dict] | Path | str | None = None,
    haplotype_roles: Path | str | None = None,
    max_reads: int = MAX_READS,
    hap_symbols: Mapping[str, str] | None = None,
) -> list[Path]:
    """Write one read-pileup JPEG per region (or per locus).

    Args:
        merged_vcf: Merged VCF; samples are HP values; haploid GT 1/0/.
        hap_dirs: HP -> directory holding ``reads.bam``.
        pathogenic_keys: ``allele_key`` strings classed P/LP.
        pseudogene_tsv: known_gene_pseudogene_variants.tsv.
        reference: GRCh38 FASTA; its ``.fai`` must list every contig of
            the VCF (checked).
        out_dir: Output directory (created).
        workers: Processes over regions (1 = run in this process).
        region_of: Optional HP -> region name; default is the HP text
            before its last underscore (``GBA_gene1hap2`` -> ``GBA``).
        allele_labels: Optional allele key -> label for P/LP alleles;
            the position is the fallback label.
        gene_models: Optional UCSC refGene table; draws the gene track.
        default_transcripts: Optional DefaultTranscript.tsv choosing the
            transcript drawn per gene.
        hap_info: Optional HP -> dict, or the path of the Paraphase JSON
            (see ``load_hap_info``); adds truncation and copy notes.
        haplotype_roles: Optional replacement for the packaged
            haplotype_roles.tsv.
        max_reads: Reads drawn per haplotype (seeded random sample; the
            label states the true number).
        hap_symbols: Optional HP -> gene symbol for known-variant markers
            where the role table has no rows (RCCX).

    Returns:
        Written JPEG paths, sorted by name.

    Raises:
        ValueError: If a VCF contig is missing from the reference index.
        FileNotFoundError: If a haplotype has no reads.bam.
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    with pysam.FastaFile(str(reference)) as fasta:
        contigs = set(fasta.references)
    with pysam.VariantFile(str(merged_vcf)) as vcf:
        missing = set(vcf.header.contigs) - contigs
    if missing:
        raise ValueError(f"VCF contigs not in reference: {sorted(missing)}")
    alleles = read_alleles(merged_vcf, pathogenic_keys, allele_labels)
    hap_region = {
        h: region_of[h] if region_of else default_region(h) for h in hap_dirs
    }
    known_rows = read_known(pseudogene_tsv)
    known = [
        Allele(
            a.chrom,
            a.pos,
            a.pos + len(a.ref) - 1,
            a.ref,
            a.alt,
            False,
            haps,
            label=row["label"],
        )
        for row, a, haps in _known_carriers(
            known_rows,
            alleles,
            hap_region,
            load_roles(haplotype_roles) if haplotype_roles else None,
            hap_symbols,
        )
    ]
    if isinstance(hap_info, (str, Path)):
        hap_info = load_hap_info(hap_info)
    by_region: dict[str, dict[str, str]] = collections.defaultdict(dict)
    for hap, folder in hap_dirs.items():
        bam = Path(folder) / "reads.bam"
        if not bam.exists():
            raise FileNotFoundError(bam)
        name = region_of[hap] if region_of else default_region(hap)
        by_region[name][hap] = str(bam)
    jobs = []
    for region, bams in sorted(by_region.items()):
        mine = [v for v in alleles if v.haps & set(bams)]
        options = {
            "gene_models": str(gene_models) if gene_models else None,
            "default_transcripts": (
                str(default_transcripts) if default_transcripts else None
            ),
            "hap_info": {
                h: i for h, i in (hap_info or {}).items() if h in bams
            },
            "roles_tsv": str(haplotype_roles) if haplotype_roles else None,
            "max_reads": max_reads,
        }
        jobs.append((region, bams, mine, known, str(out), options))
    if workers <= 1 or len(jobs) <= 1:
        results = [_region_task(*job) for job in jobs]
    else:
        with ProcessPoolExecutor(
            max_workers=workers, initializer=_init_worker
        ) as pool:
            futures = [pool.submit(_region_task, *job) for job in jobs]
            results = [f.result() for f in futures]
    return sorted(Path(p) for paths in results for p in paths)

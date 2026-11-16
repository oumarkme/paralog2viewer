# SPDX-License-Identifier: GPL-3.0-or-later
"""Draw read-pileup figures per region or locus.

One JPEG per region (or per locus when the haplotypes of a region are
projected to different loci). One row per read, grouped by ``HP`` tag.
Coordinates: VCF positions are 1-based; read spans from pysam are 0-based
half-open; the x axis is in GRCh38 bp.

Port of the former Nextflow process ``PLOT_REGION``.
"""

from __future__ import annotations

import bisect
import collections
import csv
import os
from collections.abc import Mapping, Sequence
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import pysam  # noqa: E402
from matplotlib import pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import FuncFormatter  # noqa: E402

from paralog2viewer.datalayout import allele_key  # noqa: E402
from paralog2viewer.figstyle import FONT_PT, base_rc  # noqa: E402

# Okabe-Ito (Okabe and Ito, 2008) without black; black marks overflow.
PALETTE = (
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7",
)
READ_COLOR = "#9A9A9A"
BAND_COLOR = "#F2F2F2"
LINE_COLOR = "#4D4D4D"
MM = 25.4
WIDTH_MM = 180.0
MAX_HEIGHT_MM = 247.0
DOT_PT = 2.6
MIN_SV_BP = 30
GAP_ROWS = 2


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
        label: Text shown above the plot for P/LP and known alleles.
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


def read_known(path: Path | str) -> list[Allele]:
    """Read pseudogene-type alleles from known_gene_pseudogene_variants.tsv.

    Args:
        path: Tab-separated file; ``#`` lines are comments; needs columns
            chrom, pos (1-based), ref, alt, label.

    Returns:
        One Allele per row, with ``haps`` empty (set per locus later).
    """
    known = []
    with open(path, newline="", encoding="utf-8") as handle:
        lines = (ln for ln in handle if not ln.startswith("#"))
        for row in csv.DictReader(lines, delimiter="\t"):
            pos = int(row["pos"])
            known.append(
                Allele(
                    row["chrom"],
                    pos,
                    pos,
                    row["ref"],
                    row["alt"],
                    False,
                    frozenset(),
                    label=row["label"],
                )
            )
    return known


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
    rows.sort(
        key=lambda g: (
            min(a.reference_start for a in g),
            max(a.reference_end for a in g),
        )
    )
    return HapReads(rows, chrom)


def span(reads: HapReads) -> tuple[int, int]:
    """Return the 0-based half-open reference span covered by the reads."""
    return (
        min(a.reference_start for g in reads.rows for a in g),
        max(a.reference_end for g in reads.rows for a in g),
    )


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
            "ytick.major.size": 0,
            "figure.dpi": 300,
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


def _scan(
    reads: HapReads,
    carriers: list[Allele],
    row0: int,
    hits: dict,
) -> tuple[list, int]:
    """Find, for each read row, the carried alleles it supports.

    Args:
        reads: Haplotype reads.
        carriers: Alleles called in this haplotype, sorted by pos.
        row0: First y row of this haplotype.
        hits: Accumulator, allele -> list of y rows with the ALT.

    Returns:
        Line segments for the reads and the next free row.
    """
    segments = []
    positions = [v.pos for v in carriers]
    row = row0
    for alignments in reads.rows:
        for a in alignments:
            segments.append([(a.reference_start, row), (a.reference_end, row)])
        lo = min(a.reference_start for a in alignments)
        hi = max(a.reference_end for a in alignments)
        maps: dict[int, dict[int, int]] = {}
        for v in carriers[
            bisect.bisect_right(positions, lo) : bisect.bisect_right(
                positions, hi
            )
        ]:
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
                                for q, r in a.get_aligned_pairs(
                                    matches_only=True
                                )
                            }
                        state = alt_state(a, maps[id(a)], s0, v.ref, v.alt)
                        break
            if state:
                hits.setdefault(id(v), []).append(row)
        row += 1
    return segments, row


def draw_locus(
    locus: dict,
    reads: Mapping[str, HapReads],
    alleles: Sequence[Allele],
    known: Sequence[Allele],
    target: Path,
) -> Path:
    """Render one locus to a JPEG (300 dpi, 180 mm wide).

    Args:
        locus: Dict from ``cluster_loci``.
        reads: HP -> reads.
        alleles: Called alleles (any locus; filtered here).
        known: Pseudogene-type alleles from the TSV.
        target: Output JPEG path.

    Returns:
        ``target``.
    """
    _set_style()
    chrom, x_lo, x_hi = locus["chrom"], locus["start0"], locus["end0"]
    names = sorted(locus["haps"])
    nameset = set(names)
    here = [
        v
        for v in alleles
        if v.chrom == chrom and x_lo < v.pos <= x_hi and v.haps & nameset
    ]
    # Only alleles called in some, not all, haplotypes separate them.
    discr = [v for v in here if not nameset <= v.haps] if len(names) > 1 else []
    sep_ids = {id(v) for v in discr}
    counts = collections.Counter(tuple(sorted(v.haps & nameset)) for v in discr)
    patterns = [p for p, _ in counts.most_common()]
    color_of = {
        p: PALETTE[i] if i < len(PALETTE) else "#000000"
        for i, p in enumerate(patterns)
    }
    plp = [v for v in here if v.plp]
    known_here = [
        Allele(
            k.chrom,
            k.pos,
            k.end,
            k.ref,
            k.alt,
            False,
            frozenset(names),
            label=k.label,
        )
        for k in known
        if k.chrom == chrom and x_lo < k.pos <= x_hi
    ]
    known_ids = {id(k) for k in known_here}
    relevant = sorted(
        discr + [v for v in plp if id(v) not in sep_ids] + known_here,
        key=lambda v: v.pos,
    )
    plp_label = {id(v): v.label or f"{v.pos:,}" for v in plp}

    total_rows = sum(len(reads[h].rows) for h in names) + GAP_ROWS * (
        len(names) - 1
    )
    label_items = [(v.pos, plp_label[id(v)]) for v in plp] + [
        (k.pos, k.label) for k in known_here
    ]
    longest = min(30, max([len(t) for _, t in label_items] or [0]))
    pt_mm = 0.3528
    top_mm = 4.5 + (0.5 * FONT_PT * pt_mm * longest if longest else 0.0)
    left_mm = 4.0 + 0.6 * FONT_PT * pt_mm * max(len(h) for h in names)
    bottom_mm, right_mm = 9.0, 3.0
    has_sv = any(v.symbolic for v in discr)
    legend_entries = (
        min(len(patterns), len(PALETTE) + 1)
        + has_sv
        + bool(plp)
        + bool(known_here)
        + 1
    )
    legend_mm = 5.0 + 3.0 * ((legend_entries + 2) // 3)
    max_axes_mm = MAX_HEIGHT_MM - top_mm - bottom_mm - legend_mm
    pitch_pt = min(2.6, max_axes_mm / pt_mm / max(total_rows, 1))
    dot_pt = min(DOT_PT, 1.15 * pitch_pt)
    axes_mm = max(total_rows * pitch_pt * pt_mm, 20.0)
    height_mm = top_mm + axes_mm + bottom_mm + legend_mm
    fig = plt.figure(figsize=(WIDTH_MM / MM, height_mm / MM))
    ax = fig.add_axes(
        [
            left_mm / WIDTH_MM,
            (bottom_mm + legend_mm) / height_mm,
            (WIDTH_MM - left_mm - right_mm) / WIDTH_MM,
            axes_mm / height_mm,
        ]
    )
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(total_rows + 1, -1)

    row = 0
    centers = []
    hits: dict[int, list[int]] = {}
    for index, hap in enumerate(names):
        first = row
        carriers = [v for v in relevant if hap in v.haps]
        segments, row = _scan(reads[hap], carriers, row, hits)
        ax.add_collection(
            LineCollection(
                segments,
                colors=READ_COLOR,
                linewidths=0.6,
                zorder=2,
                capstyle="butt",
            )
        )
        if index % 2 == 0:
            ax.axhspan(first - 0.5, row - 0.5, color=BAND_COLOR, zorder=0, lw=0)
        centers.append((first + row - 1) / 2)
        row += GAP_ROWS

    # Per-read marks. Each allele's rows are only valid for the haplotypes
    # that were scanned, so rows are rebuilt from the per-allele hit lists.
    dots = collections.defaultdict(lambda: ([], []))
    sv_dots = collections.defaultdict(lambda: ([], []))
    known_xy: tuple[list, list] = ([], [])
    plp_xy: tuple[list, list] = ([], [])
    for v in relevant:
        rows = hits.get(id(v), [])
        if id(v) in known_ids:
            known_xy[0].extend([v.pos] * len(rows))
            known_xy[1].extend(rows)
            continue
        if id(v) in sep_ids:
            pattern = tuple(sorted(v.haps & nameset))
            out = sv_dots if v.symbolic else dots
            out[pattern][0].extend([v.pos] * len(rows))
            out[pattern][1].extend(rows)
        if v.plp:
            plp_xy[0].extend([v.pos] * len(rows))
            plp_xy[1].extend(rows)
    for markers, shape in ((dots, "o"), (sv_dots, "v")):
        for pattern, (xs, ys) in markers.items():
            ax.scatter(
                xs,
                ys,
                s=dot_pt**2,
                marker=shape,
                c=color_of[pattern],
                edgecolors="black",
                linewidths=0.25,
                zorder=4,
            )
    if known_xy[0]:
        ax.scatter(
            *known_xy,
            s=(1.15 * dot_pt) ** 2,
            marker="s",
            c="black",
            edgecolors="white",
            linewidths=0.3,
            zorder=6,
        )
    if plp_xy[0]:
        ax.scatter(
            *plp_xy,
            s=(1.3 * dot_pt) ** 2,
            marker="D",
            c="white",
            edgecolors="black",
            linewidths=0.45,
            zorder=5,
        )
    ax.set_yticks(centers)
    ax.set_yticklabels([f"{h}\n(n={len(reads[h].rows)})" for h in names])
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:,.0f}"))
    ax.set_xlabel(f"{chrom} position (GRCh38, bp)")
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)

    # Label above the plot, left to right; skip only on overlap.
    min_gap_bp = (
        (x_hi - x_lo)
        / ((WIDTH_MM - left_mm - right_mm) / MM * 72)
        * (FONT_PT * 1.3)
    )
    placed: list[int] = []
    unlabeled = 0
    for position, text in sorted(label_items):
        if placed and position - placed[-1] < min_gap_bp:
            if position != placed[-1]:
                unlabeled += 1
            continue
        placed.append(position)
        ax.annotate(
            text[:30],
            xy=(position, -1),
            xytext=(position, -1.6),
            rotation=90,
            ha="center",
            va="bottom",
            fontsize=FONT_PT,
            annotation_clip=False,
            color="black",
        )
        ax.plot(
            [position] * 2,
            [-1.5, total_rows + 1],
            lw=0.3,
            ls=(0, (2, 2)),
            color=LINE_COLOR,
            zorder=1,
            clip_on=False,
        )

    handles = [
        _marker("o", DOT_PT, color_of[p], "black", 0.25, " + ".join(p))
        for p in patterns[: len(PALETTE)]
    ]
    if len(patterns) > len(PALETTE):
        handles.append(
            _marker("o", DOT_PT, "#000000", "black", 0.25, "other combinations")
        )
    if sv_dots:
        handles.append(
            _marker("v", DOT_PT, "white", "black", 0.25, "structural allele")
        )
    if plp:
        handles.append(
            _marker(
                "D",
                1.3 * DOT_PT,
                "white",
                "black",
                0.6,
                "pathogenic or likely pathogenic (ClinVar or InterVar)",
            )
        )
    if known_here:
        handles.append(
            _marker(
                "s",
                1.15 * DOT_PT,
                "black",
                "white",
                0.3,
                "pseudogene-type allele at a published gene-vs-pseudogene site",
            )
        )
    if unlabeled:
        handles.append(
            Line2D(
                [],
                [],
                ls="",
                label=f"{unlabeled} more such alleles not labeled",
            )
        )
    fig.legend(
        handles=handles,
        loc="lower left",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(left_mm / WIDTH_MM, 0.0),
        handletextpad=0.2,
        columnspacing=1.2,
        alignment="left",
        title_fontsize=FONT_PT,
        title="Colored dots: alleles called in only some haplotypes "
        "(color = haplotypes calling them)",
    )
    fig.savefig(
        target,
        format="jpeg",
        dpi=300,
        pil_kwargs={"quality": 95, "dpi": (300, 300)},
    )
    plt.close(fig)
    return target


# ------------------------------------------------------------ the driver


def _region_task(
    region: str,
    hap_bams: Mapping[str, str],
    alleles: Sequence[Allele],
    known: Sequence[Allele],
    out_dir: str,
) -> list[str]:
    """Draw every locus of one region (runs inside a worker process)."""
    reads = {h: load_reads(b) for h, b in hap_bams.items()}
    spans = {h: (r.chrom, *span(r)) for h, r in reads.items() if r.rows}
    loci = cluster_loci(spans)
    written = []
    for number, locus in enumerate(loci, 1):
        name = (
            f"{region}.jpeg"
            if len(loci) == 1
            else f"{region}_locus{number}.jpeg"
        )
        draw_locus(locus, reads, alleles, known, Path(out_dir) / name)
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
    known = read_known(pseudogene_tsv)
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
        jobs.append((region, bams, mine, known, str(out)))
    if workers <= 1 or len(jobs) <= 1:
        results = [_region_task(*job) for job in jobs]
    else:
        with ProcessPoolExecutor(
            max_workers=workers, initializer=_init_worker
        ) as pool:
            futures = [pool.submit(_region_task, *job) for job in jobs]
            results = [f.result() for f in futures]
    return sorted(Path(p) for paths in results for p in paths)

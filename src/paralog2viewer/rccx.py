# SPDX-License-Identifier: GPL-3.0-or-later
"""RCCX module calls and the module and gene-conversion figure.

Three parts:

* ``build_rccx_sites`` (one time, from the reference only): BLAST of each
  distal RCCX gene against its proximal paralog; every mismatching aligned
  column becomes a discriminating site.
* ``classify_rccx`` (per run): module structure per allele, CYP21A2 versus
  CYP21A1P from the Paraphase ``hap_variants`` marker hits, and for the RP,
  C4, TNX (and, as a second opinion, CYP21) pairs a vote over the
  discriminating sites using the haploid genotypes of the merged VCF.
* ``plot_rccx_module``: module schematic per allele plus the per-haplotype
  call table (the read pileup is the standard figure of ``plot.py``).

Coordinates: the gene table and the sites table are 1-based and inclusive
(GRCh38 chr6, VCF-like); Paraphase ``boundary`` values are 1-based
positions. Paraphase aligns every haplotype to the distal (reference) copy
of each gene, but the pipeline realigns each read to its best locus. A
proximal-type haplotype therefore carries its alleles in distal-gene
coordinates (ALT at the site) or, when its reads sit in the proximal copy,
at ``proximal_pos`` (reference there). Each site is voted in the frame
where the haplotype has the reads (see ``_frame_vote``).

Clean rewrite of ``test/fig_and_tab.py`` of the former RCCX project. Vote
thresholds are unchanged: >= 90 % proximal sites is proximal type, <= 10 %
is distal type, chimeric needs opposite majorities with a margin of at
least 0.65 and at least 5 informative sites in each half, fewer than 3
sites is not determined. Anything else is mixed for RP and CYP21. For C4
and TNX it is a leaning call: "<gene>-like" for the majority type with
the share of covered sites matching it (above 50 %, below 90 %), and
mixed only at exactly 50 %.
"""

from __future__ import annotations

import csv
import json
import os
import tempfile
from collections import defaultdict
from collections.abc import Callable, Mapping, Sequence
from dataclasses import astuple, dataclass
from importlib import resources
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import pysam  # noqa: E402
from matplotlib import pyplot as plt  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402

from paralog2viewer import proc  # noqa: E402
from paralog2viewer.figstyle import FIGURE_DPI, FONT_PT, base_rc  # noqa: E402

BLASTN_PARAMS = ("-task", "blastn", "-dust", "no", "-evalue", "1e-20")
SITE_COLUMNS = (
    "family",
    "chrom",
    "pos",
    "distal_gene",
    "proximal_gene",
    "distal_base",
    "proximal_base",
    "proximal_pos",
)
MIN_SITES = 3
PROXIMAL_FRACTION = 0.9
DISTAL_FRACTION = 0.1
# Alignments below this identity are repeat or self-similar hits, not the
# orthologous gene copy; their mismatches are not fixed paralog differences.
MIN_HSP_IDENTITY = 0.95
CHIMERA_MIN_MARGIN = 0.65
CHIMERA_MIN_HALF_SITES = 5
DEFAULT_MIN_DEPTH = 3
# Families that report a leaning call ("C4A-like 78 %") instead of mixed.
LIKE_FAMILIES = ("RP", "C4", "TNX")
# Haplotypes of a standard bimodular allele; fewer means deleted modules.
STANDARD_MODULES = 2
STRUCTURE_NAMES = {1: "monomodular", 2: "bimodular", 3: "trimodular"}

# Okabe-Ito (Okabe and Ito, 2008), one color per gene.
GENE_COLORS = {
    "WHR1": "#E69F00",
    "WHR1B": "#F0E442",
    "C4A": "#0072B2",
    "C4B": "#56B4E9",
    "CYP21A1P": "#D55E00",
    "CYP21A2": "#009E73",
    "TNXA": "#CC79A7",
    "TNXB": "#000000",
}
LIKE_LINEWIDTH = 1.0  # twice the normal box edge
LIKE_LEGEND_EDGE = "#555555"  # gray: a gene color would name a gene
UNRESOLVED_FILL = {"chimeric": "#FFFFFF", "mixed": "#CFCFCF", "n.d.": "#FFFFFF"}
# The conflict flag is a dotted hatch on white, not a color: all eight
# Okabe-Ito colors already name genes, so any hue would read as a gene.
CONFLICT_FILL = "#FFFFFF"
CONFLICT_HATCH = "...."
HEAD_FILL = "#DDDDDD"
BAND_FILL = "#F2F2F2"
EDGE = "#000000"
MM = 25.4
WIDTH_MM = 180.0
MAX_HEIGHT_MM = 247.0


# ------------------------------------------------------------- gene table


@dataclass(frozen=True)
class Gene:
    """One gene box (GRCh38, 1-based inclusive)."""

    name: str
    family: str
    role: str
    chrom: str
    start: int
    end: int
    display: str


@dataclass(frozen=True)
class GeneTable:
    """Gene boxes, gene pairs and CYP21 marker coordinates."""

    genes: tuple[Gene, ...]
    marker_window: tuple[str, int, int] | None
    promoter_sites: tuple[int, ...]

    def families(self) -> list[str]:
        """Return the families in 5' to 3' order of first appearance."""
        return list(dict.fromkeys(g.family for g in self.genes))

    def pair(self, family: str) -> tuple[Gene, Gene]:
        """Return (proximal, distal) genes of a family.

        Raises:
            ValueError: If the family does not have exactly one of each.
        """
        prox = [
            g for g in self.genes if g.family == family and g.role == "proximal"
        ]
        dist = [
            g for g in self.genes if g.family == family and g.role == "distal"
        ]
        if len(prox) != 1 or len(dist) != 1:
            raise ValueError(
                f"family {family!r} needs one proximal and one distal gene"
            )
        return prox[0], dist[0]

    def gene(self, name: str) -> Gene:
        """Return a gene by name."""
        for gene in self.genes:
            if gene.name == name:
                return gene
        raise KeyError(name)


def default_genes_tsv() -> Path:
    """Return the path of the packaged ``rccx_genes.tsv``."""
    return Path(
        str(resources.files("paralog2viewer") / "resources" / "rccx_genes.tsv")
    )


def _data_lines(path: Path) -> list[str]:
    with open(path, encoding="utf-8") as handle:
        return [ln for ln in handle if ln.strip() and not ln.startswith("#")]


def load_gene_table(path: Path | None = None) -> GeneTable:
    """Read the gene table.

    Args:
        path: TSV with columns kind, name, family, role, chrom, start, end,
            display (1-based inclusive); the packaged table when None.

    Returns:
        The parsed table.

    Raises:
        ValueError: If a row is malformed or start exceeds end.
    """
    genes: list[Gene] = []
    window = None
    promoters: list[int] = []
    rows = csv.DictReader(
        _data_lines(path or default_genes_tsv()), delimiter="\t"
    )
    for row in rows:
        start, end = int(row["start"]), int(row["end"])
        if start < 1 or end < start:
            raise ValueError(f"bad interval in gene table: {row}")
        kind = row["kind"]
        if kind == "gene":
            genes.append(
                Gene(
                    row["name"],
                    row["family"],
                    row["role"],
                    row["chrom"],
                    start,
                    end,
                    row["display"] or row["name"],
                )
            )
        elif kind == "marker_window":
            window = (row["chrom"], start, end)
        elif kind == "promoter_site":
            promoters.append(start)
        else:
            raise ValueError(f"unknown kind {kind!r} in gene table")
    return GeneTable(tuple(genes), window, tuple(promoters))


# ------------------------------------------------------------- precompute


@dataclass(frozen=True)
class Site:
    """One discriminating site (1-based).

    ``pos`` is in distal-gene coordinates. ``proximal_pos`` is the aligned
    position in the proximal gene (0 when unknown, e.g. a table built
    before this column existed, or a reverse-strand alignment); reads
    that were realigned to the proximal copy carry the site there.
    """

    family: str
    chrom: str
    pos: int
    distal_gene: str
    proximal_gene: str
    distal_base: str
    proximal_base: str
    proximal_pos: int = 0


def _blast_mismatches(
    blast_text: str,
    query_start: int,
    subject_start: int = 0,
    min_identity: float = 0.0,
) -> dict[int, tuple[str, str, int]]:
    """Return {query pos: (query base, subject base, subject pos)}.

    Args:
        blast_text: blastn output, ``-outfmt "6 qstart qend sstart send
            qseq sseq"``; sseq is in query orientation.
        query_start: 1-based genomic position of the first query base.
        subject_start: 1-based genomic position of the first subject
            base; 0 leaves the subject position unknown (0).
        min_identity: alignments whose fraction of identical columns
            (over all alignment columns) is below this are skipped.

    Returns:
        Mismatching aligned columns without gaps or ambiguity codes. If
        several alignments cover a position, the first reported one wins.
        The subject position is 0 for a reverse-strand alignment.
    """
    diffs: dict[int, tuple[str, str, int]] = {}
    for line in blast_text.splitlines():
        if not line.strip():
            continue
        qstart, _qend, sstart, send, qseq, sseq = line.split("\t")
        matches = sum(
            a == b for a, b in zip(qseq.upper(), sseq.upper(), strict=True)
        )
        if matches / len(qseq) < min_identity:
            continue
        qpos = query_start + int(qstart) - 1
        forward = int(sstart) <= int(send)
        spos = subject_start + int(sstart) - 1 if subject_start else 0
        for qbase, sbase in zip(qseq.upper(), sseq.upper(), strict=True):
            if qbase == "-":
                spos += 1 if forward else 0
                continue
            if (
                qbase != sbase
                and sbase != "-"
                and qbase in "ACGT"
                and sbase in "ACGT"
            ):
                known = forward and subject_start > 0
                diffs.setdefault(qpos, (qbase, sbase, spos if known else 0))
            qpos += 1
            spos += 1 if sbase != "-" and forward else 0
    return diffs


def build_rccx_sites(
    reference: Path,
    genes_tsv: Path,
    out_tsv: Path,
    *,
    threads: int = 1,
    blastn: str = "blastn",
) -> int:
    """Write the RCCX discriminating sites.

    For every gene family, the distal gene span is the BLAST query and the
    proximal gene span the subject (``-task blastn -dust no -evalue
    1e-20``); each mismatching aligned column is one site.

    Args:
        reference: GRCh38 FASTA with a ``.fai`` index.
        genes_tsv: Gene table (see ``load_gene_table``).
        out_tsv: Output TSV with ``SITE_COLUMNS``, sorted by family order
            and position; written atomically.
        threads: Threads for blastn.
        blastn: blastn executable.

    Returns:
        Number of sites written.

    Raises:
        ValueError: If a gene lies outside its contig in the reference.
        proc.ToolError: If blastn fails.
    """
    table = load_gene_table(genes_tsv)
    sites: list[Site] = []
    with (
        pysam.FastaFile(str(reference)) as fasta,
        tempfile.TemporaryDirectory(prefix="rccx_blast_") as tmp,
    ):
        lengths = dict(zip(fasta.references, fasta.lengths, strict=True))
        for family in table.families():
            prox, dist = table.pair(family)
            for gene in (prox, dist):
                if gene.chrom not in lengths or gene.end > lengths[gene.chrom]:
                    raise ValueError(
                        f"{gene.name} {gene.chrom}:{gene.start}-"
                        f"{gene.end} is outside the reference"
                    )
            query, subject = Path(tmp, "q.fa"), Path(tmp, "s.fa")
            for path, gene in ((query, dist), (subject, prox)):
                seq = fasta.fetch(gene.chrom, gene.start - 1, gene.end)
                path.write_text(f">{gene.name}\n{seq.upper()}\n")
            result = Path(tmp, f"{family}.blast.tsv")
            proc.run(
                [
                    blastn,
                    "-query",
                    query,
                    "-subject",
                    subject,
                    *BLASTN_PARAMS,
                    "-num_threads",
                    str(threads),
                    "-outfmt",
                    "6 qstart qend sstart send qseq sseq",
                ],
                stdout=result,
            )
            diffs = _blast_mismatches(
                result.read_text(),
                dist.start,
                prox.start,
                min_identity=MIN_HSP_IDENTITY,
            )
            sites += [
                Site(
                    family,
                    dist.chrom,
                    pos,
                    dist.name,
                    prox.name,
                    dbase,
                    pbase,
                    ppos,
                )
                for pos, (dbase, pbase, ppos) in sorted(diffs.items())
            ]
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    partial = out_tsv.with_name(out_tsv.name + ".partial")
    with open(partial, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(SITE_COLUMNS)
        writer.writerows(astuple(site) for site in sites)
    os.replace(partial, out_tsv)
    return len(sites)


def load_sites(path: Path) -> dict[str, list[Site]]:
    """Read the sites TSV into {family: sites sorted by position}."""
    out: dict[str, list[Site]] = defaultdict(list)
    with open(path, encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            out[row["family"]].append(
                Site(
                    row["family"],
                    row["chrom"],
                    int(row["pos"]),
                    row["distal_gene"],
                    row["proximal_gene"],
                    row["distal_base"],
                    row["proximal_base"],
                    int(row.get("proximal_pos") or 0),
                )
            )
    for sites in out.values():
        sites.sort(key=lambda s: s.pos)
    return dict(out)


# ------------------------------------------------------------------ votes


def _vcf_genotypes(
    vcf_path: Path, haps: Sequence[str], chrom: str, lo: int, hi: int
) -> tuple[dict[str, dict[int, list[tuple[str, str, str]]]], set[str]]:
    """Return ({hap: {pos: [(ref, alt, gt)]}}, VCF sample names).

    Records are those of ``chrom`` at 1-based positions in [lo, hi]; gt is
    "1", "0" or "." from the first genotype allele (haploid). Haplotypes
    absent from the VCF get an empty dict.
    """
    out: dict[str, dict[int, list[tuple[str, str, str]]]] = {
        h: defaultdict(list) for h in haps
    }
    with pysam.VariantFile(str(vcf_path)) as vcf:
        samples = list(vcf.header.samples)
        if len(set(samples)) != len(samples):
            raise ValueError("duplicate sample names in the merged VCF")
        present = [h for h in haps if h in samples]
        for rec in vcf:
            if rec.chrom != chrom or not lo <= rec.pos <= hi:
                continue
            if not rec.alts:
                continue
            if len(rec.alts) > 1:
                raise ValueError(
                    f"multi-allelic record at {rec.chrom}:{rec.pos}; the "
                    "merged VCF must have one record per allele"
                )
            for hap in present:
                call = rec.samples[hap]["GT"]
                if not call or call[0] is None:
                    gt = "."
                else:
                    gt = "1" if call[0] > 0 else "0"
                out[hap][rec.pos].append((rec.ref, rec.alts[0], gt))
    return {h: dict(v) for h, v in out.items()}, set(samples)


def _depths(
    bam_path: Path, chrom: str, positions: Sequence[int]
) -> dict[int, int]:
    """Return {1-based position: read depth} from the BAM."""
    if not positions:
        return {}
    lo, hi = min(positions), max(positions)
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        counts = bam.count_coverage(
            chrom, lo - 1, hi, quality_threshold=0, read_callback="all"
        )
    depth = [sum(col) for col in zip(*counts, strict=True)]
    return {p: depth[p - lo] for p in positions}


Records = list[tuple[str, str, str]] | None


def _vote(
    records: Records,
    elsewhere: bool,
    covered: bool | None,
    ref_vote: str,
    alt_vote: str,
    ref_base: str,
    alt_base: str,
) -> str:
    """Return the vote of one haplotype in one gene-copy frame.

    ``ref_vote`` is the vote of the reference allele in that frame and
    ``alt_vote`` that of the single expected REF/ALT change. A site without
    a VCF record counts as reference only if it is provably callable: the
    haplotype BAM has the minimum depth there (``covered``), or, without a
    BAM, another haplotype has a record there (``elsewhere``). Otherwise it
    is missing.
    """
    if not records:
        seen = covered if covered is not None else elsewhere
        return ref_vote if seen else "missing"
    gts = {gt for _ref, _alt, gt in records}
    if "1" in gts:
        hit = any(
            gt == "1" and ref == ref_base and alt == alt_base
            for ref, alt, gt in records
        )
        return alt_vote if hit else "other"
    return ref_vote if "0" in gts else "missing"


def _site_vote(
    site: Site, records: Records, callable_elsewhere: bool, covered: bool | None
) -> str:
    """Return proximal, distal, other or missing in the distal frame."""
    return _vote(
        records,
        callable_elsewhere,
        covered,
        "distal",
        "proximal",
        site.distal_base,
        site.proximal_base,
    )


def _proximal_frame_vote(
    site: Site, records: Records, callable_elsewhere: bool, covered: bool | None
) -> str:
    """Return the vote of a read realigned to the proximal gene copy.

    Mirror of ``_site_vote`` at ``site.proximal_pos``: there the proximal
    base is the reference, so reference is a proximal-type vote and the
    distal base is the ALT.
    """
    return _vote(
        records,
        callable_elsewhere,
        covered,
        "proximal",
        "distal",
        site.proximal_base,
        site.distal_base,
    )


def _has_call(records: Records) -> bool:
    return bool(records) and any(gt in ("0", "1") for _r, _a, gt in records)


def _frame_vote(
    site: Site,
    dist: tuple[Records, bool, int | None],
    prox: tuple[Records, bool, int | None],
    min_depth: int,
    use_bam: bool,
) -> str:
    """Vote one site in the frame (distal or proximal copy) of the reads.

    Paraphase aligns every haplotype to the distal gene, but the pipeline
    realigns each read to its best locus, so a proximal-type haplotype
    carries its alleles in proximal-gene coordinates. Each argument is
    (records, called by another haplotype, depth in the haplotype BAM or
    None); ``prox`` is ignored when ``site.proximal_pos`` is unknown.

    With a BAM the frame with the larger depth (at least ``min_depth``)
    decides. Without one, the frame where the haplotype has a genotype
    decides; a site without any record is reference in the single frame
    where another haplotype was called, otherwise missing.
    """
    d_rec, d_else, d_depth = dist
    p_rec, p_else, p_depth = prox
    in_distal = _site_vote(
        site, d_rec, d_else, None if d_depth is None else d_depth >= min_depth
    )
    if site.proximal_pos <= 0:
        return in_distal
    in_proximal = _proximal_frame_vote(
        site, p_rec, p_else, None if p_depth is None else p_depth >= min_depth
    )
    if use_bam:
        d_depth, p_depth = d_depth or 0, p_depth or 0
        if p_depth >= min_depth and (d_depth < min_depth or p_depth > d_depth):
            return in_proximal
        return in_distal
    if _has_call(d_rec):
        return in_distal
    if _has_call(p_rec):
        return in_proximal
    if d_else and not p_else:
        return "distal"
    if p_else and not d_else:
        return "proximal"
    return "missing"


def _side(
    gts: Mapping[str, Mapping[int, list[tuple[str, str, str]]]],
    hap: str,
    pos: int,
    depth: Mapping[int, int] | None,
) -> tuple[Records, bool, int | None]:
    """Return (records, called in another haplotype, depth) at a position.

    depth is None without a haplotype BAM.
    """
    elsewhere = any(
        pos in calls for other, calls in gts.items() if other != hap
    )
    return (
        gts[hap].get(pos),
        elsewhere,
        None if depth is None else depth.get(pos, 0),
    )


def _half_majority(votes: Sequence[str]) -> tuple[str | None, float]:
    if len(votes) < CHIMERA_MIN_HALF_SITES:
        return None, 0.0
    prox = sum(1 for v in votes if v == "proximal")
    dist = len(votes) - prox
    share = max(prox, dist) / len(votes)
    if prox == dist:
        return "tie", share
    return ("proximal" if prox > dist else "distal"), share


def vote_haplotype(
    sites: Sequence[Site],
    votes: Mapping[int, str],
    span: tuple[int, int],
    *,
    like_calls: bool = False,
) -> dict[str, object]:
    """Call one haplotype for one gene pair from per-site votes.

    Args:
        sites: Discriminating sites of the family.
        votes: {position: proximal | distal | other | missing}.
        span: Inclusive (start, end) of the haplotype boundary clipped to
            the distal gene.
        like_calls: If True, a call that is neither confident nor
            chimeric is proximal_like (50 % < share < 90 %) or
            distal_like (10 % < share < 50 %); only exactly 50 % stays
            mixed.

    Returns:
        Dict with call (proximal_type, distal_type, proximal_like,
        distal_like, chimeric, mixed or not_determined), n_sites,
        n_informative (covered sites), n_matching (covered sites matching
        the type named by the call, proximal for chimeric and mixed),
        like_pct (rounded percent of covered sites matching a like call,
        or 50 for mixed under like_calls, else empty), frac_proximal,
        left_majority, right_majority and reason.
    """
    result: dict[str, object] = {
        "call": "not_determined",
        "n_sites": 0,
        "n_informative": 0,
        "n_matching": "",
        "like_pct": "",
        "frac_proximal": "",
        "left_majority": "",
        "right_majority": "",
        "reason": "",
    }
    lo, hi = span
    if lo > hi:
        result["reason"] = "no overlap with the haplotype boundary"
        return result
    inside = [s for s in sites if lo <= s.pos <= hi]
    result["n_sites"] = len(inside)
    if len(inside) < MIN_SITES:
        result["reason"] = (
            f"only {len(inside)} discriminating site(s) in the covered span"
        )
        return result
    informative = [
        (s.pos, votes[s.pos])
        for s in inside
        if votes.get(s.pos) in ("proximal", "distal")
    ]
    result["n_informative"] = len(informative)
    if not informative:
        result["reason"] = "no informative genotype at any site"
        return result
    n_prox = sum(1 for _p, v in informative if v == "proximal")
    frac = n_prox / len(informative)
    result["frac_proximal"] = round(frac, 4)
    result["n_matching"] = n_prox
    mid = (lo + hi) // 2
    left_maj, left_margin = _half_majority(
        [v for p, v in informative if p < mid]
    )
    right_maj, right_margin = _half_majority(
        [v for p, v in informative if p >= mid]
    )
    result["left_majority"] = left_maj or ""
    result["right_majority"] = right_maj or ""
    sides = ("proximal", "distal")
    if frac >= PROXIMAL_FRACTION:
        result["call"] = "proximal_type"
    elif frac <= DISTAL_FRACTION:
        result["call"] = "distal_type"
        result["n_matching"] = len(informative) - n_prox
    elif (
        left_maj in sides
        and right_maj in sides
        and left_maj != right_maj
        and min(left_margin, right_margin) >= CHIMERA_MIN_MARGIN
    ):
        result["call"] = "chimeric"
    elif not like_calls:
        result["call"] = "mixed"
    elif 2 * n_prox == len(informative):
        result["call"] = "mixed"
        result["like_pct"] = 50
    elif frac > 0.5:
        result["call"] = "proximal_like"
        result["like_pct"] = round(100 * frac)
    else:
        result["call"] = "distal_like"
        result["n_matching"] = len(informative) - n_prox
        result["like_pct"] = round(100 * (1 - frac))
    return result


def _prefix(family: str) -> str:
    return "cyp21_blast" if family == "CYP21" else family.lower()


def cyp21_symbols(
    paraphase_json: Path,
    genes_tsv: Path | None = None,
    region_key: str = "rccx",
) -> dict[str, str]:
    """Return HP -> CYP21 gene name of each RCCX haplotype.

    Same rule as ``cyp21_paraphase`` of ``classify_rccx``: the proximal
    gene when Paraphase lists marker hits for the haplotype, else the
    distal gene.

    Args:
        paraphase_json: Paraphase JSON with the ``region_key`` entry.
        genes_tsv: Gene table; the packaged one when None.
        region_key: Key of the RCCX region.

    Returns:
        Empty dict if the region is absent.
    """
    with open(paraphase_json, encoding="utf-8") as handle:
        region = json.load(handle).get(region_key)
    if not region or not region.get("hap_variants"):
        return {}
    prox, dist = load_gene_table(genes_tsv).pair("CYP21")
    return {
        hap: prox.name if markers else dist.name
        for hap, markers in region["hap_variants"].items()
    }


def _infer_unlinked(
    groups: list[list[str]],
    pool: list[str],
    is_distal: Callable[[str], bool],
) -> tuple[list[list[str]], set[str], list[str]]:
    """Place haplotypes that Paraphase did not link to an allele.

    Every chromosome carries exactly one distal-end module (CYP21A2 with
    TNXB). An unlinked distal haplotype therefore belongs to the one
    allele that still lacks a distal end. Unlinked proximal haplotypes
    cannot be placed by this rule: either allele could carry them.

    Args:
        groups: Alleles listed by Paraphase (lists of haplotype names).
        pool: Haplotypes in no listed allele.
        is_distal: True for a haplotype that is a distal-end module.

    Returns:
        (alleles with the placed haplotypes added, placed haplotypes,
        haplotypes that stay unlinked).
    """
    alleles = [list(g) for g in groups]
    pool_d = [h for h in pool if is_distal(h)]
    pool_p = [h for h in pool if not is_distal(h)]
    placed: set[str] = set()
    if len(alleles) == 2 and len(pool_d) == 1:
        lacking = [a for a in alleles if not any(is_distal(h) for h in a)]
        if len(lacking) == 1:
            lacking[0].append(pool_d[0])
            placed.add(pool_d[0])
            pool_d = []
    elif (
        len(alleles) == 1
        and len(pool_d) == 1
        and not pool_p
        and any(is_distal(h) for h in alleles[0])
    ):
        # One linked allele with its distal end, one lone distal module.
        alleles.append([pool_d[0]])
        placed.add(pool_d[0])
        pool_d = []
    return alleles, placed, pool_d + pool_p


def _collect_votes(
    haps: Sequence[str],
    families: Sequence[str],
    sites_by_family: Mapping[str, Sequence[Site]],
    merged_vcf: Path,
    hap_bam_dir: Path | None,
    min_depth: int,
    default_chrom: str,
) -> tuple[dict[tuple[str, str], dict[int, str]], dict[str, str]]:
    """Vote every discriminating site for every haplotype and family.

    Returns:
        ({(hap, family): {position: vote}}, {hap: "bam" or "vcf"}).
    """
    all_sites = [s for f in families for s in sites_by_family.get(f, [])]
    chrom = all_sites[0].chrom if all_sites else default_chrom
    coords = [p for s in all_sites for p in (s.pos, s.proximal_pos) if p > 0]
    lo = min(coords, default=1)
    hi = max(coords, default=1)
    gts, vcf_samples = _vcf_genotypes(merged_vcf, haps, chrom, lo, hi)
    votes_by: dict[tuple[str, str], dict[int, str]] = {}
    coverage_source: dict[str, str] = {}
    for hap in haps:
        bam = (hap_bam_dir / hap / "reads.bam") if hap_bam_dir else None
        use_bam = bam is not None and bam.exists()
        coverage_source[hap] = "bam" if use_bam else "vcf"
        # A haplotype that is neither a VCF sample nor has a BAM failed or
        # was skipped upstream: it must not vote reference (distal).
        unobserved = hap not in vcf_samples and not use_bam
        for family in families:
            sites = sites_by_family.get(family, [])
            if unobserved:
                votes_by[(hap, family)] = {s.pos: "missing" for s in sites}
                continue
            d_depth = p_depth = None
            if use_bam:
                d_depth = _depths(bam, chrom, [s.pos for s in sites])
                p_depth = _depths(
                    bam,
                    chrom,
                    [s.proximal_pos for s in sites if s.proximal_pos],
                )
            votes_by[(hap, family)] = {
                site.pos: _frame_vote(
                    site,
                    _side(gts, hap, site.pos, d_depth),
                    _side(gts, hap, site.proximal_pos, p_depth),
                    min_depth,
                    use_bam,
                )
                for site in sites
            }
    return votes_by, coverage_source


def classify_tnx_haplotypes(
    paraphase_json: Path,
    merged_vcf: Path,
    sites_tsv: Path,
    *,
    genes_tsv: Path | None = None,
    hap_bam_dir: Path | None = None,
    min_depth: int = DEFAULT_MIN_DEPTH,
    region_key: str = "TNXB",
) -> dict[str, str]:
    """Return HP -> TNXB or TNXA for the haplotypes of the TNXB region.

    Paraphase phases the duplicated 3' end of TNXB, so a haplotype of this
    region is a TNXB copy (distal) or a TNXA copy (proximal). Each one is
    called with the TNX discriminating sites, exactly like the TNX module
    call of ``classify_rccx``. A haplotype with a like-call keeps its
    leaning gene; a chimeric, mixed or undetermined one is omitted.

    Args:
        paraphase_json: Paraphase JSON with a ``region_key`` entry.
        merged_vcf: Merged VCF of the run.
        sites_tsv: Output of ``build_rccx_sites``.
        genes_tsv: Gene table; the packaged one when None.
        hap_bam_dir: Folder with ``<HP>/reads.bam`` per haplotype.
        min_depth: Minimum depth for a reference vote.
        region_key: Key of the TNXB region in the Paraphase JSON.

    Returns:
        {HP: gene symbol}; empty when the region or the sites are absent.
    """
    with open(paraphase_json, encoding="utf-8") as handle:
        region = json.load(handle).get(region_key)
    table = load_gene_table(genes_tsv)
    sites_by_family = load_sites(sites_tsv)
    family = "TNX"
    if not region or not region.get("haplotype_details"):
        return {}
    if not sites_by_family.get(family):
        return {}
    details = region["haplotype_details"]
    prox, dist = table.pair(family)
    votes_by, _ = _collect_votes(
        list(details),
        [family],
        sites_by_family,
        merged_vcf,
        hap_bam_dir,
        min_depth,
        dist.chrom,
    )
    genes = {
        "proximal_type": prox.name,
        "proximal_like": prox.name,
        "distal_type": dist.name,
        "distal_like": dist.name,
    }
    out = {}
    for hap, detail in details.items():
        lo_b, hi_b = detail["boundary"]
        res = vote_haplotype(
            sites_by_family[family],
            votes_by[(hap, family)],
            (max(dist.start, lo_b), min(dist.end, hi_b)),
            like_calls=True,
        )
        gene = genes.get(str(res["call"]))
        if gene:
            out[hap] = gene
    return out


def classify_rccx(
    paraphase_json: Path,
    merged_vcf: Path,
    sites_tsv: Path,
    *,
    genes_tsv: Path | None = None,
    hap_bam_dir: Path | None = None,
    min_depth: int = DEFAULT_MIN_DEPTH,
    region_key: str = "rccx",
) -> list[dict[str, object]]:
    """Return the module and gene-conversion calls per haplotype.

    Args:
        paraphase_json: Paraphase JSON; uses ``haplotype_details``
            (boundary), ``alleles_final``, ``annotated_alleles`` and
            ``hap_variants`` of the ``region_key`` entry. Haplotype names
            are the exact HP values.
        merged_vcf: Merged VCF, one record per allele, one haploid sample
            column per HP (GT 1, 0 or .). Records at a site are matched on
            REF and ALT against the precomputed bases.
        sites_tsv: Output of ``build_rccx_sites``. Without the
            ``proximal_pos`` column (older tables) only the distal frame
            is voted.
        genes_tsv: Gene table; the packaged one when None.
        hap_bam_dir: Folder with ``<HP>/reads.bam`` per haplotype. When
            given, a site without a VCF record votes distal only if the
            haplotype BAM has depth >= ``min_depth`` there; otherwise (or
            if that BAM is missing) it votes distal only if another
            haplotype has a record there.
        min_depth: Minimum depth, see ``hap_bam_dir``.
        region_key: Key of the RCCX region in the Paraphase JSON.

    Returns:
        One dict per haplotype, ordered by allele and then by module
        position (proximal first): allele, allele_tag, n_modules,
        structure, deleted_modules, module_position, module_order
        (5' to 3' slot, a range for modules of one type whose order is
        not resolved), assignment (paraphase, inferred or unlinked), hap,
        cyp21_paraphase (CYP21A1P or CYP21A2), cyp21_marker_hits,
        coverage_source (bam or vcf) and, per family prefix (rp, c4,
        cyp21_blast, tnx), ``<prefix>_call``, ``_gene`` (resolved gene or
        empty), ``_n_sites``, ``_n_informative``, ``_frac_proximal``,
        ``_left_majority``, ``_right_majority`` and ``_reason``; at the
        end of the row, for RP, C4 and TNX only, ``<prefix>_like_gene``
        (leaning gene of a proximal_like or distal_like call, else
        empty), ``_like_pct``, ``_matching_sites`` and ``_covered_sites``
        (see ``vote_haplotype``).

    Raises:
        KeyError: If a haplotype is missing from ``hap_variants`` or
            ``haplotype_details``.
        ValueError: On a multi-allelic VCF record or an allele that lists
            an unknown haplotype.
    """
    with open(paraphase_json, encoding="utf-8") as handle:
        region = json.load(handle)[region_key]
    table = load_gene_table(genes_tsv)
    sites_by_family = load_sites(sites_tsv)
    details = region["haplotype_details"]
    haps = list(details)
    if len(set(haps)) != len(haps):
        raise ValueError("duplicate haplotype names")
    cyp_prox, cyp_dist = table.pair("CYP21")
    votes_by, coverage_source = _collect_votes(
        haps,
        table.families(),
        sites_by_family,
        merged_vcf,
        hap_bam_dir,
        min_depth,
        cyp_dist.chrom,
    )

    def cyp21_gene(hap: str) -> str:
        return cyp_prox.name if region["hap_variants"][hap] else cyp_dist.name

    listed_groups = region["alleles_final"]
    for group in listed_groups:
        for hap in group:
            if hap not in details:
                raise ValueError(f"allele lists unknown haplotype {hap!r}")
    tags = list(region.get("annotated_alleles") or [])
    assigned = {h for g in listed_groups for h in g}
    pool = [h for h in haps if h not in assigned]

    def is_distal(hap: str) -> bool:
        return cyp21_gene(hap) != cyp_prox.name

    groups, placed, leftover = _infer_unlinked(
        [list(g) for g in listed_groups], pool, is_distal
    )
    n_listed = len(groups)
    groups += [[h] for h in leftover]
    rows: list[dict[str, object]] = []
    for index, group in enumerate(groups):
        listed = index < n_listed
        # Chromosome order, 5' to 3': proximal modules (CYP21A1P side),
        # then the distal end (CYP21A2 with TNXB). Modules of one type keep
        # the Paraphase order; their mutual order is not resolved.
        ordered = sorted(group, key=is_distal)
        n_mod = len(group)
        n_distal = sum(1 for h in group if is_distal(h))
        for order, hap in enumerate(ordered, start=1):
            lo_b, hi_b = details[hap]["boundary"]
            if n_mod == 1 and listed:
                position = "sole"
            elif is_distal(hap):
                position = "distal"
            else:
                position = "proximal"
            tied = n_distal if is_distal(hap) else n_mod - n_distal
            first_slot = n_mod - n_distal + 1 if is_distal(hap) else 1
            slots = (
                str(order)
                if tied == 1 or not listed
                else f"{first_slot}-{first_slot + tied - 1}"
            )
            row: dict[str, object] = {
                "allele": index + 1,
                "allele_tag": tags[index]
                if listed and index < len(tags)
                else "",
                "n_modules": n_mod,
                "structure": (
                    STRUCTURE_NAMES.get(n_mod, f"{n_mod}-modular")
                    if listed
                    else "unassigned"
                ),
                "deleted_modules": max(0, STANDARD_MODULES - n_mod)
                if listed
                else "",
                "module_position": position,
                "module_order": slots if listed else "",
                "assignment": (
                    "inferred"
                    if hap in placed
                    else "paraphase"
                    if listed
                    else "unlinked"
                ),
                "hap": hap,
                "cyp21_paraphase": cyp21_gene(hap),
                "cyp21_marker_hits": len(region["hap_variants"][hap]),
                "coverage_source": coverage_source[hap],
            }
            extra: dict[str, object] = {}
            for family in table.families():
                prox, dist = table.pair(family)
                span = (max(dist.start, lo_b), min(dist.end, hi_b))
                res = vote_haplotype(
                    sites_by_family.get(family, []),
                    votes_by[(hap, family)],
                    span,
                    like_calls=family in LIKE_FAMILIES,
                )
                pre = _prefix(family)
                call = str(res["call"])
                row[f"{pre}_call"] = call
                row[f"{pre}_gene"] = {
                    "proximal_type": prox.name,
                    "distal_type": dist.name,
                }.get(call, "")
                for key in (
                    "n_sites",
                    "n_informative",
                    "frac_proximal",
                    "left_majority",
                    "right_majority",
                    "reason",
                ):
                    row[f"{pre}_{key}"] = res[key]
                if family in LIKE_FAMILIES:
                    extra[f"{pre}_like_gene"] = {
                        "proximal_like": prox.name,
                        "distal_like": dist.name,
                    }.get(call, "")
                    extra[f"{pre}_like_pct"] = res["like_pct"]
                    extra[f"{pre}_matching_sites"] = res["n_matching"]
                    extra[f"{pre}_covered_sites"] = res["n_informative"]
            row.update(extra)
            rows.append(row)
    return rows


def write_rccx_table(rows: Sequence[Mapping[str, object]], path: Path) -> None:
    """Write ``classify_rccx`` rows as TSV (keys of the first row)."""
    if not rows:
        raise ValueError("no rows to write")
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


# ----------------------------------------------------------------- figure


def _contrast(hex_a: str, hex_b: str) -> float:
    def lum(hex_color: str) -> float:
        chans = [int(hex_color[i : i + 2], 16) / 255 for i in (1, 3, 5)]
        lin = [
            c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4
            for c in chans
        ]
        return 0.2126 * lin[0] + 0.7152 * lin[1] + 0.0722 * lin[2]

    hi, lo = sorted((lum(hex_a), lum(hex_b)), reverse=True)
    return (hi + 0.05) / (lo + 0.05)


def _text_on(fill: str) -> str:
    return (
        "#FFFFFF"
        if _contrast(fill, "#FFFFFF") > _contrast(fill, "#000000")
        else "#000000"
    )


def _short_name(gene: str, family: str) -> str:
    """Return the label without the family prefix, e.g. A1P for CYP21A1P."""
    rest = gene[len(family) :] if gene.startswith(family) else gene
    return rest or gene


def _label_for(row: Mapping[str, object], family: str) -> tuple[str, str, str]:
    """Return (kind, gene, percent) of a family call.

    kind is "gene" with the resolved gene name, "like" with the leaning
    gene and its percent, or chimeric, mixed (percent only when given) or
    n.d. with an empty gene.
    """
    pre = _prefix(family)
    call = str(row[f"{pre}_call"])
    gene = str(row[f"{pre}_gene"])
    pct = str(row.get(f"{pre}_like_pct", ""))
    if gene:
        return "gene", gene, ""
    like = str(row.get(f"{pre}_like_gene", ""))
    if like:
        return "like", like, pct
    kind = {"chimeric": "chimeric", "mixed": "mixed"}.get(call, "n.d.")
    return kind, "", pct if kind == "mixed" else ""


def _set_font() -> dict[str, object]:
    return {**base_rc(), "hatch.linewidth": 0.5}


def _conflict(row: Mapping[str, object]) -> bool:
    """True if CYP21 BLAST resolves to the other gene than Paraphase."""
    gene = str(row.get("cyp21_blast_gene", ""))
    return bool(gene) and gene != str(row["cyp21_paraphase"])


def _module_text(row: Mapping[str, object]) -> str:
    """Return "<order>/<n> <position>", for example "1-2/3 proximal"."""
    order = str(row.get("module_order", ""))
    if not order or row["module_position"] == "sole":
        return str(row["module_position"])
    return f"{order}/{row['n_modules']} {row['module_position']}"


def _table_cells(
    row: Mapping[str, object], table: GeneTable, first: bool
) -> list[str]:
    def fam_text(family: str) -> str:
        kind, gene, pct = _label_for(row, family)
        pre = _prefix(family)
        detail = ""
        if pct:
            detail = (
                f" ({pct}%, {row.get(f'{pre}_matching_sites', '')}/"
                f"{row.get(f'{pre}_covered_sites', '')})"
            )
        if kind == "like":
            return f"{table.gene(gene).display}-like{detail}"
        if gene:
            return table.gene(gene).display
        return kind + detail

    cyp_blast = fam_text("CYP21")
    cyp = f"{table.gene(str(row['cyp21_paraphase'])).display}/{cyp_blast}" + (
        " *" if _conflict(row) else ""
    )
    return [
        str(row["allele"]) if first else "",
        str(row["allele_tag"]) if first else "",
        str(row["structure"]) if first else "",
        str(row["hap"])
        + (" \u2020" if row["assignment"] == "inferred" else ""),
        _module_text(row),
        fam_text("RP"),
        fam_text("C4"),
        cyp,
        fam_text("TNX"),
    ]


TABLE_HEAD = [
    "Allele",
    "Tag",
    "Structure",
    "Haplotype",
    "Module",
    "RP",
    "C4",
    "CYP21 (markers/BLAST)",
    "TNX",
]


def _draw_box(
    ax: plt.Axes,
    x: float,
    y: float,
    width: float,
    height: float,
    row: Mapping[str, object],
    family: str,
    table: GeneTable,
    fits: bool,
    like_fits: bool,
) -> None:
    kind, gene, pct = _label_for(row, family)
    prox, dist = table.pair(family)
    edge = {"linewidth": 0.5, "edgecolor": EDGE}

    def box(left: float = x, w: float = width, **kwargs: object) -> None:
        ax.add_patch(Rectangle((left, y), w, height, **kwargs))

    bbox = None
    if kind == "like":
        box(
            facecolor="#FFFFFF",
            linewidth=LIKE_LINEWIDTH,
            edgecolor=GENE_COLORS.get(gene, "#999999"),
        )
        name = table.gene(gene).display
        short = name if like_fits else _short_name(name, family)
        label, color = f"{short}-like\n{pct}%", "#000000"
    elif gene:
        fill = GENE_COLORS.get(gene, "#999999")
        box(facecolor=fill, **edge)
        name = table.gene(gene).display
        label = name if fits else _short_name(name, family)
        color = _text_on(fill)
    else:
        if kind == "chimeric":
            for k, g in enumerate((prox, dist)):
                box(
                    x + k * width / 2,
                    width / 2,
                    facecolor=GENE_COLORS.get(g.name, "#999999"),
                    linewidth=0,
                )
            box(fill=False, hatch="////", **edge)
        elif kind == "mixed":
            box(facecolor=UNRESOLVED_FILL[kind], hatch="xxxx", **edge)
        else:
            box(facecolor=UNRESOLVED_FILL[kind], linestyle="--", **edge)
        label, color = kind, "#000000"
        if kind == "mixed" and pct:
            label = f"mixed\n{pct}%"
        bbox = {"facecolor": "white", "edgecolor": "none", "pad": 0.6}
    ax.text(
        x + width / 2,
        y + height / 2,
        label,
        ha="center",
        va="center",
        color=color,
        bbox=bbox,
        zorder=5,
    )


def _fits(ax: plt.Axes, text: str, width_mm: float) -> bool:
    probe = ax.text(0, 0, text)
    ax.figure.canvas.draw()
    extent = probe.get_window_extent().width / ax.figure.dpi * MM
    probe.remove()
    return extent <= width_mm - 0.8


def _build_figure(
    rows: Sequence[Mapping[str, object]], table: GeneTable
) -> plt.Figure:
    families = table.families()
    alleles: dict[int, list[Mapping[str, object]]] = defaultdict(list)
    for row in rows:
        alleles[int(str(row["allele"]))].append(row)
    slots = max(STANDARD_MODULES, max(len(v) for v in alleles.values()))

    left, right, gap = 26.0, 3.0, 3.0
    row_h, box_h, head_h = 12.0, 7.0, 5.0
    cell_h = 4.6
    n_boxes = len(families)
    draw_w = WIDTH_MM - left - right
    mod_gap, box_gap = 3.0, 0.7
    box_w = (
        draw_w - (slots - 1) * mod_gap - slots * (n_boxes - 1) * box_gap
    ) / (slots * n_boxes)
    schema_h = head_h + len(alleles) * row_h
    table_h = cell_h * (len(rows) + 1)
    legend_h = 8.3 * 3.2 + 2.0
    height = 2.0 + schema_h + gap + table_h + gap + legend_h + 1.0
    if height > MAX_HEIGHT_MM:
        raise ValueError(
            f"figure would be {height:.0f} mm high "
            f"(limit {MAX_HEIGHT_MM:.0f} mm)"
        )
    # Whole pixels, so the JPEG is exactly 180 mm (1063 px at 150 dpi).
    width_in = round(WIDTH_MM / MM * FIGURE_DPI) / FIGURE_DPI
    fig = plt.figure(figsize=(width_in, height / MM))
    fig.set_facecolor("white")

    def axes_mm(x: float, top: float, w: float, h: float) -> plt.Axes:
        axis = fig.add_axes(
            (x / WIDTH_MM, 1 - (top + h) / height, w / WIDTH_MM, h / height)
        )
        axis.set_xlim(0, w)
        axis.set_ylim(0, h)
        axis.axis("off")
        return axis

    ax = axes_mm(left, 2.0, draw_w, schema_h)
    # Full gene names only if every gene of the family fits its box.
    fits = {
        f: all(_fits(ax, g.display, box_w) for g in table.pair(f))
        for f in families
    }
    # The "-like" label needs room for the thick outline as well.
    like_fits = {
        f: all(
            _fits(ax, f"{g.display}-like", box_w - 0.6) for g in table.pair(f)
        )
        for f in families
    }
    ax.text(-1.5, schema_h - 1.0, "Gene", ha="right", va="top")
    for k in range(slots):
        for j, family in enumerate(families):
            x = k * (n_boxes * (box_w + box_gap) - box_gap + mod_gap) + j * (
                box_w + box_gap
            )
            ax.text(
                x + box_w / 2, schema_h - 1.0, family, ha="center", va="top"
            )
    for i, (index, members) in enumerate(sorted(alleles.items())):
        top = schema_h - head_h - i * row_h
        y = top - row_h + 1.2
        first = members[0]
        ax.text(
            -1.5,
            y + box_h / 2 + 2.0,
            f"Allele {index}\n{first['allele_tag']}\n{first['structure']}",
            ha="right",
            va="center",
            linespacing=1.15,
        )
        offset = 0
        if len(members) == 1 and first["structure"] != "unassigned":
            width = n_boxes * box_w + (n_boxes - 1) * box_gap
            ax.add_patch(
                Rectangle(
                    (0, y),
                    width,
                    box_h,
                    fill=False,
                    linestyle="--",
                    linewidth=0.5,
                    edgecolor="#666666",
                )
            )
            ax.text(
                width / 2,
                y + box_h / 2,
                "deleted module",
                ha="center",
                va="center",
                color="#444444",
                style="italic",
            )
            offset = 1
        for k, row in enumerate(members, start=offset):
            x0 = k * (n_boxes * (box_w + box_gap) - box_gap + mod_gap)
            ax.text(
                x0 + (n_boxes * box_w + (n_boxes - 1) * box_gap) / 2,
                y + box_h + 0.8,
                str(row["hap"]),
                ha="center",
                va="bottom",
            )
            for j, family in enumerate(families):
                _draw_box(
                    ax,
                    x0 + j * (box_w + box_gap),
                    y,
                    box_w,
                    box_h,
                    row,
                    family,
                    table,
                    fits[family],
                    like_fits[family],
                )

    # Table of calls, one row per haplotype.
    cells: list[list[str]] = []
    seen: set[object] = set()
    for row in rows:
        cells.append(_table_cells(row, table, row["allele"] not in seen))
        seen.add(row["allele"])
    top_t = 2.0 + schema_h + gap
    axt = axes_mm(gap, top_t, WIDTH_MM - 2 * gap, table_h)
    lengths = [
        max(len(c[k]) for c in [TABLE_HEAD, *cells]) + 2
        for k in range(len(TABLE_HEAD))
    ]
    widths = [n / sum(lengths) for n in lengths]
    tbl = axt.table(
        cellText=cells,
        colLabels=TABLE_HEAD,
        colWidths=widths,
        loc="center",
        cellLoc="center",
        bbox=(0, 0, 1, 1),
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(FONT_PT)
    shade = {
        a: (BAND_FILL if n % 2 else "#FFFFFF") for n, a in enumerate(alleles)
    }
    for (r, c), cell in tbl.get_celld().items():
        cell.set_linewidth(0.4)
        cell.set_edgecolor("#999999")
        if r == 0:
            cell.set_facecolor(HEAD_FILL)
            cell.set_text_props(fontweight="bold")
            continue
        row = rows[r - 1]
        cell.set_facecolor(shade[int(str(row["allele"]))])
        text = cells[r - 1][c]
        if c == 7 and _conflict(row):
            cell.set_facecolor(CONFLICT_FILL)
            cell.set_hatch(CONFLICT_HATCH)
        elif text in ("chimeric", "n.d.") or text.startswith("mixed"):
            cell.set_text_props(fontweight="bold")

    # Legend and thresholds.
    top_l = top_t + table_h + gap
    axl = axes_mm(gap, top_l, WIDTH_MM - 2 * gap, legend_h)
    handles = [
        Patch(
            facecolor="white",
            edgecolor=EDGE,
            hatch="////",
            linewidth=0.5,
            label="chimeric",
        ),
        Patch(
            facecolor=UNRESOLVED_FILL["mixed"],
            edgecolor=EDGE,
            hatch="xxxx",
            linewidth=0.5,
            label="mixed",
        ),
        Patch(
            facecolor="white",
            edgecolor=EDGE,
            linestyle="--",
            linewidth=0.5,
            label="n.d. (not determined)",
        ),
        Patch(
            facecolor="white",
            edgecolor=LIKE_LEGEND_EDGE,
            linewidth=LIKE_LINEWIDTH,
            label="-like: majority type below 90 %; % = share of covered\n"
            "sites matching that type (outline = gene color)",
        ),
        Patch(
            facecolor=CONFLICT_FILL,
            edgecolor=EDGE,
            hatch=CONFLICT_HATCH,
            linewidth=0.5,
            label="* CYP21 markers and BLAST disagree",
        ),
    ]
    axl.legend(
        handles=handles,
        loc="upper left",
        ncol=3,
        frameon=False,
        handlelength=1.6,
        columnspacing=1.2,
        borderaxespad=0,
    )
    axl.text(
        0,
        0.0,
        "Calls: >= 90 % proximal-type sites = proximal gene, <= 10 % = "
        "distal gene. RP, C4, TNX: 50-90 % = proximal-like, 10-50 % =\n"
        "distal-like (% of covered sites matching the named gene), exactly "
        "50 % = mixed; CYP21 BLAST: anything else = mixed. chimeric = "
        "opposite majorities\n(>= 65 %, >= 5 sites per half); n.d. = fewer "
        "than 3 sites. CYP21: Paraphase marker hits / BLAST vote.\n"
        "Module = slot 5' to 3' / modules in the allele; a range = same-type "
        "modules, order not resolved.\n\u2020 = haplotype not linked by "
        "Paraphase, placed on the allele that lacked a distal (CYP21A2) "
        "end.\n"
        "In-house pipeline: types are voted by paralog2viewer, not taken "
        "from the Paraphase result.",
        ha="left",
        va="bottom",
        linespacing=1.3,
    )
    return fig


def plot_rccx_module(
    rows: Sequence[Mapping[str, object]],
    out_path: Path,
    *,
    genes_tsv: Path | None = None,
) -> Path:
    """Draw the RCCX module and gene-conversion figure.

    One schematic row per allele (modules left to right, proximal first;
    boxes RP, C4, CYP21, TNX colored by the resolved gene, hatched when
    chimeric or mixed, dashed when not determined; a dashed placeholder
    for a deleted module), then a table with one row per haplotype.
    JPEG, 150 dpi, quality 95, 180 mm wide, Liberation Sans 7 pt (DejaVu
    Sans if missing), Okabe-Ito colors.

    Args:
        rows: Output of ``classify_rccx``.
        out_path: JPEG path, for example ``figures/RCCX_module.jpeg``.
        genes_tsv: Gene table; the packaged one when None.

    Returns:
        ``out_path``.

    Raises:
        ValueError: If rows is empty or the figure would exceed 247 mm.
    """
    if not rows:
        raise ValueError("no haplotype rows to draw")
    table = load_gene_table(genes_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with plt.rc_context(_set_font()):
        fig = _build_figure(rows, table)
        fig.savefig(
            out_path,
            format="jpeg",
            dpi=FIGURE_DPI,
            pil_kwargs={"quality": 95, "dpi": (FIGURE_DPI, FIGURE_DPI)},
        )
        plt.close(fig)
    return out_path

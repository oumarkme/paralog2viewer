# SPDX-License-Identifier: GPL-3.0-or-later
"""Build and verify the locus-limited annotation data directory.

``build_data`` copies or subsets every annotation source to the padded
Paraphase loci and writes ``manifest.json`` with a sha256 for every file.
``verify_data`` is called by ``run`` to check the manifest and that the
calling regions of a run lie inside the subsets. The directory layout is
defined in ``datalayout``.

Coordinates: ``Region`` is ``(chrom, start0, end0)``, 0-based half-open, like
BED. Annotation tables use 1-based positions; each subset function converts
once, through ``RegionSet.overlaps``. A table row is kept when its start
position lies inside a region (padding makes this equivalent to an overlap
test for the query alleles, which lie inside the regions).

Contig names: ClinVar and HGMD omit the ``chr`` prefix, while the FASTA,
Taiwan Biobank and AlphaMissense keep it. Matching uses the name without a
leading ``chr`` (``bare_chrom``); output rows keep the source spelling.
"""

from __future__ import annotations

import bisect
import concurrent.futures as futures
import gzip
import hashlib
import json
import logging
import os
import platform
import re
import shutil
from collections.abc import (
    Collection,
    Iterable,
    Iterator,
    Mapping,
    Sequence,
)
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

import pysam

from paralog2viewer import __version__
from paralog2viewer.datalayout import DataLayout, bare_chrom
from paralog2viewer.split import load_native_regions, parse_phase_region

LOGGER = logging.getLogger("paralog2viewer")

Region = tuple[str, int, int]

MANIFEST_SCHEMA = 1
GENOME_BUILD = "GRCh38"
DEFAULT_PADDING = 50_000

# Primary-contig lengths of GRCh38 (checked against the lab FASTA index).
GRCH38_PRIMARY_LENGTHS: Mapping[str, int] = {
    "chr1": 248956422,
    "chr2": 242193529,
    "chr3": 198295559,
    "chr4": 190214555,
    "chr5": 181538259,
    "chr6": 170805979,
    "chr7": 159345973,
    "chr8": 145138636,
    "chr9": 138394717,
    "chr10": 133797422,
    "chr11": 135086622,
    "chr12": 133275309,
    "chr13": 114364328,
    "chr14": 107043718,
    "chr15": 101991189,
    "chr16": 90338345,
    "chr17": 83257441,
    "chr18": 80373285,
    "chr19": 58617616,
    "chr20": 64444167,
    "chr21": 46709983,
    "chr22": 50818468,
    "chrX": 156040895,
    "chrY": 57227415,
    "chrM": 16569,
}

# ANNOVAR hg38 files behind the table_annovar.pl protocol list of the
# Nextflow version.
# The 1000g2015aug_* protocols read the *.sites.2015_08 files.
ANNOVAR_FILTER_DATABASES: tuple[str, ...] = (
    "ALL.sites.2015_08",
    "AFR.sites.2015_08",
    "AMR.sites.2015_08",
    "EAS.sites.2015_08",
    "EUR.sites.2015_08",
    "SAS.sites.2015_08",
    "exac03",
    "esp6500siv2_all",
    "esp6500siv2_aa",
    "esp6500siv2_ea",
    "nci60",
    "avsnp147",
    "cosmic70",
    "clinvar_20240611",
    "dbnsfp47a",
    "gnomad_exome",
    "gnomad_genome",
    "dbscsnv11",
    "kaviar_20150923",
    "hrcr1",
    "regsnpintron",
    "gene4denovo201907",
    "icgc28",
    "gnomad41_exome",
    "gnomad41_genome",
)
# Region-based database: (stem, chrom, start, end) 0-based column numbers of
# the UCSC rmsk table; start is 0-based, end is exclusive.
ANNOVAR_REGION_DATABASES: tuple[tuple[str, int, int, int], ...] = (
    ("rmsk", 5, 6, 7),
)
# Gene-based databases (refGene, ensGene, knownGene, refGeneWithVer) are small
# and annotate by nearest gene, so they are copied whole with their mRNA FASTA.
ANNOVAR_GENE_FILES: tuple[str, ...] = (
    "hg38_refGene.txt",
    "hg38_refGeneMrna.fa",
    "hg38_refGeneVersion.txt",
    "hg38_refGeneWithVer.txt",
    "hg38_refGeneWithVerMrna.fa",
    "hg38_ensGene.txt",
    "hg38_ensGeneMrna.fa",
    "hg38_knownGene.txt",
    "hg38_knownGeneMrna.fa",
    "hg38_kgXref.txt",
)

_PROFILE_COPY = (
    "genemap2.txt",
    "PanelList.tsv",
    "DefaultTranscript.tsv",
    "Comments.tsv",
    "Methylated.tsv",
    "variant_catalog.json",
)
# ClinVar.tsv.gz holds GRCh37 and GRCh38 rows; only the latter are kept.
CLINVAR_ROW_FILTER: Mapping[str, tuple[str, ...]] = {"Assembly": ("GRCh38",)}
# Present while build-data runs. It lets a rerun after an interrupted build
# (no manifest yet) tell its own directory from an unrelated one.
BUILD_SENTINEL = ".build-in-progress"

# File names of other builds; never copied into the GRCh38 data directory.
_OTHER_BUILD = re.compile(r"hg19|hg18|grch37|hs37", re.IGNORECASE)


@dataclass(frozen=True)
class DataSources:
    """Source paths for ``build_data``; the five directories are required.

    Attributes:
        reference: GRCh38 FASTA (with ``.fai`` and ``.dict`` next to it).
        vep_cache: VEP cache root holding ``homo_sapiens_merged`` and
            ``plugins``.
        annovar_humandb: ANNOVAR ``humandb`` directory.
        intervar_db: InterVar ``intervardb`` directory.
        profile: MViewer ``Profiles/Default`` directory (ClinVar, Taiwan
            Biobank, HGMD, OMIM genemap2 and the small tables).
        alphamissense: AlphaMissense hg38 table; default is under
            ``vep_cache/plugins/alphamissense``.
        maxentscan: MaxEntScan data directory; default is under
            ``vep_cache/plugins/maxentscan``.
        clinvar_name: ClinVar file name inside ``profile``.
        taiwan_biobank_name: Taiwan Biobank GRCh38 file name in ``profile``.
        hgmd_name: HGMD hg38 VCF file name in ``profile``.
    """

    reference: Path
    vep_cache: Path
    annovar_humandb: Path
    intervar_db: Path
    profile: Path
    alphamissense: Path | None = None
    maxentscan: Path | None = None
    clinvar_name: str = "ClinVar.tsv.gz"
    taiwan_biobank_name: str = (
        "TaiwanBioBank_population_allelefreq_GRCh38.illumina1000"
        ".freq-only.tsv.gz"
    )
    hgmd_name: str = "HGMD_Pro_hg38.vcf.gz"

    @property
    def alphamissense_path(self) -> Path:
        """AlphaMissense table, falling back to the VEP plugins folder."""
        return self.alphamissense or (
            self.vep_cache
            / "plugins"
            / "alphamissense"
            / "AlphaMissense_hg38.tsv.gz"
        )

    @property
    def maxentscan_path(self) -> Path:
        """MaxEntScan data directory, falling back to the VEP plugins."""
        return self.maxentscan or (self.vep_cache / "plugins" / "maxentscan")


# --------------------------------------------------------------------------
# Regions
# --------------------------------------------------------------------------


def read_fai(path: Path) -> dict[str, int]:
    """Reads contig lengths from a FASTA index, in file order."""
    lengths: dict[str, int] = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                lengths[fields[0]] = int(fields[1])
    return lengths


def read_bed(path: Path) -> list[Region]:
    """Reads a BED file (0-based half-open); ignores comments and headers."""
    regions: list[Region] = []
    with open(path, encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split()
            try:
                start, end = int(fields[1]), int(fields[2])
            except (IndexError, ValueError) as err:
                raise ValueError(f"{path}:{number}: bad BED line") from err
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{number}: empty or negative span")
            regions.append((fields[0], start, end))
    return regions


def write_bed(path: Path, regions: Iterable[Region]) -> None:
    """Writes regions as a 3-column BED file."""
    with open(path, "w", encoding="utf-8") as handle:
        for chrom, start, end in regions:
            handle.write(f"{chrom}\t{start}\t{end}\n")


def regions_from_paraphase_json(path: Path) -> list[Region]:
    """Returns the unpadded loci of every gene in a Paraphase JSON.

    Uses ``phase_region`` (``38:chr:start-end``, 1-based inclusive, several
    comma-separated loci for clustered genes) and adds the native locus of
    each gene that has one in the packaged native-region table
    (``split.load_native_regions``).

    Args:
        path: Paraphase JSON.

    Returns:
        Regions, 0-based half-open, not merged.

    Raises:
        ValueError: If a ``phase_region`` is not GRCh38 or is malformed.
    """
    with open(path, encoding="utf-8") as handle:
        data = json.load(handle)
    native = load_native_regions()
    regions: list[Region] = []
    for gene, result in data.items():
        text = result.get("phase_region")
        if not text:
            if result.get("final_haplotypes"):
                raise ValueError(f"{path}: {gene} has haplotypes but no region")
            continue
        try:
            loci = parse_phase_region(text)
        except ValueError as err:
            raise ValueError(f"{path}: {err}") from err
        for chrom, lo, hi in loci:
            if not 0 < lo <= hi:
                raise ValueError(f"{path}: invalid phase_region {text!r}")
            regions.append((chrom, lo - 1, hi))
        if gene in native:
            chrom, lo, hi = native[gene]
            regions.append((chrom, lo - 1, hi))
    return regions


def merge_regions(
    regions: Iterable[Region], contig_order: Sequence[str] | None = None
) -> list[Region]:
    """Sorts regions and merges overlapping or touching ones.

    Args:
        regions: 0-based half-open regions.
        contig_order: Contig order for sorting; names not listed sort last,
            alphabetically.
    """
    rank = {name: i for i, name in enumerate(contig_order or ())}
    ordered = sorted(
        regions, key=lambda r: (rank.get(r[0], len(rank)), r[0], r[1], r[2])
    )
    merged: list[Region] = []
    for chrom, start, end in ordered:
        if merged and merged[-1][0] == chrom and start <= merged[-1][2]:
            merged[-1] = (chrom, merged[-1][1], max(merged[-1][2], end))
        else:
            merged.append((chrom, start, end))
    return merged


def build_regions(
    inputs: Path | Sequence[Path],
    padding: int,
    contig_lengths: Mapping[str, int],
) -> list[Region]:
    """Derives the padded, merged region list.

    Args:
        inputs: One BED file, or one or more Paraphase JSON files (decided
            by the ``.json`` suffix; JSON and BED are not mixed).
        padding: Bases added on each side before merging.
        contig_lengths: Reference contig lengths; padded ends are clipped to
            them and the contig order sets the output order.

    Returns:
        Merged regions, 0-based half-open.

    Raises:
        ValueError: On a contig missing from the reference, a region outside
            its contig, or mixed or empty input.
    """
    paths = (
        [Path(inputs)]
        if isinstance(inputs, (str, os.PathLike))
        else [Path(p) for p in inputs]
    )
    if not paths:
        raise ValueError("no region input given")
    is_json = [p.suffix.lower() == ".json" for p in paths]
    if any(is_json) and not all(is_json):
        raise ValueError("give either one BED or only Paraphase JSON files")
    if not all(is_json) and len(paths) != 1:
        raise ValueError("give exactly one BED file")
    raw: list[Region] = []
    for path in paths:
        raw.extend(
            regions_from_paraphase_json(path)
            if all(is_json)
            else read_bed(path)
        )
    padded: list[Region] = []
    for chrom, start, end in raw:
        if chrom not in contig_lengths:
            raise ValueError(f"region contig {chrom} is not in the reference")
        if end > contig_lengths[chrom]:
            raise ValueError(
                f"region {chrom}:{start}-{end} exceeds contig length "
                f"{contig_lengths[chrom]}"
            )
        padded.append(
            (
                chrom,
                max(0, start - padding),
                min(contig_lengths[chrom], end + padding),
            )
        )
    if not padded:
        raise ValueError("no regions found in the input")
    return merge_regions(padded, list(contig_lengths))


class RegionSet:
    """Disjoint regions with fast position and overlap lookup.

    Regions are stored 0-based half-open; queries use 1-based closed
    positions, as in VCF and ANNOVAR tables.
    """

    def __init__(self, regions: Iterable[Region]) -> None:
        """Merges the regions and indexes them by ``bare_chrom``."""
        self._starts: dict[str, list[int]] = {}
        self._ends: dict[str, list[int]] = {}
        self._names: dict[str, list[str]] = {}
        for chrom, start, end in merge_regions(regions):
            key = bare_chrom(chrom)
            self._starts.setdefault(key, []).append(start + 1)
            self._ends.setdefault(key, []).append(end)
            self._names.setdefault(key, []).append(chrom)

    def keys(self) -> set[str]:
        """Returns the contig keys (``chr`` removed) that have regions."""
        return set(self._starts)

    def intervals(self, key: str) -> list[tuple[int, int]]:
        """Returns (start0, end0) intervals of one contig key."""
        return [
            (s - 1, e)
            for s, e in zip(
                self._starts.get(key, []), self._ends.get(key, []), strict=True
            )
        ]

    def overlaps(self, key: str, low: int, high: int) -> bool:
        """Tests overlap with the 1-based closed span ``low..high``."""
        starts = self._starts.get(key)
        if not starts:
            return False
        i = bisect.bisect_right(starts, high) - 1
        return i >= 0 and self._ends[key][i] >= low


def uncovered_intervals(
    regions: Iterable[Region], intervals: Iterable[Region]
) -> list[Region]:
    """Returns the parts of ``intervals`` that no region covers.

    Contig names are compared exactly. Both inputs are 0-based half-open.
    """
    merged = merge_regions(regions)
    by_chrom: dict[str, list[tuple[int, int]]] = {}
    for chrom, start, end in merged:
        by_chrom.setdefault(chrom, []).append((start, end))
    missing: list[Region] = []
    for chrom, start, end in intervals:
        cursor = start
        for r_start, r_end in by_chrom.get(chrom, []):
            if r_end <= cursor:
                continue
            if r_start >= end:
                break
            if r_start > cursor:
                missing.append((chrom, cursor, r_start))
            cursor = max(cursor, r_end)
            if cursor >= end:
                break
        if cursor < end:
            missing.append((chrom, cursor, end))
    return missing


# --------------------------------------------------------------------------
# Subsetting
# --------------------------------------------------------------------------


PARTIAL = ".partial"  # Marks temporary files; build_data removes leftovers.


def _partial(path: Path) -> Path:
    """Returns the temporary sibling of ``path`` (same directory)."""
    return path.with_name(path.name + PARTIAL)


def _bgzip_and_index(
    plain: Path,
    dst: Path,
    *,
    preset: str | None,
    columns: tuple[int, int, int],
    line_skip: int,
) -> None:
    """Compresses ``plain`` to ``dst`` (BGZF) and writes the tabix index.

    Both files are built under temporary names next to ``dst`` and renamed
    when complete (index first), so an interrupted run never leaves a
    truncated ``dst``; existing files are overwritten.
    """
    tmp = dst.with_name(dst.name + PARTIAL + ".gz")
    pysam.tabix_compress(str(plain), str(tmp), force=True)
    if preset:
        pysam.tabix_index(str(tmp), preset=preset, force=True)
    else:
        pysam.tabix_index(
            str(tmp),
            seq_col=columns[0],
            start_col=columns[1],
            end_col=columns[2],
            line_skip=line_skip,
            force=True,
        )
    index = Path(str(dst) + ".tbi")
    os.replace(str(tmp) + ".tbi", index)
    os.replace(tmp, dst)


def _filter_columns(
    src: Path, header_line: str, keep_values: Mapping[str, Collection[str]]
) -> list[tuple[int, frozenset[str]]]:
    """Resolves ``keep_values`` column names to (index, lowered values)."""
    names = [n.strip().lstrip("#").lower() for n in header_line.split("\t")]
    resolved: list[tuple[int, frozenset[str]]] = []
    for column, values in keep_values.items():
        try:
            index = names.index(column.lower())
        except ValueError as err:
            raise ValueError(f"{src}: no column {column!r} to filter") from err
        resolved.append((index, frozenset(v.lower() for v in values)))
    return resolved


def subset_table(
    src: Path,
    dst: Path,
    regions: RegionSet,
    *,
    chrom_col: int = 0,
    pos_col: int = 1,
    header_lines: int = 1,
    keep_values: Mapping[str, Collection[str]] | None = None,
    comment_prefix: str | None = None,
    preset: str | None = None,
) -> int:
    """Subsets a (b)gzip TSV with 1-based positions and writes a tabix copy.

    Streams the whole source (these tables have no index), keeps the header
    lines and every row whose position lies in a region (and that passes
    ``keep_values``), sorts rows by contig (first-seen order) and position,
    then bgzips and indexes.

    Args:
        src: Source ``.gz`` table.
        dst: Output ``.gz``; the index is ``dst`` + ``.tbi``.
        regions: Regions to keep; contigs are matched with ``bare_chrom``.
        chrom_col: 0-based contig column.
        pos_col: 0-based position column (1-based coordinate).
        header_lines: Header lines to copy and to skip in the index.
        keep_values: Optional row filter: column name (last header line,
            compared case-insensitively) -> accepted cell values (compared
            case-insensitively). A row is kept only if every named column
            holds an accepted value.
        comment_prefix: If set, header lines are those starting with this
            prefix (for example ``#`` for VCF) and ``header_lines`` is
            ignored; incompatible with ``keep_values``.
        preset: Tabix preset for the index (``vcf``), or None to index
            ``chrom_col``/``pos_col``.

    Returns:
        Number of data rows kept.

    Raises:
        ValueError: If a position is not an integer, a ``keep_values``
            column is not in the header, or both ``keep_values`` and
            ``comment_prefix`` are given.
    """
    if comment_prefix and keep_values:
        raise ValueError("keep_values needs a fixed header_lines count")
    rank: dict[str, int] = {}
    wanted: list[tuple[int, frozenset[str]]] = []
    kept: list[tuple[int, int, str]] = []
    header: list[str] = []
    with gzip.open(src, "rt", encoding="utf-8", newline="") as handle:
        for number, line in enumerate(handle):
            if (
                line.startswith(comment_prefix)
                if comment_prefix
                else number < header_lines
            ):
                header.append(line)
                if keep_values and number == header_lines - 1:
                    wanted = _filter_columns(src, line, keep_values)
                continue
            if wanted:
                cells = line.rstrip("\r\n").split("\t")
                if any(
                    idx >= len(cells) or cells[idx].strip().lower() not in ok
                    for idx, ok in wanted
                ):
                    continue
            fields = line.split("\t", pos_col + 1)
            chrom = fields[chrom_col]
            rank.setdefault(chrom, len(rank))
            try:
                pos = int(fields[pos_col])
            except (IndexError, ValueError) as err:
                raise ValueError(f"{src}:{number + 1}: bad position") from err
            if regions.overlaps(bare_chrom(chrom), pos, pos):
                kept.append((rank[chrom], pos, line))
    kept.sort(key=lambda item: item[:2])
    plain = _partial(dst.with_suffix(dst.suffix + ".plain"))
    with open(plain, "w", encoding="utf-8", newline="") as out:
        out.writelines(header)
        out.writelines(item[2] for item in kept)
    _bgzip_and_index(
        plain,
        dst,
        preset=preset,
        columns=(chrom_col, pos_col, pos_col),
        line_skip=0 if comment_prefix else header_lines,
    )
    plain.unlink()
    return len(kept)


def subset_tabix(
    src: Path,
    dst: Path,
    regions: RegionSet,
    *,
    preset: str | None = None,
    columns: tuple[int, int, int] = (0, 1, 1),
) -> int:
    """Extracts regions from a bgzip table and re-indexes the result.

    Keeps the full header and the rows whose position (second column) lies
    in a region. Contigs follow the source order, so the output is sorted.
    With a ``.tbi``/``.csi`` index beside ``src`` the regions are fetched;
    without one the file is streamed once (header = ``#`` lines), as for
    the other profile tables. Both paths give the same rows.

    Args:
        src: Source bgzip file; a tabix index is optional (fast path).
        dst: Output bgzip file; the index is ``dst`` + ``.tbi``.
        regions: Regions to extract.
        preset: Tabix preset for the new index (``vcf``), or None to use
            ``columns``.
        columns: 0-based (contig, start, end) columns when no preset.

    Returns:
        Number of data rows kept.
    """
    if not any(Path(str(src) + ext).exists() for ext in (".tbi", ".csi")):
        return subset_table(
            src,
            dst,
            regions,
            chrom_col=columns[0],
            pos_col=columns[1],
            comment_prefix="#",
            preset=preset,
        )
    count = 0
    plain = _partial(dst.with_suffix(dst.suffix + ".plain"))
    with (
        pysam.TabixFile(str(src)) as tbx,
        open(plain, "w", encoding="utf-8", newline="") as out,
    ):
        header = list(tbx.header)
        for line in header:
            out.write(line + "\n")
        for contig in tbx.contigs:
            key = bare_chrom(contig)
            for start0, end0 in regions.intervals(key):
                for row in tbx.fetch(contig, start0, end0):
                    pos = int(row.split("\t", 2)[1])
                    # fetch() returns overlapping rows; keep the start rule.
                    if start0 < pos <= end0:
                        out.write(row + "\n")
                        count += 1
    _bgzip_and_index(plain, dst, preset=preset, columns=columns, line_skip=0)
    plain.unlink()
    return count


class _IdxBuilder:
    """Collects an ANNOVAR ``.idx`` index while a table is written.

    Mirrors ``index_annovar.pl``: each (contig, bin) keeps the byte offset of
    its first line and the end of its last line. ANNOVAR seeks to the first
    offset and tests every line up to the second by position, so the rows of
    a bin need not be contiguous (several public tables are not sorted) and
    the range may then also cover rows of other bins.
    """

    def __init__(self, binsize: int) -> None:
        """Creates a builder for bins of ``binsize`` bases."""
        self.binsize = binsize
        self._entries: dict[tuple[str, int], list[int]] = {}
        self._current: tuple[str, int] | None = None
        self.reentries = 0

    def add(self, chrom: str, start: int, offset: int, length: int) -> None:
        """Registers one line at byte ``offset`` of the table."""
        key = (chrom, start - start % self.binsize)
        entry = self._entries.get(key)
        if entry is None:
            self._entries[key] = [offset, offset + length]
        else:
            if key != self._current:
                self.reentries += 1
            entry[1] = offset + length
        self._current = key

    def write(self, path: Path, table_size: int) -> None:
        """Writes the index: ``#BIN``, bin size, table size, then bins."""
        with open(path, "w", encoding="utf-8", newline="") as out:
            out.write(f"#BIN\t{self.binsize}\t{table_size}\n")
            for (chrom, first), (low, high) in self._entries.items():
                out.write(f"{chrom}\t{first}\t{low}\t{high}\n")

    def warn_if_unsorted(self, table: Path) -> None:
        """Logs a warning naming ``table`` when a bin was seen again."""
        if self.reentries:
            LOGGER.warning(
                "%s is not position-sorted (%d bin re-entries); row order "
                "is kept and its .idx bins span interleaved rows, as "
                "ANNOVAR's index_annovar.pl writes them",
                table.name,
                self.reentries,
            )


def write_annovar_idx(table: Path, binsize: int) -> Path:
    """Writes ``table`` + ``.idx`` for an ANNOVAR text table.

    The bin of a row is ``start - start % binsize`` (column 2). Offsets are
    bytes of the first line and of the end of the last line of the bin. An
    unsorted table is indexed too (with a warning), as ANNOVAR does.
    """
    builder = _IdxBuilder(binsize)
    offset = 0
    with open(table, "rb") as handle:
        for line in handle:
            if line.strip() and not line.startswith(b"#"):
                fields = line.split(b"\t", 2)
                builder.add(
                    fields[0].decode(), int(fields[1]), offset, len(line)
                )
            offset += len(line)
    builder.warn_if_unsorted(table)
    idx = table.with_name(table.name + ".idx")
    builder.write(idx, offset)
    return idx


def _read_idx_ranges(
    idx: Path, table_size: int, regions: RegionSet
) -> tuple[int, list[tuple[int, int]]] | None:
    """Returns (bin size, merged byte ranges) of bins touching the regions.

    Returns None when the index is missing, malformed or out of date (table
    size differs), as ANNOVAR itself ignores such an index.
    """
    if not idx.exists():
        return None
    with open(idx, encoding="utf-8") as handle:
        match = re.fullmatch(r"#BIN\t(\d+)\t(\d+)", handle.readline().strip())
        if not match or int(match.group(2)) != table_size:
            LOGGER.warning("%s is malformed or out of date; full scan", idx)
            return None
        binsize = int(match.group(1))
        ranges: list[tuple[int, int]] = []
        for line in handle:
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 4:
                continue
            first = int(fields[1])
            # Bin covers 1-based positions first .. first + binsize - 1.
            if regions.overlaps(
                bare_chrom(fields[0]), first, first + binsize - 1
            ):
                ranges.append((int(fields[2]), int(fields[3])))
    ranges.sort()
    merged: list[tuple[int, int]] = []
    for low, high in ranges:
        if merged and low <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], high))
        else:
            merged.append((low, high))
    return binsize, merged


def _lines_in_ranges(
    handle, ranges: list[tuple[int, int]] | None
) -> Iterator[bytes]:
    """Yields lines of an open binary file, restricted to byte ranges."""
    if ranges is None:
        handle.seek(0)
        yield from handle
        return
    for low, high in ranges:
        handle.seek(low)
        position = low
        while position < high:
            line = handle.readline()
            if not line:
                break
            position += len(line)
            yield line


def subset_annovar_table(
    src: Path,
    dst: Path,
    regions: RegionSet,
    *,
    chrom_col: int = 0,
    start_col: int = 1,
    zero_based_start: bool = False,
    end_col: int | None = None,
) -> int:
    """Subsets an ANNOVAR text table and rebuilds its ``.idx``.

    Leading ``#`` lines (header) are kept. When ``src`` has an up-to-date
    ``.idx``, only the byte ranges of bins that touch the regions are read
    (the full tables are up to 100 GB); otherwise the file is scanned. The
    new ``.idx`` (same bin size) is written only if ``src`` had a usable one.
    Rows keep the source order, so an unsorted source (a warning is logged)
    gives the same rows as a sorted one would.

    Args:
        src: Source table.
        dst: Output table.
        regions: Regions to keep.
        chrom_col: 0-based contig column.
        start_col: 0-based start column.
        zero_based_start: True when the start column is 0-based (UCSC);
            the row is kept when ``[start, end)`` overlaps a region.
        end_col: End column for the overlap test; None tests the start only.

    Returns:
        Number of data rows kept.
    """
    size = src.stat().st_size
    idx = _read_idx_ranges(src.with_name(src.name + ".idx"), size, regions)
    builder = _IdxBuilder(idx[0]) if idx else None
    ranges = idx[1] if idx else None
    count = 0
    offset = 0
    tmp = _partial(dst)
    with open(src, "rb") as handle, open(tmp, "wb") as out:
        while True:
            here = handle.tell()
            line = handle.readline()
            if not line.startswith(b"#"):
                handle.seek(here)
                break
            out.write(line)
            offset += len(line)
        for line in _lines_in_ranges(handle, ranges):
            if not line.strip() or line.startswith(b"#"):
                continue
            fields = line.rstrip(b"\r\n").split(b"\t")
            chrom = fields[chrom_col].decode()
            start = int(fields[start_col]) + (1 if zero_based_start else 0)
            end = int(fields[end_col]) if end_col is not None else start
            if not regions.overlaps(bare_chrom(chrom), start, max(start, end)):
                continue
            out.write(line)
            if builder:
                builder.add(
                    chrom,
                    start - (1 if zero_based_start else 0),
                    offset,
                    len(line),
                )
            offset += len(line)
            count += 1
    if builder:
        builder.warn_if_unsorted(src)
        idx_dst = dst.with_name(dst.name + ".idx")
        builder.write(_partial(idx_dst), offset)
        os.replace(_partial(idx_dst), idx_dst)
    os.replace(tmp, dst)
    return count


def _read_synonyms(path: Path) -> dict[str, set[str]]:
    """Reads the VEP ``chr_synonyms.txt`` pairs as an undirected map."""
    synonyms: dict[str, set[str]] = {}
    if path.exists():
        with open(path, encoding="utf-8") as handle:
            for line in handle:
                fields = line.split()
                if len(fields) >= 2:
                    synonyms.setdefault(fields[0], set()).add(fields[1])
                    synonyms.setdefault(fields[1], set()).add(fields[0])
    return synonyms


def select_vep_contigs(
    cache_dir: Path, contigs: Iterable[str]
) -> tuple[list[str], list[str]]:
    """Maps region contigs to the contig directories of a VEP cache.

    Args:
        cache_dir: ``homo_sapiens_merged/<release>_GRCh38`` directory.
        contigs: Region contig names (``chr1``, ...).

    Returns:
        (directories present in the cache, contigs with no directory).
    """
    synonyms = _read_synonyms(cache_dir / "chr_synonyms.txt")
    found: set[str] = set()
    missing: list[str] = []
    for contig in sorted(set(contigs)):
        names = {contig, bare_chrom(contig)}
        for name in list(names):
            names |= synonyms.get(name, set())
        hits = {n for n in names if (cache_dir / n).is_dir()}
        if hits:
            found |= hits
        else:
            missing.append(contig)
    return sorted(found), missing


# --------------------------------------------------------------------------
# Copy, hash, manifest
# --------------------------------------------------------------------------


def _copy_file(src: Path, dst: Path) -> None:
    """Copies file contents (no metadata) and creates parent directories."""
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = _partial(dst)
    shutil.copyfile(src, tmp)
    os.replace(tmp, dst)


def _copy_tree(src: Path, dst: Path) -> list[tuple[Path, Path]]:
    """Lists (source, target) file pairs of a directory tree.

    Files whose name marks another genome build (hg19, GRCh37; for example
    InterVar's ``*.hg19`` tables) are left out.
    """
    return [
        (path, dst / path.relative_to(src))
        for path in sorted(src.rglob("*"))
        if path.is_file() and not _OTHER_BUILD.search(path.name)
    ]


def _check_source_names(sources: DataSources) -> None:
    """Rejects profile table names that are not GRCh38.

    Raises:
        ValueError: If a name marks another build, or the Taiwan Biobank
            name does not say GRCh38/hg38 (its hg19 file has no build tag).
    """
    for name in (
        sources.clinvar_name,
        sources.taiwan_biobank_name,
        sources.hgmd_name,
    ):
        if _OTHER_BUILD.search(name):
            raise ValueError(f"{name}: not a GRCh38 table")
    if not re.search(r"grch38|hg38", sources.taiwan_biobank_name, re.I):
        raise ValueError(f"{sources.taiwan_biobank_name}: not marked GRCh38")


def _run_copies(pairs: Sequence[tuple[Path, Path]], threads: int) -> None:
    """Copies file pairs with a thread pool; re-raises the first error."""
    with futures.ThreadPoolExecutor(max_workers=threads) as pool:
        for job in [pool.submit(_copy_file, s, d) for s, d in pairs]:
            job.result()


def sha256_file(path: Path) -> str:
    """Returns the sha256 hex digest of a file."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _file_record(root: Path, path: Path) -> tuple[str, dict[str, object]]:
    """Returns (relative path, {size, sha256}) for one file."""
    return path.relative_to(root).as_posix(), {
        "size": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def _list_files(root: Path) -> list[Path]:
    """Lists all files below ``root`` except the manifest, sorted."""
    return sorted(
        p
        for p in root.rglob("*")
        if p.is_file() and p.relative_to(root).as_posix() != "manifest.json"
    )


def _check_reference(
    fasta: Path, expected: Mapping[str, int]
) -> dict[str, int]:
    """Asserts the FASTA index is GRCh38 and matches the dictionary.

    Returns:
        Contig lengths from the FASTA index.

    Raises:
        ValueError: If a primary contig is missing or has another length,
            or ``.dict`` and ``.fai`` disagree.
    """
    fai = Path(str(fasta) + ".fai")
    if not fai.exists():
        pysam.faidx(str(fasta))
    lengths = read_fai(fai)
    for name, length in expected.items():
        if lengths.get(name) != length:
            raise ValueError(
                f"{fasta} is not {GENOME_BUILD}: {name} length "
                f"{lengths.get(name)} != {length}"
            )
    dictionary = fasta.with_suffix(".dict")
    if not dictionary.exists():
        raise FileNotFoundError(f"missing sequence dictionary {dictionary}")
    with open(dictionary, encoding="utf-8") as handle:
        dict_lengths = {}
        for line in handle:
            if line.startswith("@SQ"):
                tags = dict(f.split(":", 1) for f in line.split()[1:])
                dict_lengths[tags["SN"]] = int(tags["LN"])
    if dict_lengths != lengths:
        raise ValueError(f"{dictionary} and {fai} list different contigs")
    return lengths


def _tool_versions() -> dict[str, str]:
    """Returns versions of the tools that make the data."""
    return {
        "paralog2viewer": __version__,
        "python": platform.python_version(),
        "pysam": pysam.__version__,
        "htslib": getattr(pysam.version, "__htslib_version__", "unknown"),
    }


def _vep_info(cache_dir: Path) -> dict[str, str]:
    """Reads ``info.txt`` of the VEP cache (release and source versions)."""
    info: dict[str, str] = {}
    path = cache_dir / "info.txt"
    if path.exists():
        with open(path, encoding="utf-8") as handle:
            for line in handle:
                if "\t" in line and not line.startswith("#"):
                    key, value = line.rstrip("\n").split("\t", 1)
                    if key in {"assembly", "species"} or key.startswith(
                        "source_"
                    ):
                        info[key] = value
    return info


def _subset_job(
    kind: str, src: str, dst: str, regions: RegionSet, options: dict
) -> tuple[str, int]:
    """Runs one subset in a worker process; returns (name, rows kept)."""
    source, target = Path(src), Path(dst)
    if kind == "table":
        rows = subset_table(source, target, regions, **options)
    elif kind == "tabix":
        rows = subset_tabix(source, target, regions, **options)
    else:
        rows = subset_annovar_table(source, target, regions, **options)
    return target.name, rows


def build_data(
    out_dir: Path,
    sources: DataSources,
    regions: Path | Sequence[Path],
    padding: int = DEFAULT_PADDING,
    threads: int = 4,
    *,
    expected_lengths: Mapping[str, int] = GRCH38_PRIMARY_LENGTHS,
) -> dict:
    """Builds the annotation data directory and its manifest.

    Steps: derive the padded regions; copy reference, InterVar, MaxEntScan,
    small profile files, VEP plugin code and the VEP cache contigs the
    regions touch; subset ClinVar (GRCh38 rows only), Taiwan Biobank, HGMD,
    AlphaMissense and the ANNOVAR databases of the protocol list
    (rebuilding their indexes); hash every file and write ``manifest.json``.
    RCCX discriminating sites are not made here.

    Args:
        out_dir: Data directory. It may hold an earlier or interrupted
            build: temporary files and the old manifest are removed first,
            and every output is rewritten under a temporary name and
            renamed when complete.
        sources: Source paths.
        regions: One BED file or one or more Paraphase JSON files.
        padding: Bases added on both sides of each locus.
        threads: Worker count (copy, subset and hash); keep at most 4.
        expected_lengths: Contig lengths the reference must have; the
            default is GRCh38 primary contigs.

    Returns:
        The manifest as a dict.

    Raises:
        ValueError: If the reference is not GRCh38, a table name is not
            GRCh38, or regions are invalid.
    """
    out_dir = Path(out_dir)
    _check_source_names(sources)
    layout = DataLayout(out_dir)
    lengths = _check_reference(sources.reference, expected_lengths)
    region_list = build_regions(regions, padding, lengths)
    region_set = RegionSet(region_list)
    if (
        out_dir.is_dir()
        and any(out_dir.iterdir())
        and not layout.manifest.exists()
        and not (out_dir / BUILD_SENTINEL).exists()
    ):
        raise ValueError(
            f"{out_dir} is not empty and has no {layout.manifest.name}; "
            "refusing to build into a directory this tool did not make"
        )
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / BUILD_SENTINEL).write_text("build in progress\n")
    for stale in out_dir.rglob("*"):
        if stale.is_file() and PARTIAL in stale.name:
            stale.unlink()
    layout.manifest.unlink(missing_ok=True)
    write_bed(_partial(layout.regions_bed), region_list)
    os.replace(_partial(layout.regions_bed), layout.regions_bed)

    copies: list[tuple[Path, Path]] = [
        (sources.reference, layout.reference_fasta),
        (
            Path(str(sources.reference) + ".fai"),
            Path(str(layout.reference_fasta) + ".fai"),
        ),
        (
            sources.reference.with_suffix(".dict"),
            layout.reference_fasta.with_suffix(".dict"),
        ),
    ]
    copies += _copy_tree(sources.intervar_db, layout.intervardb)
    copies += _copy_tree(sources.maxentscan_path, layout.maxentscan)
    for name in _PROFILE_COPY:
        copies.append((sources.profile / name, layout.profile / name))
    plugin_src = sources.vep_cache / "plugins"
    for plugin in sorted(plugin_src.glob("*.pm")):
        copies.append((plugin, layout.vep_plugins / plugin.name))
    for name in ANNOVAR_GENE_FILES:
        copies.append(
            (sources.annovar_humandb / name, layout.annovar_humandb / name)
        )

    cache_src = next(
        iter(
            sorted((sources.vep_cache / "homo_sapiens_merged").glob("*_GRCh38"))
        ),
        None,
    )
    if cache_src is None:
        raise FileNotFoundError("no homo_sapiens_merged/*_GRCh38 in VEP cache")
    contig_dirs, missing_contigs = select_vep_contigs(
        cache_src, [r[0] for r in region_list]
    )
    if missing_contigs:
        LOGGER.warning("no VEP cache directory for %s", missing_contigs)
    cache_dst = layout.vep_cache / "homo_sapiens_merged" / cache_src.name
    for name in ("info.txt", "chr_synonyms.txt"):
        if (cache_src / name).exists():
            copies.append((cache_src / name, cache_dst / name))
    for name in contig_dirs:
        copies += _copy_tree(cache_src / name, cache_dst / name)

    jobs: list[tuple[str, str, str, RegionSet, dict]] = []
    profile, humandb = sources.profile, sources.annovar_humandb
    jobs.append(
        (
            "table",
            str(profile / sources.clinvar_name),
            str(layout.clinvar),
            region_set,
            {"keep_values": CLINVAR_ROW_FILTER},
        )
    )
    jobs.append(
        (
            "table",
            str(profile / sources.taiwan_biobank_name),
            str(layout.taiwan_biobank),
            region_set,
            {},
        )
    )
    jobs.append(
        (
            "tabix",
            str(profile / sources.hgmd_name),
            str(layout.hgmd),
            region_set,
            {"preset": "vcf"},
        )
    )
    jobs.append(
        (
            "tabix",
            str(sources.alphamissense_path),
            str(layout.alphamissense),
            region_set,
            {"columns": (0, 1, 1)},
        )
    )
    for stem in ANNOVAR_FILTER_DATABASES:
        name = f"hg38_{stem}.txt"
        jobs.append(
            (
                "annovar",
                str(humandb / name),
                str(layout.annovar_humandb / name),
                region_set,
                {},
            )
        )
    for stem, chrom_col, start_col, end_col in ANNOVAR_REGION_DATABASES:
        name = f"hg38_{stem}.txt"
        jobs.append(
            (
                "annovar",
                str(humandb / name),
                str(layout.annovar_humandb / name),
                region_set,
                {
                    "chrom_col": chrom_col,
                    "start_col": start_col,
                    "end_col": end_col,
                    "zero_based_start": True,
                },
            )
        )

    for _, _, dst, _, _ in jobs:
        Path(dst).parent.mkdir(parents=True, exist_ok=True)
    _run_copies(copies, threads)
    row_counts: dict[str, int] = {}
    with futures.ProcessPoolExecutor(max_workers=threads) as pool:
        pending = [pool.submit(_subset_job, *job) for job in jobs]
        for job in pending:
            name, rows = job.result()
            row_counts[name] = rows
            LOGGER.info("subset %s: %d rows", name, rows)

    (out_dir / BUILD_SENTINEL).unlink(missing_ok=True)
    records: dict[str, dict[str, object]] = {}
    with futures.ThreadPoolExecutor(max_workers=threads) as pool:
        for rel, record in pool.map(
            lambda p: _file_record(out_dir, p), _list_files(out_dir)
        ):
            records[rel] = record
    digest = hashlib.sha256()
    for rel, record in sorted(records.items()):
        digest.update(f"{rel}\t{record['sha256']}\n".encode())
    manifest = {
        "schema_version": MANIFEST_SCHEMA,
        "data_version": digest.hexdigest()[:16],
        "genome_build": GENOME_BUILD,
        "created_utc": datetime.now(UTC).isoformat(timespec="seconds"),
        "padding": padding,
        "regions": [list(r) for r in region_list],
        "contig_lengths_checked": dict(expected_lengths),
        "tool_versions": _tool_versions(),
        "vep_cache": _vep_info(cache_src),
        "vep_contigs_copied": contig_dirs,
        "vep_contigs_missing": missing_contigs,
        "annovar_databases": [f"hg38_{s}.txt" for s in ANNOVAR_FILTER_DATABASES]
        + [f"hg38_{s[0]}.txt" for s in ANNOVAR_REGION_DATABASES]
        + list(ANNOVAR_GENE_FILES),
        "subset_rows": dict(sorted(row_counts.items())),
        "sources": {
            "reference": str(sources.reference),
            "vep_cache": str(sources.vep_cache),
            "annovar_humandb": str(sources.annovar_humandb),
            "intervar_db": str(sources.intervar_db),
            "profile": str(sources.profile),
            "alphamissense": str(sources.alphamissense_path),
            "maxentscan": str(sources.maxentscan_path),
            "regions_input": [
                str(p)
                for p in (
                    [regions]
                    if isinstance(regions, (str, os.PathLike))
                    else regions
                )
            ],
        },
        "files": records,
    }
    tmp_manifest = _partial(layout.manifest)
    with open(tmp_manifest, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(tmp_manifest, layout.manifest)
    return manifest


def verify_data(
    data_dir: Path,
    calling_bed: Path | Sequence[Path] | None = None,
    *,
    checksums: bool = True,
) -> list[str]:
    """Checks a data directory and optionally calling regions against it.

    Checks: manifest present and readable; genome build is GRCh38; the
    FASTA index matches the contig lengths recorded at build time; every
    manifest file exists with the recorded size and, if ``checksums``, the
    recorded sha256; no unlisted file; every interval of ``calling_bed``
    lies inside ``regions.bed``.

    Args:
        data_dir: Data directory made by ``build_data``.
        calling_bed: Calling BED file(s), 0-based half-open.
        checksums: Hash all files (reads the whole directory, slow for the
            VEP cache); when False only existence and size are checked.

    Returns:
        Problem descriptions; empty when everything is fine. Intervals that
        regions.bed does not cover are reported as
        ``outside regions: chrom:start-end`` (0-based half-open).
    """
    layout = DataLayout(Path(data_dir))
    if not layout.manifest.exists():
        return [f"manifest missing: {layout.manifest}"]
    try:
        manifest = layout.read_manifest()
    except (OSError, json.JSONDecodeError) as err:
        return [f"manifest unreadable: {err}"]
    problems: list[str] = []
    if manifest.get("schema_version") != MANIFEST_SCHEMA:
        problems.append(f"manifest schema {manifest.get('schema_version')}")
    if manifest.get("genome_build") != GENOME_BUILD:
        problems.append(f"genome build is {manifest.get('genome_build')}")
    fai = Path(str(layout.reference_fasta) + ".fai")
    if fai.exists():
        lengths = read_fai(fai)
        for name, length in manifest.get("contig_lengths_checked", {}).items():
            if lengths.get(name) != length:
                problems.append(f"reference {name} length differs")
    else:
        problems.append(f"missing {fai}")
    files = manifest.get("files", {})
    for rel, record in sorted(files.items()):
        path = layout.root / rel
        if not path.is_file():
            problems.append(f"missing file: {rel}")
        elif path.stat().st_size != record["size"]:
            problems.append(f"size differs: {rel}")
    if checksums:
        present = [
            rel
            for rel in files
            if (layout.root / rel).is_file()
            and (layout.root / rel).stat().st_size == files[rel]["size"]
        ]
        with futures.ThreadPoolExecutor(max_workers=4) as pool:
            digests = pool.map(
                lambda rel: sha256_file(layout.root / rel), present
            )
            for rel, digest in zip(present, digests, strict=True):
                if digest != files[rel]["sha256"]:
                    problems.append(f"checksum differs: {rel}")
    listed = set(files)
    for path in _list_files(layout.root):
        if path.relative_to(layout.root).as_posix() not in listed:
            problems.append(
                f"unlisted file: {path.relative_to(layout.root).as_posix()}"
            )
    if calling_bed is not None and layout.regions_bed.exists():
        beds = (
            [Path(calling_bed)]
            if isinstance(calling_bed, (str, os.PathLike))
            else [Path(p) for p in calling_bed]
        )
        calling = [iv for bed in beds for iv in read_bed(bed)]
        for chrom, start, end in uncovered_intervals(
            read_bed(layout.regions_bed), calling
        ):
            problems.append(f"outside regions: {chrom}:{start}-{end}")
    elif calling_bed is not None:
        problems.append(f"missing {layout.regions_bed}")
    return problems

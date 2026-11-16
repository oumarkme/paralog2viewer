# SPDX-License-Identifier: GPL-3.0-or-later
"""Layout of the data directory made by ``build-data`` and read by ``run``.

This module is the single contract between the data builder and the modules
that read annotation data. Only paths are defined here.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

MANIFEST_NAME = "manifest.json"


@dataclass(frozen=True)
class DataLayout:
    """Paths inside one data directory.

    Attributes:
        root: The data directory.
    """

    root: Path

    @property
    def manifest(self) -> Path:
        """JSON file with data version, genome build and file checksums."""
        return self.root / MANIFEST_NAME

    @property
    def reference_fasta(self) -> Path:
        """GRCh38 FASTA; ``.fai`` and ``.dict`` sit next to it."""
        return self.root / "reference" / "genome.fa"

    @property
    def regions_bed(self) -> Path:
        """Padded regions the subsets cover (0-based half-open BED)."""
        return self.root / "regions.bed"

    @property
    def vep_cache(self) -> Path:
        """VEP cache directory (holds ``homo_sapiens_merged``)."""
        return self.root / "vep_cache"

    @property
    def vep_plugins(self) -> Path:
        """VEP plugin code directory."""
        return self.root / "vep_cache" / "plugins"

    @property
    def maxentscan(self) -> Path:
        """MaxEntScan plugin data directory."""
        return self.root / "maxentscan"

    @property
    def alphamissense(self) -> Path:
        """AlphaMissense table (bgzip, tabix-indexed)."""
        return self.root / "alphamissense" / "AlphaMissense_hg38.tsv.gz"

    @property
    def annovar_humandb(self) -> Path:
        """ANNOVAR database directory with only the databases in use."""
        return self.root / "annovar_humandb"

    @property
    def intervardb(self) -> Path:
        """InterVar database directory."""
        return self.root / "intervardb"

    @property
    def profile(self) -> Path:
        """Directory with the former MViewer profile tables."""
        return self.root / "profile"

    @property
    def clinvar(self) -> Path:
        """ClinVar table (bgzip, tabix-indexed)."""
        return self.profile / "ClinVar.tsv.gz"

    @property
    def taiwan_biobank(self) -> Path:
        """Taiwan Biobank GRCh38 allele-frequency table (bgzip, indexed)."""
        return self.profile / "TaiwanBioBank_GRCh38.freq-only.tsv.gz"

    @property
    def hgmd(self) -> Path:
        """HGMD Professional hg38 VCF (bgzip, tabix-indexed)."""
        return self.profile / "HGMD_Pro_hg38.vcf.gz"

    @property
    def genemap2(self) -> Path:
        """OMIM genemap2 table."""
        return self.profile / "genemap2.txt"

    @property
    def panel_list(self) -> Path:
        """Gene panel list."""
        return self.profile / "PanelList.tsv"

    @property
    def default_transcript(self) -> Path:
        """Default transcript per gene."""
        return self.profile / "DefaultTranscript.tsv"

    @property
    def comments(self) -> Path:
        """Curated comments table."""
        return self.profile / "Comments.tsv"

    @property
    def methylated(self) -> Path:
        """Methylated loci table."""
        return self.profile / "Methylated.tsv"

    @property
    def variant_catalog(self) -> Path:
        """Variant catalog JSON."""
        return self.profile / "variant_catalog.json"

    @property
    def rccx_sites(self) -> Path:
        """Precomputed RCCX discriminating sites (TSV)."""
        return self.root / "rccx" / "discriminating_sites.tsv"

    def read_manifest(self) -> dict:
        """Return the parsed manifest.

        Raises:
            FileNotFoundError: If the manifest does not exist.
        """
        with open(self.manifest, encoding="utf-8") as handle:
            return json.load(handle)


def allele_key(chrom: str, pos: int, ref: str, alt: str) -> str:
    """Return the allele key used to join annotation to calls.

    The key uses normalized VCF coordinates (1-based POS), for example
    ``chr6:32038400:G:A``.
    """
    return f"{chrom}:{pos}:{ref}:{alt}"


def bare_chrom(name: str) -> str:
    """Return a contig name without a leading ``chr`` (for matching only).

    A name that is just ``chr`` is returned unchanged.
    """
    return name[3:] if name.startswith("chr") and len(name) > 3 else name

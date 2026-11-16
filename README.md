# paralog2viewer

Per-haplotype variant calls, annotation, and read-pileup figures for paralogous genes, starting from an Oxford Nanopore (ONT) Paraphase result.

Version 0.9.0 (pre-release). License GPL-3.0-or-later. Research use only, not a clinical diagnostic.

For each sample you get one annotated variant call format (VCF) file, one 365-column MViewer-format table, and one figure per gene region. A figure shows every read of a region as one row, grouped by haplotype (`HP` tag). Colored dots mark alleles that some haplotypes carry and others do not. White diamonds mark ClinVar or InterVar Pathogenic and Likely pathogenic alleles. Black squares mark pseudogene-type alleles from a curated table. For the RCCX region (RP-C4-CYP21-TNX modules), a second figure shows the module structure and the gene-conversion call per haplotype.

## What it does

Input: an existing Paraphase BAM, its `.bai` index, and its JSON. Paraphase is not rerun.

```
Paraphase BAM + .bai + JSON
  [0] check    BAM dictionary vs GRCh38 FASTA, data manifest, unique HP values
  [1] split    confident primary reads (RN, HP, YC tags) per gene region and haplotype
  [2] realign  minimap2 (map-ont) from raw sequence -> calling BED per haplotype
  [3] call     Clair3 (SNV/indel) and Sawfish (SV) per haplotype
  [4] normalize bcftools norm; GT = 1, 0 or . (haploid)
  [5] merge    one VCF, one sample column per HP value
  [6] annotate unique alleles once: VEP, ANNOVAR, InterVar, ClinVar,
               Taiwan Biobank, HGMD, OMIM (results cached in SQLite)
  [7] write    annotated VCF, MViewer table, figures (+ RCCX figure), timings.tsv
```

In the merged VCF `GT` is `1` when the haplotype carries the allele, `0` when its reads cover the site and show the reference, and `.` when its reads do not cover the site.

## Quick start

Requirements: Apptainer (or Docker), a GRCh38 FASTA with `.fai` and `.dict`, and about tens of GB of disk for the data directory.

The image holds open-source tools only (Clair3, Sawfish, minimap2, samtools, bcftools, VEP, BLAST). These cannot be redistributed, so you must obtain them yourself: ANNOVAR (with its `humandb`), InterVar (with `intervardb`), and the profile tables (ClinVar, Taiwan Biobank, HGMD, OMIM genemap2, and the small tables PanelList, DefaultTranscript, Comments, Methylated). Mount ANNOVAR at `/opt/annovar` and InterVar at `/opt/intervar` inside the container.

**1. Pull the image.** The GitHub package is private, so Apptainer needs a token with `read:packages`.

```bash
export APPTAINER_DOCKER_USERNAME=<github-user>
export APPTAINER_DOCKER_PASSWORD=<token>
apptainer pull paralog2viewer.sif docker://ghcr.io/oumarkme/paralog2viewer:0.9.0
```

**2. Build the data directory (one time).** `build-data` subsets every annotation source to the Paraphase loci plus padding (default 50 kb), rebuilds the indexes, precomputes RCCX discriminating sites, and writes `manifest.json` with checksums. The first build produced roughly 30 GB.

```bash
apptainer exec \
  -B /path/to/annovar:/opt/annovar,/path/to/intervar:/opt/intervar \
  -B /path/to/sources,/path/to/data \
  paralog2viewer.sif paralog2viewer build-data \
    --out /path/to/data \
    --reference /path/to/GRCh38.fa \
    --vep-cache /path/to/vep_cache \
    --annovar-db /path/to/annovar/humandb \
    --intervar-db /path/to/intervar/intervardb \
    --mviewer-profile /path/to/mviewer_profile \
    --regions /path/to/sample.paraphase.json \
    --threads 8
```

`--regions` takes one BED file or one or more Paraphase JSON files. The VEP cache root must contain `homo_sapiens_merged/` and `plugins/` (VEP plugin code, `alphamissense/`, `maxentscan/`).

**3. Run a sample.**

```bash
apptainer exec \
  -B /path/to/annovar:/opt/annovar,/path/to/intervar:/opt/intervar \
  -B /path/to/input,/path/to/data,/path/to/results \
  paralog2viewer.sif paralog2viewer run \
    --bam /path/to/input/SAMPLE.paraphase.bam \
    --data /path/to/data \
    --outdir /path/to/results
```

Use `--dry-run` first to check inputs and print the planned steps. `paralog2viewer check --data /path/to/data` verifies the data directory against its manifest (`--quick` compares sizes only).

## Outputs

Written to `<outdir>/<SAMPLE_ID>/`:

| File | Content |
|---|---|
| `{SAMPLE_ID}.paralog.vcf.gz`, `.tbi` | Merged, annotated VCF, one sample column per `HP` value |
| `{SAMPLE_ID}.paralog.mviewer.gz` | 365-column MViewer-format table |
| `figures/` | Read-pileup JPEG per region (or per locus when a region has haplotypes projected to different loci); RCCX module figure |
| `timings.tsv` | Wall time, CPU time, and peak memory per step |

`--format vcf,mviewer` (default) writes both. `--format vcf` or `--format mviewer` writes one. Annotation runs in every case; a single-format run skips only the file writing. `--no-figures` skips plotting.

`SAMPLE_ID` is the BAM file name without `.paraphase.bam`. The BAM must end in `.paraphase.bam` unless you pass `--sample`.

## Options of `run`

| Option | Default | Meaning |
|---|---|---|
| `--bam` | required | Paraphase BAM (adjacent `.bai` and `.json` required) |
| `--data` | `$PARALOG2VIEWER_DATA` | Data directory from `build-data` |
| `--reference` | copy in the data directory | GRCh38 FASTA with `.fai` and `.dict` |
| `--outdir` | `results` | Output directory |
| `--sample` | from BAM name | Sample ID |
| `--tmp` | none | Scratch directory. The tool makes its own subdirectory, uses only that, and removes it after a successful run. Without `--tmp`, task files stay in `./work/<sample>` |
| `--work-dir` | none | Reuse this work directory to resume a run (never deleted) |
| `--resume`, `--no-resume` | resume on | Skip steps with a valid done-marker |
| `--cache` | `<data>/annotation_cache.sqlite` | Annotation cache, keyed by allele and data manifest checksum |
| `--format` | `vcf,mviewer` | Comma list of output formats |
| `--figures`, `--no-figures` | figures on | Draw figures |
| `--vcf-annotation` | `default` | `full` also adds selected ANNOVAR fields to the VCF |
| `--pseudogene-tsv` | table shipped with the package | Pseudogene allele table |
| `--threads` | 40 | CPU budget in threads |
| `--cpus` | affinity unchanged | Pin the run to these logical CPUs, for example `0-19,48-67` |
| `--max-memory-gb` | 800 | Memory budget for task estimates (estimates, not measurements) |
| `--gpu` | off | Accepted, but currently has no effect (see below) |
| `--verify-checksums` | sizes only | Hash every data file at start (slow) |
| `--dry-run` | off | Check inputs, print the planned steps, exit |

`--gpu`: the Clair3 GPU option has not been verified against Clair3 2.0.3, so `GPU_FLAGS` in `call.py` is empty. With `--gpu` the run logs a warning and uses the CPU.

## Things to know about the table

- The columns are the standard 365 MViewer columns, in the order of `resources/mviewer_header.tsv`. The file starts with two comment lines (`#keywords=`, `#sex=`).
- One row per allele and per haplotype with GT 1. `Genotype` is the exact `HP` tag value.
- 58 columns are always empty: legacy columns, `Methylated loci`, CentoMD, PharmGKB, MaxEntScan, SpliceAI, family genotypes, Author Interpretation, Comment, Het_Percent, ExpansionHunter, AutoMap, `am_class`, `am_pathogenicity`.
- `Severity Score` is empty because its formula is not decoded.
- Nine columns are rebuilt by rule from VEP and `DefaultTranscript.tsv`: `Consequence`, `Gene`, `Variant`, `Type`, `Max Allele Freq`, `Transcript`, `Exon`, `Nucleotide` (HGVS c.), `AAChange` (HGVS p.).
- Transcript: a `DefaultTranscript.tsv` transcript is used whenever VEP reports one for the allele (its accession version is ignored). Otherwise the tool falls back to the MANE Select transcript, then the canonical RefSeq transcript, then any canonical transcript, then the first VEP transcript row.
- Symbolic and breakend alleles (Sawfish) keep their true ALT, coordinates, and VEP consequence. ANNOVAR-derived columns are empty for them.
- `Otherinfo1-13` come from the haplotype's own call. The `ID` field is `.`.
- The VCF carries the full VEP `CSQ` field and INFO tags `INTERVAR` and `CLNSIG`. It never contains HGMD or Taiwan Biobank values, because those data are licensed or restricted.
- Not verified: the HGMD, ClinVar, and Taiwan Biobank join rules were inferred from the output values of two earlier validation tables (9,901 rows), not from the vendors' documentation.

## RCCX notes

- Variants come from the same Clair3 and Sawfish path as every other gene.
- The module call (C4A/C4B, TNXA/TNXB, WHR1/WHR1B by a vote over precomputed BLAST sites; CYP21A2 versus CYP21A1P from Paraphase marker hits) depends on Paraphase phasing. If Paraphase phasing failed, alleles show as `unassigned`.
- The calls have had no clinical-geneticist review on a known CYP21A2 chimera. Treat them as research output.

## Resources and speed

- Default budget: 40 logical CPUs (`--threads`) and 800 GB of memory estimates (`--max-memory-gb`). Set lower values on a shared machine.
- The first full run of a validation sample took about 10 minutes on a 40-CPU budget. The previous Nextflow pipeline took about 33 minutes on the same sample.
- `--tmp` should point to a fast local disk. I have not benchmarked network storage against local disk for this tool.
- A finished step leaves a done-marker keyed by an input hash. A rerun with the same `--work-dir` skips it. The annotation cache also speeds up later samples.

## Troubleshooting

| Message or symptom | Action |
|---|---|
| Data manifest mismatch (build, tool versions, checksum) | Rebuild with `build-data`, or run `check --data DIR` to see which file differs. |
| Calling regions lie outside the data subsets | The data directory was built from other loci. Rebuild with `--regions` that include this sample's JSON, or a larger `--padding`. |
| Reference or dictionary missing | Give `--reference` a FASTA with `.fai` and `.dict` next to it. Its contigs must match the Paraphase BAM header. |
| Run failed | The work directory is kept. Start again with `--work-dir <same directory>` to resume. |
| ANNOVAR or InterVar not found | The tool looks for `/opt/annovar/table_annovar.pl` and `/opt/intervar/Intervar.py`. Check the `-B` mounts. |
| Input not accepted | The BAM must end in `.paraphase.bam` (or pass `--sample`) and have `.bai` and `.json` next to it. |

## Development

```bash
pip install -e '.[dev]'
pytest -q tests
ruff check && ruff format --check
```

The version lives in the `VERSION` file only. GitHub Actions (`.github/workflows/ci.yml`) builds the test image from `docker/Dockerfile`, runs ruff and pytest inside it, and on a push to `main` or a `v*` tag pushes the runtime image to `ghcr.io/oumarkme/paralog2viewer`. A tag must equal `v` plus the content of `VERSION`.

## License and third-party tools

paralog2viewer is GPL-3.0-or-later (see `LICENSE`). Third-party tools keep their own licenses. The table lists what a web lookup on 2026-11-13 found; "not verified" means no license text could be fetched, so check upstream before redistributing.

| Tool | Role | License | Citation |
|---|---|---|---|
| [Paraphase](https://github.com/PacificBiosciences/paraphase) | Haplotype phasing (input) | PacBio Software License Agreement (not OSI) | Chen et al., Nat Commun 2025, doi:10.1038/s41467-025-57505-2 |
| [Clair3](https://github.com/HKU-BAL/Clair3) | SNV/indel calling | BSD-3-Clause (README badge only) | Journal reference not verified (preprint doi:10.1101/2021.12.29.474431) |
| [Sawfish](https://github.com/PacificBiosciences/sawfish) | SV calling | not verified | Saunders et al., Bioinformatics 2025, doi:10.1093/bioinformatics/btaf136 |
| [minimap2](https://github.com/lh3/minimap2) | Alignment | MIT | Li, Bioinformatics 2018, doi:10.1093/bioinformatics/bty191 |
| [samtools, bcftools](https://github.com/samtools) | BAM and VCF tools | MIT (from a citation abstract) | Danecek et al., GigaScience 2021, doi:10.1093/gigascience/giab008 |
| [VEP](https://www.ensembl.org/info/docs/tools/vep/) | Consequence annotation | not verified | McLaren et al., Genome Biol 2016, doi:10.1186/s13059-016-0974-4 |
| [ANNOVAR](https://annovar.openbioinformatics.org/) | Annotation (user-supplied) | not verified; check the ANNOVAR site | Wang et al., Nucleic Acids Res 2010;38:e164 |
| [InterVar](https://github.com/WGLab/InterVar) | ACMG-style classification (user-supplied) | Free for non-commercial use | Li and Wang, Am J Hum Genet 2017, doi:10.1016/j.ajhg.2017.01.004 |
| BLAST+ | RCCX site design | not verified | Camacho et al., BMC Bioinformatics 2009, doi:10.1186/1471-2105-10-421 |

Citation: no paper describes paralog2viewer yet. If you publish results, cite the tools above. PubMed IDs were not checked.

History: v0.9.0 replaces an earlier Nextflow pipeline with this single Python tool; the old files are kept under `legacy/`.

## Author

Jen-Hsiang Ou, Ph.D., Department of Medical Genetics, National Taiwan University Hospital.

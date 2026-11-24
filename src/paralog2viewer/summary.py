# SPDX-License-Identifier: GPL-3.0-or-later
"""One-page HTML summary of a run (``paralog2viewer_summary.html``).

Standard library only. The page is self-contained: figures are embedded as
base64 (the exact JPEG files of ``figures/``) and set on the first toggle
of a gene fold. Nothing here changes the VCF, the table or the figures.

Inputs, all produced by the run: ``split_audit.json`` (genes, haplotypes,
reads), the Paraphase JSON (``total_cn``), the 365-column table, the
known-variant hits of ``plot.known_hits``, ``rccx_calls.tsv``, the figure
files and the failure messages of the run.
"""

from __future__ import annotations

import base64
import collections
import csv
import datetime
import gzip
import html
import json
import re
from collections.abc import Mapping, Sequence
from pathlib import Path

from paralog2viewer.annotate import PLP_CLINVAR, PLP_INTERVAR

SUMMARY_NAME = "paralog2viewer_summary.html"

# Flag rule (Broad). Every threshold of the page is defined here and
# printed in the footer.
RARE_AF = 0.01  # Max Allele Freq below this is rare; missing counts as 0
REVEL_MIN = 0.7  # VUS flagged above this REVEL score
CADD_MIN = 25.0  # VUS flagged above this CADD phred
HIGH_IMPACT = (
    "frameshift_variant",
    "stop_gained",
    "start_lost",
    "splice_donor_variant",
    "splice_acceptor_variant",
)
MAX_FLAG_ROWS = 300  # flagged rows listed per gene

_INTERVAR = re.compile(r"InterVar:\s*([A-Za-z ]+?)\s+PVS1", re.I)
_LABEL = {"ok": "OK", "review": "REVIEW", "problem": "PROBLEM"}
_SHAPE = {"ok": "●", "review": "▲", "problem": "■"}
# Evidence tiers of a flagged variant, most severe first. P is pure red
# (PI decision 2026-11-23, the one exception to Okabe-Ito); LP and note are
# Okabe-Ito orange and yellow. The tier is always also written.
TIERS = ("P", "LP", "note")
TIER_NAME = {"P": "Pathogenic", "LP": "Likely pathogenic", "note": "Note"}
_P_CLINVAR = frozenset({"Pathogenic", "Pathogenic/Likely_pathogenic"})
_LP_CLINVAR = frozenset({"Likely_pathogenic"})
# Locus role of the haplotype that carries the variant (plot.ROLE_ORDER).
ROLES = (
    "gene",
    "paralog",
    "pseudogene",
    "fusion",
    "other",
    "unresolved",
    "copy number",
)
_ROLE_OF = {"paralog_gene": "paralog"}
RCCX_METHOD = (
    "This section is not the Paraphase result itself. It comes from the "
    "in-house paralog2viewer RCCX module: Paraphase haplotypes and phasing "
    "are the input; gene type per haplotype is voted at the discriminating "
    "paralog sites (BLAST-derived); a haplotype Paraphase left unlinked is "
    "placed on the allele that lacks a distal (CYP21A2) end; modules are "
    "ordered 5' to 3' by type, and same-type neighbours stay unordered."
)


# ------------------------------------------------------------- flag rule


def _number(text: str | None) -> float | None:
    """Largest number in a cell (``0.5,0.7`` gives 0.7), None if none."""
    values = []
    for part in re.split(r"[,;|]", text or ""):
        try:
            values.append(float(part))
        except ValueError:
            continue
    return max(values) if values else None


def intervar_class(text: str) -> str:
    """Return the ACMG class text of an InterVar cell ("" if absent)."""
    match = _INTERVAR.search(text or "")
    return match.group(1).strip() if match else ""


def clinvar_plp(row: Mapping[str, str]) -> bool:
    """True if ClinVar says P/LP; same sets as ``annotate.is_pathogenic``."""
    values = (row.get("ClinVar", ""), row.get("ClinVar_SIG", ""))
    return any(v.strip().replace(" ", "_") in PLP_CLINVAR for v in values)


def intervar_plp(row: Mapping[str, str]) -> bool:
    """True if InterVar says Pathogenic or Likely pathogenic."""
    return intervar_class(row.get("InterVar", "")) in PLP_INTERVAR


def plp_tier(row: Mapping[str, str]) -> str:
    """Return "P", "LP" or "" from ClinVar and InterVar (higher wins)."""
    values = {
        v.strip().replace(" ", "_")
        for v in (row.get("ClinVar", ""), row.get("ClinVar_SIG", ""))
    }
    inter = intervar_class(row.get("InterVar", ""))
    if values & _P_CLINVAR or inter == "Pathogenic":
        return "P"
    if values & _LP_CLINVAR or inter == "Likely pathogenic":
        return "LP"
    return ""


def flag_reasons(row: Mapping[str, str]) -> list[str]:
    """Return why a table row is flagged (empty list: not flagged).

    Known-variant hits are added by the caller. Rules: ClinVar or InterVar
    Pathogenic/Likely pathogenic; high-impact consequence with Max Allele
    Freq < RARE_AF; VUS (InterVar or ClinVar) with REVEL > REVEL_MIN or
    CADD_phred > CADD_MIN and Max Allele Freq < RARE_AF.
    """
    reasons = []
    clinvar = row.get("ClinVar", "")
    inter = intervar_class(row.get("InterVar", ""))
    if clinvar_plp(row):
        reasons.append(
            f"ClinVar {clinvar.strip() or row.get('ClinVar_SIG', '')}"
        )
    if intervar_plp(row):
        reasons.append(f"InterVar {inter}")
    freq = _number(row.get("Max Allele Freq"))
    freq = 0.0 if freq is None else freq
    rare = freq < RARE_AF
    terms = set(re.split(r"[;,&]", row.get("Consequence", "")))
    impact = sorted(terms & set(HIGH_IMPACT))
    if impact and rare:
        reasons.append(f"high impact {', '.join(impact)}; AF {freq:g}")
    vus = "uncertain" in inter.lower() or "uncertain" in clinvar.lower()
    if vus and rare:
        revel = _number(row.get("REVEL_score"))
        cadd = _number(row.get("CADD_phred"))
        if revel is not None and revel > REVEL_MIN:
            reasons.append(f"VUS REVEL {revel:g}; AF {freq:g}")
        if cadd is not None and cadd > CADD_MIN:
            reasons.append(f"VUS CADD {cadd:g}; AF {freq:g}")
    return reasons


def _minimal(pos: int, ref: str, alt: str) -> tuple[int, str, str]:
    """Normalise ANNOVAR-style and VCF-style alleles to one form.

    ANNOVAR order: shared leading bases are trimmed first (the position
    moves right), then shared trailing bases, so REF AGA ALT A and the
    ANNOVAR row GA/- give the same result. ``-`` is empty; an empty REF
    (ANNOVAR insertion) moves POS one base right.
    """
    if ref == "-":
        return pos + 1, "", alt.upper()
    ref = ref.upper()
    alt = "" if alt == "-" else alt.upper()
    shift = 0
    while ref and alt and ref[0] == alt[0]:
        ref, alt, shift = ref[1:], alt[1:], shift + 1
    while ref and alt and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    return pos + shift, ref, alt


def _site(chrom: str, pos: int, ref: str, alt: str) -> tuple:
    return (chrom.removeprefix("chr"), *_minimal(pos, ref, alt))


# ----------------------------------------------------------------- inputs


def load_table(path: Path | str | None) -> list[dict[str, str]]:
    """Read the 365-column table (gzip or plain; ``#`` lines skipped)."""
    if not path or not Path(path).is_file():
        return []
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        lines = (ln for ln in handle if not ln.startswith("#"))
        return list(csv.DictReader(lines, delimiter="\t"))


def _read_json(path: Path) -> dict:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return data if isinstance(data, dict) else {}


def _paraphase_cn(paraphase_json: Path | None) -> dict[str, tuple]:
    """Region -> (expected total_cn or None, haplotype copies resolved).

    Copies resolved count a haplotype listed in ``two_copy_haplotypes``
    twice.
    """
    data = _read_json(paraphase_json) if paraphase_json else {}
    out = {}
    for region, rec in data.items():
        if not isinstance(rec, dict):
            continue
        final = rec.get("final_haplotypes")
        final = final if isinstance(final, dict) else {}
        names = set(final.values())
        two = {final.get(k, k) for k in rec.get("two_copy_haplotypes") or ()}
        expected = rec.get("total_cn")
        expected = (
            expected
            if isinstance(expected, int) and not isinstance(expected, bool)
            else None
        )
        out[region] = (expected, len(names) + len(two & names))
    return out


def _region_of(hap: str, region_of: Mapping[str, str]) -> str:
    return region_of.get(hap) or (hap.rsplit("_", 1)[0] if "_" in hap else hap)


# ------------------------------------------------------------- the model


def build_model(
    *,
    sample: str,
    version: str,
    work_dir: Path,
    failures: Mapping[str, str],
    table_path: Path | str | None = None,
    paraphase_json: Path | None = None,
    known_hits: Sequence[tuple] = (),
    hap_symbol: Mapping[str, str] | None = None,
    hap_role: Mapping[str, str] | None = None,
    figures_dir: Path | None = None,
    generated: datetime.datetime | None = None,
) -> dict:
    """Collect everything the page shows (no HTML yet).

    Args:
        sample: Sample label.
        version: Tool version.
        work_dir: Run work directory (``split_audit.json``,
            ``rccx_calls.tsv``).
        failures: ``RunResult.failures`` (step or ``haplotype <HP>`` ->
            message).
        table_path: Published 365-column table; None gives no variants.
        paraphase_json: Paraphase JSON (``total_cn`` per region).
        known_hits: ``plot.known_hits`` tuples (gene, chrom, pos, ref, alt,
            label, hgvs, reference).
        hap_symbol: HP -> gene symbol; a hit applies to rows of
            haplotypes whose symbol equals the hit gene (case-insensitive).
        hap_role: HP -> role (``plot.ROLE_ORDER`` value) of the carrying
            haplotype; a haplotype without an entry is "unresolved", except
            that a haplotype in ``hap_symbol`` is a "gene".
        figures_dir: ``figures/`` of the sample.
        generated: Timestamp shown in the header (default: now, local time).

    Returns:
        Dict with ``genes`` (sorted, problems first), ``problems`` (run
        level lines), ``overall``, ``rccx`` and header fields.
    """
    audit = _read_json(Path(work_dir) / "split_audit.json")
    haps = [h for h in audit.get("haplotypes", []) if isinstance(h, dict)]
    region_of = {h.get("hap"): h.get("region") for h in haps}
    cn = _paraphase_cn(paraphase_json)
    hap_symbol = hap_symbol or {}
    roles = {
        hap: _ROLE_OF.get(role, role) for hap, role in (hap_role or {}).items()
    }
    for hap in hap_symbol:
        roles.setdefault(hap, "gene")
    problems: list[str] = []
    if not audit:
        problems.append("split_audit.json missing or unreadable")

    genes: dict[str, dict] = {}
    no_haps: list[str] = []
    for rec in audit.get("genes", []):
        name = rec.get("gene", "")
        if rec.get("status") != "resolved":
            no_haps.append(name)
            continue
        expected, copies = cn.get(name, (None, 0))
        genes[name] = {
            "name": name,
            "resolved": rec.get("resolved_haplotypes", 0),
            "copies": copies,
            "expected": expected,
            "reads": 0,
            "notes": [],
            "status": "ok",
            "variants": 0,
            "flagged": [],
            "tiers": collections.Counter(),
            "figures": [],
            "haps": [],
        }
        if expected is not None and expected != copies:
            genes[name]["status"] = "review"
            genes[name]["notes"].append(
                f"CN {copies} resolved vs {expected} expected"
            )
            genes[name]["tiers"][("note", "copy number")] += 1
    for h in haps:
        gene = genes.get(h.get("region"))
        if gene is None:
            continue
        gene["haps"].append(h.get("hap"))
        reads = h.get("matched_primary_read_ids")
        gene["reads"] += reads if isinstance(reads, int) else 0
        if h.get("status") != "ready":
            gene["notes"].append(f"{h.get('hap')}: {h.get('status')}")
            if gene["status"] == "ok":
                gene["status"] = "review"

    steps = []
    for key, message in failures.items():
        text = f"{key}: {str(message).strip().splitlines()[0][:200]}"
        hap = (
            key.removeprefix("haplotype ")
            if key.startswith("haplotype ")
            else ""
        )
        gene = genes.get(_region_of(hap, region_of)) if hap else None
        if gene is not None:
            gene["status"] = "problem"
            gene["notes"].append(f"failed: {text}")
        steps.append(text)

    _attach_variants(
        genes, table_path, region_of, known_hits, hap_symbol, roles
    )
    _attach_figures(genes, figures_dir)

    def _n(gene: Mapping, tier: str) -> int:
        return sum(n for (t, _r), n in gene["tiers"].items() if t == tier)

    # P, then LP, then notes (flag and CN notes, failures), then name.
    ordered = sorted(
        genes.values(),
        key=lambda g: (
            -_n(g, "P"),
            -_n(g, "LP"),
            -(
                _n(g, "note")
                + sum(not n.startswith("CN ") for n in g["notes"])
            ),
            g["name"],
        ),
    )
    rccx = _read_rccx(Path(work_dir) / "rccx_calls.tsv")
    if steps:
        problems.append("failed: " + "; ".join(steps))
    if any(g["status"] == "problem" for g in ordered) or problems:
        overall = "problem"
    elif any(g["status"] == "review" for g in ordered):
        overall = "review"
    else:
        overall = "ok"
    reviews = [
        f"{g['name']}: {n}"
        for g in ordered
        for n in g["notes"]
        if n.startswith("CN ")
    ]
    return {
        "sample": sample,
        "version": version,
        "generated": (
            generated or datetime.datetime.now().astimezone()
        ).strftime("%Y-%m-%d %H:%M"),
        "genes": ordered,
        "problems": problems,
        "reviews": reviews,
        "no_haps": sorted(no_haps),
        "overall": overall,
        "rccx": rccx,
        "figures_dir": figures_dir,
    }


def _attach_variants(
    genes: dict[str, dict],
    table_path: Path | str | None,
    region_of: Mapping[str, str],
    known_hits: Sequence[tuple],
    hap_symbol: Mapping[str, str],
    roles: Mapping[str, str],
) -> None:
    """Count variants and list flagged ones per gene (in place).

    A flagged row gets ``_tier`` (P, LP or note) and ``_role`` (locus role
    of its haplotype). ``tiers`` counts distinct (variant, haplotype) pairs, so
    a variant on three haplotypes counts three times.
    """
    hits = collections.defaultdict(list)
    for gene, chrom, pos, ref, alt, label, hgvs, reference in known_hits:
        hits[_site(chrom, int(pos), ref, alt)].append(
            (gene.upper(), f"known variant {label} ({hgvs}; {reference})")
        )
    counted: set[tuple] = set()
    for row in load_table(table_path):
        hap = row.get("Genotype", "")
        name = _region_of(hap, region_of)
        gene = genes.get(name)
        if gene is None:
            continue
        gene["variants"] += 1
        reasons = flag_reasons(row)
        try:
            site = _site(
                row.get("Chr", ""),
                int(row.get("Start", "")),
                row.get("Ref", ""),
                row.get("Alt", ""),
            )
        except ValueError:
            site = None
        symbol = hap_symbol.get(hap, "").upper()
        reasons += [t for g, t in hits.get(site, ()) if g == symbol]
        if not reasons:
            continue
        tier = plp_tier(row) or "note"
        role = roles.get(hap, "unresolved")
        gene["flagged"].append(
            {**row, "_reasons": reasons, "_tier": tier, "_role": role}
        )
        key = (name, site or hap, hap, tier, role)
        if key not in counted:
            counted.add(key)
            gene["tiers"][(tier, role)] += 1
    for gene in genes.values():
        gene["flagged"].sort(
            key=lambda r: (
                TIERS.index(r["_tier"]),
                ROLES.index(r["_role"]) if r["_role"] in ROLES else 99,
            )
        )


def _attach_figures(genes: dict[str, dict], figures_dir: Path | None) -> None:
    """Assign ``figures/*.jpeg`` to genes by file-name prefix (in place)."""
    if not figures_dir or not Path(figures_dir).is_dir():
        return
    names = sorted(genes, key=len, reverse=True)
    for path in sorted(Path(figures_dir).glob("*.jpeg")):
        stem = path.stem
        for name in names:
            if stem == name or stem.startswith(name + "_"):
                genes[name]["figures"].append(path)
                break


def _read_rccx(path: Path) -> list[dict[str, str]] | None:
    """Return the RCCX call rows, or None if the run had no RCCX table."""
    if not path.is_file():
        return None
    with open(path, encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


# ------------------------------------------------------------------ HTML

_CSS = """
:root{color-scheme:light;--bg:#fff;--fg:#1a1a1a;--mut:#595959;
--line:#d0d0d0;--card:#f5f5f5;--ok:#009E73;--review:#E69F00;
--problem:#D55E00;--info:#0072B2;--on:#fff;--red:#D55E00}
*{box-sizing:border-box}
body{margin:0 auto;width:96%;padding:.75rem 0;background:var(--bg);
color:var(--fg);font:15px/1.4 system-ui,sans-serif}
h1{font-size:1.25rem;margin:.2rem 0}h2{font-size:1.05rem;margin:1.2rem 0 .4rem}
.meta{color:var(--mut);font-size:.85rem}
.bar{border:1px solid var(--line);border-left:.5rem solid var(--c);
background:var(--card);padding:.6rem .8rem;margin:.8rem 0;border-radius:.3rem}
.bar ul{margin:.3rem 0 0;padding-left:1.1rem}
.b{display:inline-block;padding:.05rem .45rem;border-radius:.8rem;
font-weight:600;font-size:.8rem;background:var(--c);color:#000}
.red{color:var(--red);font-weight:600}
.b.review,.b.problem{background:transparent;color:var(--red);
border:1px solid var(--red)}
.t{display:inline-block;padding:.05rem .5rem;border-radius:.8rem;
font:inherit;font-weight:600;font-size:.8rem;color:#000;border:1px solid #000;
white-space:nowrap;margin-right:.25rem}
button.t{cursor:pointer}button.t[aria-pressed=true]{outline:2px solid var(--fg);
outline-offset:1px}
.t.P{background:#CC0000;color:#fff}.t.LP{background:#E69F00}.t.note{background:#F0E442}
td.fl,td.nw{white-space:nowrap}th.cg{width:11rem}th.cf{width:12rem}
th.cn{width:8rem}td.fi{color:var(--fg)}tr.g.none td{color:var(--mut)}tr.g.none td.red{color:var(--red)}tr.g{cursor:pointer}tr.g:hover{background:var(--card)}
tr.d>td{background:var(--card);padding:.4rem .6rem}
.ok{--c:var(--ok)}.review{--c:var(--review)}.problem{--c:var(--problem)}
.info{--c:var(--info)}
.tw{overflow-x:auto}
table{border-collapse:collapse;width:100%;font-size:.85rem}
th,td{border-bottom:1px solid var(--line);padding:.25rem .4rem;
overflow-wrap:break-word;
text-align:left;vertical-align:top}
th{position:sticky;top:0;background:var(--card)}
td.n,th.n{text-align:right;font-variant-numeric:tabular-nums}
details{border:1px solid var(--line);border-radius:.3rem;margin:.4rem 0}
summary{cursor:pointer;padding:.4rem .6rem;display:flex;gap:.5rem;
flex-wrap:wrap;align-items:center}
details>div{padding:0 .6rem .6rem}
figure{margin:.6rem 0}figure img{max-width:100%;height:auto;background:#fff}
.btn{display:inline-block;border:1px solid var(--info);color:var(--info);
padding:.1rem .6rem;border-radius:.3rem;text-decoration:none;font-size:.8rem}
.note{color:var(--mut);font-size:.8rem}
[hidden]{display:none!important}
footer{margin-top:1.5rem;border-top:1px solid var(--line);padding-top:.6rem}
"""

_JS = """
(function(){
function load(d){d.querySelectorAll('img[data-src]').forEach(function(i){
if(!i.getAttribute('src')){i.src=i.dataset.src;}});}
document.querySelectorAll('details.gene').forEach(function(d){
d.addEventListener('toggle',function(){if(d.open){load(d);}});});
document.querySelectorAll('a.dl').forEach(function(a){
a.addEventListener('click',function(){
a.href=document.getElementById(a.dataset.img).dataset.src;});});
var rows=[].slice.call(document.querySelectorAll('tr.g'));
var tops=[].slice.call(document.querySelectorAll('.bar button.t'));
var global=null;
function show(r,tier,open){
var d=r.nextElementSibling;d.hidden=!open;r.dataset.open=open?'1':'';
r.dataset.tier=tier||'';
d.querySelectorAll('tr.v').forEach(function(v){
v.hidden=!!tier&&v.dataset.tier!==tier;});
r.querySelectorAll('button.t').forEach(function(b){
b.setAttribute('aria-pressed',open&&b.dataset.tier===tier);});
if(open){load(d);}}
function apply(){
tops.forEach(function(b){b.setAttribute('aria-pressed',b.dataset.tier===global);});
rows.forEach(function(r){
var has=!global||r.dataset.tiers.split(' ').indexOf(global)>=0;
r.hidden=!has;
if(!has){r.nextElementSibling.hidden=true;}
else{show(r,global,!!global);}});}
tops.forEach(function(b){b.addEventListener('click',function(){
global=global===b.dataset.tier?null:b.dataset.tier;apply();});});
rows.forEach(function(r){
r.addEventListener('click',function(e){
var b=e.target.closest('button.t');
if(b){var same=r.dataset.open&&r.dataset.tier===b.dataset.tier;
show(r,b.dataset.tier,!same);return;}
if(e.target.closest('a'))return;
show(r,global,!r.dataset.open);});});
})();
"""


def _e(value: object) -> str:
    return html.escape(str(value), quote=True)


def _red(text: object) -> str:
    """Escape ``text`` and mark it red (one CSS token, ``--red``)."""
    return f'<span class="red">{_e(text)}</span>'


def _cell(text: object, *, red: bool = False, cls: str = "") -> str:
    """Return a table cell; ``red`` marks an abnormal or to-check value."""
    classes = " ".join(c for c in (cls, "red" if red else "") if c)
    attr = f' class="{classes}"' if classes else ""
    return f"<td{attr}>{_e(text)}</td>"


def _badge(status: str) -> str:
    return f'<span class="b {status}">{_SHAPE[status]} {_LABEL[status]}</span>'


def _tier_counts(counter: Mapping[tuple[str, str], int]) -> dict[str, int]:
    return {
        tier: sum(n for (t, _r), n in counter.items() if t == tier)
        for tier in TIERS
    }


def _role_text(counter: Mapping[tuple[str, str], int], tier: str) -> str:
    parts = [
        f"{n} {role}"
        for role in ROLES
        if (n := counter.get((tier, role), 0))
    ]
    return ", ".join(parts)


def _tags(counter: Mapping[tuple[str, str], int]) -> str:
    """Return clickable tags ("3 P", "2 LP", "1 note"; "none" if empty).

    The tooltip gives the tier name and the locus roles that make it up.
    """
    totals = _tier_counts(counter)
    chips = [
        f'<button type="button" class="t {t}" data-tier="{t}" '
        f'aria-pressed="false" title="{TIER_NAME[t]}: '
        f'{_e(_role_text(counter, t))}">{n} {t}</button>'
        for t, n in totals.items()
        if n
    ]
    return "".join(chips) or '<span class="note">none</span>'


def _allele_text(allele: str, keep: int = 8) -> str:
    """Shorten a long allele: ``CCTCTTGC...(95 bp)``."""
    if len(allele) <= keep + 4:
        return allele
    return f"{allele[:keep]}…({len(allele)} bp)"


def _variant_label(row: Mapping[str, str]) -> str:
    return (
        f"chr{row.get('Chr', '')}:{row.get('Start', '')} "
        f"{_allele_text(row.get('Ref', ''))}>{_allele_text(row.get('Alt', ''))}"
    )


def _aa_short(row: Mapping[str, str]) -> str:
    """First transcript's ``c.`` and ``p.`` change ("" if none listed)."""
    text = row.get("AAChange") or row.get("AAChange.refgene", "")
    first = re.split(r"[,;]", text)[0]
    parts = [x for x in first.split(":") if x.startswith(("c.", "p."))]
    return " ".join(parts)


def _flag_table(gene: Mapping) -> str:
    """Flagged variants of one gene; one line per variant and locus.

    The same variant on several haplotypes of one locus role is one row
    with the haplotype names joined.
    """
    cols = (
        "Tier",
        "Locus",
        "Variant",
        "Gene",
        "Change",
        "Haplotypes",
        "AF",
        "REVEL/CADD",
        "Evidence",
    )
    groups: dict[tuple, dict] = {}
    for row in gene["flagged"]:
        key = (
            row["_tier"],
            row["_role"],
            _variant_label(row),
            tuple(row["_reasons"]),
        )
        group = groups.setdefault(key, {"row": row, "haps": []})
        group["haps"].append(row.get("Genotype", ""))
    rows = []
    for (tier, role, label, reasons), group in list(groups.items())[
        :MAX_FLAG_ROWS
    ]:
        row = group["row"]
        full = row.get("AAChange") or row.get("AAChange.refgene", "")
        cells = (
            f'<td class="nw"><span class="t {tier}">{tier}</span></td>',
            _cell(role, cls="nw"),
            _cell(label, red=True, cls="nw"),
            _cell(row.get("Gene", ""), cls="nw"),
            f'<td title="{_e(full)}">{_e(_aa_short(row))}</td>',
            _cell(
                ", ".join(
                    h.removeprefix(gene["name"] + "_") for h in group["haps"]
                ),
                cls="nw",
            ),
            _cell(row.get("Max Allele Freq", ""), cls="nw"),
            _cell(
                f"{row.get('REVEL_score', '')}/{row.get('CADD_phred', '')}",
                cls="nw",
            ),
            _cell("; ".join(reasons)),
        )
        rows.append(
            f'<tr class="v" data-tier="{tier}">' + "".join(cells) + "</tr>"
        )
    more = len(groups) - MAX_FLAG_ROWS
    tail = (
        f'<p class="note">{more} more flagged rows in the table.</p>'
        if more > 0
        else ""
    )
    head = "".join(f"<th>{c}</th>" for c in cols)
    body = "".join(rows)
    return f'<div class="tw"><table><tr>{head}</tr>{body}</table></div>{tail}'


def _findings(gene: Mapping, limit: int = 4) -> str:
    """Short P/LP variant list for the gene row ("+n" beyond ``limit``)."""
    seen: dict[str, None] = {}
    for row in gene["flagged"]:
        if row["_tier"] in ("P", "LP"):
            seen.setdefault(_aa_short(row) or _variant_label(row))
    names = list(seen)
    text = ", ".join(names[:limit])
    return text + (f" +{len(names) - limit}" if len(names) > limit else "")


def _tiers_attr(gene: Mapping) -> str:
    totals = _tier_counts(gene["tiers"])
    return " ".join(t for t in TIERS if totals[t])


def _figure(path: Path, ident: str, lazy: bool = True) -> str:
    data = base64.b64encode(path.read_bytes()).decode("ascii")
    attr = "data-src" if lazy else "src"
    return (
        f'<figure><img id="{ident}" alt="{_e(path.stem)}" loading="lazy" '
        f'{attr}="data:image/jpeg;base64,{data}">'
        f'<figcaption><a class="btn dl" href="#" download="{_e(path.name)}" '
        f'data-img="{ident}">Download {_e(path.name)}</a></figcaption></figure>'
    )


def _gene_rows(index: int, gene: Mapping) -> str:
    """One table row per gene plus its hidden detail row."""
    body = []
    if gene["notes"]:
        body.append(
            "<ul>"
            + "".join(f"<li>{_red(n)}</li>" for n in gene["notes"])
            + "</ul>"
        )
    if gene["flagged"]:
        body.append(_flag_table(gene))
    other = gene["variants"] - len(gene["flagged"])
    body.append(f'<p class="note">{other} other variants in the table.</p>')
    for number, path in enumerate(gene["figures"]):
        body.append(_figure(path, f"f{index}_{number}"))
    if not gene["figures"]:
        body.append('<p class="note">No figure.</p>')
    expected = "-" if gene["expected"] is None else gene["expected"]
    cn_notes = [n for n in gene["notes"] if n.startswith("CN ")]
    cn_bad = bool(cn_notes)
    cn_text = _e(f"{gene['copies']} / {expected}")
    cn_pill = (
        f'<span class="t note" title="{_e(cn_notes[0])}">{cn_text}</span>'
        if cn_bad
        else cn_text
    )
    problem = _badge("problem") + " " if gene["status"] == "problem" else ""
    return (
        f'<tr class="g{"" if gene["flagged"] else " none"}" '
        f'data-tiers="{_tiers_attr(gene)}">'
        f"<td>{problem}<strong>{_e(gene['name'])}</strong></td>"
        f"<td class=\"fl\">{_tags(gene['tiers'])}</td>"
        f"<td class=\"fi\">{_e(_findings(gene))}</td>"
        f"<td class=\"n\">{cn_pill}</td>"
        f"{_cell(gene['variants'], cls='n')}</tr>"
        f'<tr class="d" hidden><td colspan="5">{"".join(body)}</td></tr>'
    )


def _no_haps_section(model: Mapping) -> str:
    names = model["no_haps"]
    if not names:
        return ""
    return (
        f"<details><summary>{_red(f'{len(names)} genes without haplotypes')}"
        f'</summary><div class="note">{_e(", ".join(names))}</div></details>'
    )


def _rccx_section(model: Mapping) -> str:
    rows = model["rccx"]
    if rows is None:
        return ""
    out = [
        "<h2>RCCX</h2>",
        f'<div class="bar info"><span class="b info">IN-HOUSE PIPELINE</span> '
        f"{_e(RCCX_METHOD)}</div>",
    ]
    figure = Path(model["figures_dir"] or "") / "RCCX_module.jpeg"
    if model["figures_dir"] and figure.is_file():
        out.append(_figure(figure, "frccx", lazy=False))
    return "".join(out)


def render_html(model: Mapping) -> str:
    """Render the page from ``build_model`` output."""
    genes = model["genes"]
    total = collections.Counter()
    for g in genes:
        total.update(g["tiers"])
    overall = "problem" if model["overall"] == "problem" else "info"
    seen = list(dict.fromkeys(model["problems"] + model["reviews"]))
    if model["no_haps"]:
        seen.append(
            f"no haplotypes: {len(model['no_haps'])} genes (listed below)"
        )
    bar_items = "".join(f"<li>{_red(x)}</li>" for x in seen)
    rows = "".join(_gene_rows(i, g) for i, g in enumerate(genes))
    const = (
        "Tier: P and LP = ClinVar or InterVar Pathogenic or Likely "
        "pathogenic (the higher of the two); note = any other flag: known "
        "gene/pseudogene variant on a carrier haplotype of the functional "
        f"gene; high impact ({', '.join(HIGH_IMPACT)}) with Max Allele Freq "
        f"&lt; {RARE_AF:g}; VUS with REVEL &gt; {REVEL_MIN:g} or CADD_phred "
        f"&gt; {CADD_MIN:g} and Max Allele Freq &lt; {RARE_AF:g}. Missing "
        "Max Allele Freq counts as 0. Locus = role of the haplotype that "
        "carries the variant (gene, pseudogene, ...; unresolved if the "
        "haplotype could not be assigned). Counts are distinct variants per "
        f"locus. Max {MAX_FLAG_ROWS} flagged rows per gene. CN: haplotypes "
        "resolved (two-copy haplotypes counted twice) vs Paraphase total_cn."
    )
    sample = _e(model["sample"])
    return (
        '<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width,initial-scale=1">'
        f"<title>paralog2viewer {sample}</title><style>{_CSS}</style></head>"
        f"<body><h1>{sample}</h1>"
        f'<div class="meta">paralog2viewer v{_e(model["version"])} · '
        f"{_e(model['generated'])} · GRCh38</div>"
        f'<div class="bar {overall}">{len(genes)} genes: {_tags(total)} '
        '<span class="note">Click a tag to list those variants; click a gene '
        "for its variants and figures.</span>"
        f"{'<ul>' + bar_items + '</ul>' if bar_items else ''}</div>"
        '<div class="tw"><table><tr><th class="cg">Gene</th>'
        '<th class="cf">Flags</th><th>P / LP findings</th>'
        '<th class="n cn">CN / expected</th>'
        '<th class="n cn">Variants</th></tr>'
        f"{rows}</table></div>"
        f"{_no_haps_section(model)}{_rccx_section(model)}"
        f'<footer class="note">{const}<br>Files: {sample}.paralog.vcf.gz, '
        f"{sample}.paralog.tsv.gz. Figures are embedded as JPEG (150 dpi); "
        "each has a Download link."
        f"</footer><script>{_JS}</script></body></html>"
    )


def write_summary(out_path: Path | str, **kwargs: object) -> Path:
    """Write the summary page; arguments as in ``build_model``.

    Returns:
        ``out_path``.
    """
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name(out.name + ".part")
    tmp.write_text(render_html(build_model(**kwargs)), encoding="utf-8")
    tmp.replace(out)
    return out

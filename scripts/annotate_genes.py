#!/usr/bin/env python3
"""Gene-level rNMP / hotspot / mutation overlay for sacCer3.

Does NOT claim human FTD or tau biology from yeast. Overlays are ortholog
lookups plus gene-set membership flags.

Inputs (optional, skip if missing):
  results/<GSE>/<SRR>/rnmp_hotspots.tsv
  results/<GSE>/<SRR>/rnmp_signal.tsv
  catalogs/ftd_aging_tau.tsv
  catalogs/one_carbon_energy.tsv
  data/refs/sacCer3.sgd.gff  (downloaded if absent)

Outputs:
  results/_cross/gene_rnmp.tsv
  results/_cross/failed_excision.tsv
  results/_cross/stall_mutation_genes.tsv
  results/_cross/ftd_tau_overlay.tsv
  results/_cross/gene_questions.html
"""
from __future__ import annotations

import csv
import os
import urllib.request
from collections import defaultdict
from pathlib import Path

REPO = Path(os.environ.get("REPO", "/home/noahtjones/hu-rnmp-pipeline"))
RESULTS = Path(os.environ.get("RESULTS", str(REPO / "results")))
REF = Path(os.environ.get("REF_DIR", str(REPO / "data/refs")))
CAT = REPO / "catalogs"
OUT = RESULTS / "_cross"
GFF_URL = "https://ftp.ensembl.org/pub/release-113/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.113.gff3.gz"

HYDEN_NT = "SRR8503261"
HYDEN_HII = "SRR8503270"
HU = "SRR24871232"
C0 = "SRR24871236"

CHROM_ALIAS = {
    "i": "chrI", "ii": "chrII", "iii": "chrIII", "iv": "chrIV", "v": "chrV",
    "vi": "chrVI", "vii": "chrVII", "viii": "chrVIII", "ix": "chrIX", "x": "chrX",
    "xi": "chrXI", "xii": "chrXII", "xiii": "chrXIII", "xiv": "chrXIV",
    "xv": "chrXV", "xvi": "chrXVI", "m": "chrmt", "mt": "chrmt",
    "chrmito": "chrmt", "mito": "chrmt", "mitochondrion": "chrmt",
    "2-micron": "2micron", "2micron": "2micron",
}


def norm_chrom(c: str) -> str:
    s = (c or "").strip()
    if not s:
        return s
    key = s.lower().replace("chromosome_", "").replace("chromosome", "").replace("chr", "")
    key = key.replace("_", "").strip()
    if s.lower() in ("chrmt", "chrM", "chrmito"):
        return "chrmt"
    return CHROM_ALIAS.get(key, s if s.startswith("chr") else f"chr{s}")


def ensure_gff() -> Path:
    dest = REF / "sacCer3.ensembl.gff3"
    gz = dest.with_suffix(dest.suffix + ".gz")
    if dest.exists() and dest.stat().st_size > 100_000:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"==> download Ensembl GFF3 → {gz}")
    req = urllib.request.Request(GFF_URL, headers={"User-Agent": "hu-rnmp-pipeline/wave3"})
    with urllib.request.urlopen(req, timeout=180) as r, gz.open("wb") as f:
        f.write(r.read())
    import gzip
    import shutil
    with gzip.open(gz, "rb") as src, dest.open("wb") as out:
        shutil.copyfileobj(src, out)
    return dest


def parse_gff(path: Path) -> list[dict]:
    genes = []
    with path.open() as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                if line.startswith("##FASTA"):
                    break
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] not in {"gene", "ncRNA_gene", "rRNA_gene", "tRNA_gene", "snRNA_gene", "snoRNA_gene"}:
                continue
            attrs = {}
            for kv in p[8].split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    attrs[k] = v
            sys = (attrs.get("gene_id") or attrs.get("ID") or "").replace("gene:", "")
            name = attrs.get("Name") or attrs.get("gene") or attrs.get("locus_tag") or sys
            genes.append({
                "chrom": norm_chrom(p[0]),
                "start": int(p[3]) - 1,
                "end": int(p[4]),
                "strand": p[6],
                "systematic": sys,
                "symbol": name,
                "biotype": p[2],
            })
    return genes


def load_catalog(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def iter_result_files(name: str):
    for p in RESULTS.glob(f"*/*/{name}"):
        if "report_bundle" in p.parts:
            continue
        yield p.parent.name, p


def load_hotspots():
    out = []
    for acc, p in iter_result_files("rnmp_hotspots.tsv"):
        with p.open() as f:
            r = csv.DictReader(f, delimiter="\t")
            for row in r:
                try:
                    out.append({
                        "run": acc,
                        "chrom": norm_chrom(row.get("chrom") or row.get("chr") or ""),
                        "start": int(float(row.get("start") or 0)),
                        "end": int(float(row.get("end") or 0)),
                        "total": float(row.get("total") or 0),
                        "fold": float(row.get("fold_over_median") or 0),
                    })
                except (TypeError, ValueError):
                    continue
    return out


def load_signals():
    """WGS rNMP-transition candidates (chrom, pos)."""
    out = []
    for acc, p in iter_result_files("rnmp_signal.tsv"):
        with p.open() as f:
            r = csv.DictReader(f, delimiter="\t")
            fields = r.fieldnames or []
            for row in r:
                chrom = norm_chrom(row.get("chrom") or row.get("chr") or "")
                pos = row.get("pos") or row.get("position") or row.get("start")
                try:
                    pos_i = int(float(pos))
                except (TypeError, ValueError):
                    continue
                out.append({
                    "run": acc,
                    "chrom": chrom,
                    "start": pos_i,
                    "end": pos_i + 1,
                    "kind": row.get("type") or row.get("alt") or "signal",
                })
        # cap per-run for memory if someone dumped genome-wide pileup
        if len(out) > 2_000_000:
            break
    return out


def load_mutations():
    out = []
    for acc, p in iter_result_files("rnmp_mutations.tsv"):
        with p.open() as f:
            r = csv.DictReader(f, delimiter="\t")
            for row in r:
                chrom = norm_chrom(row.get("chrom") or row.get("chr") or "")
                pos = row.get("pos") or row.get("position") or row.get("start")
                try:
                    pos_i = int(float(pos))
                except (TypeError, ValueError):
                    continue
                out.append({
                    "run": acc,
                    "chrom": chrom,
                    "start": pos_i,
                    "end": pos_i + 1,
                    "ref": row.get("ref") or "",
                    "alt": row.get("alt") or "",
                })
    return out


def interval_index(genes):
    by = defaultdict(list)
    for g in genes:
        by[g["chrom"]].append(g)
    for c in by:
        by[c].sort(key=lambda x: x["start"])
    return by


def overlapping_genes(idx, chrom, start, end):
    hits = []
    for g in idx.get(chrom, []):
        if g["end"] <= start:
            continue
        if g["start"] >= end:
            break
        hits.append(g)
    return hits


def tally(features, genes_idx, weight_key=None):
    counts = defaultdict(lambda: defaultdict(float))
    for feat in features:
        hits = overlapping_genes(genes_idx, feat["chrom"], feat["start"], feat["end"])
        w = float(feat.get(weight_key) or 1) if weight_key else 1.0
        if not hits:
            counts[feat["run"]][("__intergenic__", "__intergenic__")] += w
            continue
        share = w / len(hits)
        for g in hits:
            counts[feat["run"]][(g["systematic"], g["symbol"])] += share
    return counts


def write_tsv(path: Path, headers, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=headers, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    gff = ensure_gff()
    genes = parse_gff(gff)
    print(f"genes={len(genes)} from {gff}")
    idx = interval_index(genes)

    hotspots = load_hotspots()
    signals = load_signals()
    muts = load_mutations()
    print(f"hotspots={len(hotspots)} signals={len(signals)} mutations={len(muts)}")

    hot_t = tally(hotspots, idx, "total")
    sig_t = tally(signals, idx)
    mut_t = tally(muts, idx)

    all_keys = set()
    for d in (hot_t, sig_t, mut_t):
        for run, m in d.items():
            all_keys.update(m)

    runs = sorted({*hot_t, *sig_t, *mut_t})
    rows = []
    for sys, sym in sorted(all_keys, key=lambda x: (-sum(hot_t[r].get(x, 0) for r in runs), x[0])):
        rec = {"systematic": sys, "symbol": sym}
        for run in runs:
            rec[f"hotspot_{run}"] = round(hot_t[run].get((sys, sym), 0), 2)
            rec[f"signal_{run}"] = round(sig_t[run].get((sys, sym), 0), 2)
            rec[f"mut_{run}"] = round(mut_t[run].get((sys, sym), 0), 2)
        rec["hotspot_sum"] = round(sum(rec.get(f"hotspot_{r}", 0) for r in runs), 2)
        rec["signal_sum"] = round(sum(rec.get(f"signal_{r}", 0) for r in runs), 2)
        rec["mut_sum"] = round(sum(rec.get(f"mut_{r}", 0) for r in runs), 2)
        rows.append(rec)

    gene_headers = ["systematic", "symbol", "hotspot_sum", "signal_sum", "mut_sum"]
    for run in runs:
        gene_headers += [f"hotspot_{run}", f"signal_{run}", f"mut_{run}"]
    write_tsv(OUT / "gene_rnmp.tsv", gene_headers, rows)

    # Failed excision: HydEn HII vs NT (same rnh201 background)
    excision = []
    for rec in rows:
        nt = rec.get(f"hotspot_{HYDEN_NT}", 0) + rec.get(f"signal_{HYDEN_NT}", 0)
        hii = rec.get(f"hotspot_{HYDEN_HII}", 0) + rec.get(f"signal_{HYDEN_HII}", 0)
        if nt == 0 and hii == 0:
            continue
        ratio = (hii + 1) / (nt + 1)
        excision.append({
            "systematic": rec["systematic"],
            "symbol": rec["symbol"],
            "nt_score": round(nt, 2),
            "hii_score": round(hii, 2),
            "hii_over_nt": round(ratio, 3),
            "class": (
                "HII_revealed" if ratio >= 2 and hii >= 10
                else "NT_biased" if ratio <= 0.5 and nt >= 10
                else "shared"
            ),
        })
    excision.sort(key=lambda r: (-r["hii_over_nt"], -r["hii_score"]))
    write_tsv(
        OUT / "failed_excision.tsv",
        ["systematic", "symbol", "nt_score", "hii_score", "hii_over_nt", "class"],
        excision,
    )

    # Stall + mutation: genes with hotspot AND WGS rNMP-transition/mutation
    stall = []
    for rec in rows:
        yeast_hot = rec["hotspot_sum"]
        wgs_sig = rec.get(f"signal_{HU}", 0) + rec.get(f"signal_{C0}", 0)
        wgs_mut = rec.get(f"mut_{HU}", 0) + rec.get(f"mut_{C0}", 0)
        if yeast_hot <= 0 and wgs_sig <= 0 and wgs_mut <= 0:
            continue
        if yeast_hot >= 5 and (wgs_sig > 0 or wgs_mut > 0):
            stall.append({
                "systematic": rec["systematic"],
                "symbol": rec["symbol"],
                "hotspot_sum": rec["hotspot_sum"],
                "hu_signal": rec.get(f"signal_{HU}", 0),
                "c0_signal": rec.get(f"signal_{C0}", 0),
                "hu_mut": rec.get(f"mut_{HU}", 0),
                "c0_mut": rec.get(f"mut_{C0}", 0),
                "hu_over_c0_signal": round((rec.get(f"signal_{HU}", 0) + 1) / (rec.get(f"signal_{C0}", 0) + 1), 3),
            })
    stall.sort(key=lambda r: (-r["hotspot_sum"], -r["hu_signal"]))
    write_tsv(
        OUT / "stall_mutation_genes.tsv",
        ["systematic", "symbol", "hotspot_sum", "hu_signal", "c0_signal", "hu_mut", "c0_mut", "hu_over_c0_signal"],
        stall,
    )

    ftd = load_catalog(CAT / "ftd_aging_tau.tsv")
    metab = load_catalog(CAT / "one_carbon_energy.tsv")
    by_yeast = {}
    for row in ftd + metab:
        ys = (row.get("yeast_symbol") or "").strip()
        sys = (row.get("yeast_systematic") or "").strip()
        if ys in {"", "—", "-"} and sys in {"", "—", "-"}:
            continue
        key = sys or ys
        by_yeast.setdefault(key, []).append(row)

    overlay = []
    gene_by_sys = {g["systematic"]: g for g in genes}
    gene_by_sym = {g["symbol"]: g for g in genes}
    rec_by = {(r["systematic"], r["symbol"]): r for r in rows}
    rec_by_sys = {r["systematic"]: r for r in rows}
    rec_by_sym = {r["symbol"]: r for r in rows}
    for key, members in by_yeast.items():
        rec = rec_by_sys.get(key) or rec_by_sym.get(key)
        g = gene_by_sys.get(key) or gene_by_sym.get(key)
        for m in members:
            overlay.append({
                "human_symbol": m.get("human_symbol"),
                "class": m.get("class"),
                "role": m.get("role"),
                "yeast_symbol": m.get("yeast_symbol"),
                "yeast_systematic": m.get("yeast_systematic") or (g or {}).get("systematic", ""),
                "hotspot_sum": rec["hotspot_sum"] if rec else 0,
                "signal_sum": rec["signal_sum"] if rec else 0,
                "mut_sum": rec["mut_sum"] if rec else 0,
                "hii_over_nt": next((e["hii_over_nt"] for e in excision if e["systematic"] == (rec or {}).get("systematic")), ""),
                "licensed": "yeast_ortholog_only — not human FTD/tau evidence",
            })
    overlay.sort(key=lambda r: (-float(r["hotspot_sum"] or 0), r["human_symbol"] or ""))
    write_tsv(
        OUT / "ftd_tau_overlay.tsv",
        ["human_symbol", "class", "role", "yeast_symbol", "yeast_systematic",
         "hotspot_sum", "signal_sum", "mut_sum", "hii_over_nt", "licensed"],
        overlay,
    )

    top_pref = [r for r in rows if r["systematic"] != "__intergenic__"][:25]
    top_fail = [e for e in excision if e["systematic"] != "__intergenic__"][:20]
    top_stall = stall[:20]
    top_overlay = overlay[:25]

    def tbl(headers, recs, keys):
        th = "".join(f"<th>{h}</th>" for h in headers)
        body = []
        for r in recs:
            body.append("<tr>" + "".join(f"<td>{r.get(k,'')}</td>" for k in keys) + "</tr>")
        return f"<table><thead><tr>{th}</tr></thead><tbody>{''.join(body)}</tbody></table>"

    html = f"""<!DOCTYPE html><html><head><meta charset=utf-8>
<title>Gene questions — incorporation, failed excision, stall, FTD/tau overlay</title>
<style>
body{{font-family:-apple-system,Segoe UI,sans-serif;margin:1.2rem;max-width:1100px;color:#222;font-size:14px;line-height:1.45}}
table{{border-collapse:collapse;width:100%;margin:.8rem 0 1.3rem;font-size:11.5px}}
td,th{{border:1px solid #ccc;padding:.35rem .45rem;text-align:left;vertical-align:top}}
th{{background:#f4f4f4}} .warn{{background:#fff6d6;padding:.8rem 1rem;border:1px solid #e6d27a}}
@media print{{@page{{size:letter landscape;margin:12mm}} body{{margin:0;max-width:none}} thead{{display:table-header-group}} tr{{page-break-inside:avoid}}}}
</style></head><body>
<h1>Gene-level rNMP questions</h1>
<p>Research prototype. No clinical claims. Yeast <code>sacCer3</code> only for counts below.</p>
<div class="warn"><strong>Licensed reading.</strong> Preferential incorporation and failed excision
are yeast RER-null measurements (Koh ribose-seq + Zhou HydEn-seq). HydEn RNase HII vs NT
is the only clean within-assay “failed excision / remaining rNMP” contrast.
<strong>S. cerevisiae has no MAPT/tau.</strong> FTD / aging / tau-phospho rows are
ortholog-set membership, not a human disease result. HU vs C0 is n=1 vs 1 WGS
(~1.15× after depth). Human patient / 1-carbon-drug tissue is a <em>dataset hunt</em>,
not a result in this table.</div>
<h2>Preferential incorporation (top genes by hotspot mass)</h2>
{tbl(["Systematic","Symbol","Hotspot sum","Signal sum","Mut sum"], top_pref, ["systematic","symbol","hotspot_sum","signal_sum","mut_sum"])}
<h2>Failed excision proxy (HydEn RNase HII / NT)</h2>
<p>Ratio (HII+1)/(NT+1). <code>HII_revealed</code> = more unique/hotspot mass after RNase HII cleavage.</p>
{tbl(["Systematic","Symbol","NT","HII","HII/NT","Class"], top_fail, ["systematic","symbol","nt_score","hii_score","hii_over_nt","class"])}
<h2>Stall + mutation candidates (hotspot ∩ WGS signal)</h2>
{tbl(["Systematic","Symbol","Hotspot","HU sig","C0 sig","HU mut","C0 mut","HU/C0"], top_stall, ["systematic","symbol","hotspot_sum","hu_signal","c0_signal","hu_mut","c0_mut","hu_over_c0_signal"])}
<h2>FTD / aging / tau-kinase / 1-carbon overlay (yeast orthologs only)</h2>
{tbl(["Human","Class","Yeast","Hotspot","Signal","Mut","HII/NT"], top_overlay, ["human_symbol","class","yeast_symbol","hotspot_sum","signal_sum","mut_sum","hii_over_nt"])}
<p>Full tables: <code>gene_rnmp.tsv</code>, <code>failed_excision.tsv</code>,
<code>stall_mutation_genes.tsv</code>, <code>ftd_tau_overlay.tsv</code>.</p>
</body></html>
"""
    (OUT / "gene_questions.html").write_text(html)
    print(f"wrote {OUT}/ gene_rnmp={len(rows)} excision={len(excision)} stall={len(stall)} overlay={len(overlay)}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
 haplo_phylo_v11.py

 Summary
 - Input: bgzip-compressed, tabix-indexed VCF (.vcf.gz), read via cyvcf2.
 - Pivot-based derived inference:
   * Pivot-1 = site with maximum heterozygosity; infer ancestor/descendant bits against pivot-1.
   * Choose Pivot-2 = site with maximum derived count from pass-1; re-infer; pass-2 descendant bits = final derived.
   * Exclude derived==0; optionally exclude singletons (derived==1) via --exclude-singletons / --min-derived 2.
   * Sort sites by derived count (older first). Build wavelet-tree-like Newick over unique haplotypes.
 - Haplotype naming change:
   * Use the first owner in input order: pattern -> "<sample_id>_0" for phase-0 hap, "<sample_id>_1" for phase-1 hap.
   * If the same pattern appears later, reuse the earlier name (do not create a new name).
 - hapcount_by_sample.tsv uses these names and keeps **phase order** (no lexicographic sort).
 - Tree deduplication: wavelet tree is built over **unique haplotypes only** (no duplicates).
 - Optional timing: --time-report writes {out}.timings.tsv with step-wise durations.

 Comments/docstrings are in English.
"""

import argparse
import sys
from collections import Counter, OrderedDict
from dataclasses import dataclass

try:
    from cyvcf2 import VCF
except Exception:
    VCF = None

@dataclass
class SiteInfo:
    chrom: str
    pos: int
    ref: str
    alt: str  # comma-joined ALT alleles
    id: str = ""

# ------------------------------
# Utils
# ------------------------------

def parse_region(region: str):
    if region is None:
        return None
    chrom, span = region.split(":")
    start, end = map(int, span.replace(",", "").split("-"))
    return chrom, start, end

# ------------------------------
# VCF.gz loader
# ------------------------------

def load_vcfgz(path: str, region: str = None):
    if VCF is None:
        raise RuntimeError("cyvcf2 is required. Install it and ensure the VCF has a .tbi index.")
    vcf = VCF(path)
    samples = list(vcf.samples)
    iterator = vcf(region) if region else vcf
    records = []
    for var in iterator:
        chrom = var.CHROM
        pos = var.POS
        vid = var.ID or ""
        ref = var.REF
        alt = ",".join(var.ALT) if isinstance(var.ALT, (list, tuple)) else (var.ALT or "")
        gts = []
        for gt in var.genotypes:
            if len(gt) >= 3:
                a1, a2, phased_flag = gt[0], gt[1], bool(gt[2])
            elif len(gt) == 2:
                a1, a2 = gt[0], gt[1]
                phased_flag = False
            else:
                a1, a2, phased_flag = -1, -1, False
            gts.append((a1, a2, phased_flag))
        records.append({"chrom": chrom, "pos": pos, "id": vid, "ref": ref, "alt": alt, "gts": gts})
    return samples, records

# ------------------------------
# Sample filtering
# ------------------------------

def read_fam_ids(path: str):
    ids = set()
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            cols = line.split()
            if len(cols) >= 2:
                ids.add(cols[1])  # IID
    return ids


def read_simple_id_list(path: str):
    ids = set()
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            s = line.strip()
            if s:
                ids.add(s)
    return ids


def filter_samples(samples, records, include_ids=None, exclude_ids=None):
    include_ids = include_ids or set()
    exclude_ids = exclude_ids or set()
    keep_indices = []
    for i, sid in enumerate(samples):
        if include_ids and sid not in include_ids:
            continue
        if exclude_ids and sid in exclude_ids:
            continue
        keep_indices.append(i)
    new_samples = [samples[i] for i in keep_indices]
    new_records = []
    for rec in records:
        gts_sliced = [rec["gts"][i] for i in keep_indices]
        new_rec = dict(rec)
        new_rec["gts"] = gts_sliced
        new_records.append(new_rec)
    return new_samples, new_records

# ------------------------------
# Heterozygosity & haplotype reconstruction
# ------------------------------

def compute_site_heterozygosity(records):
    het_counts = []
    for rec in records:
        cnt = 0
        for (a1, a2, _phased) in rec["gts"]:
            if a1 == -1 or a2 == -1:
                continue
            if a1 != a2:
                cnt += 1
        het_counts.append(cnt)
    return het_counts


def reconstruct_two_haps_all_sites(samples, records):
    """Build per-sample two haplotype bit strings across all sites in original order."""
    per_sample = {}
    n_sites = len(records)
    if n_sites == 0:
        for sid in samples:
            per_sample[sid] = ('NA','NA')
        return per_sample
    for s_idx, sid in enumerate(samples):
        h1 = []
        h2 = []
        phased_ok = True
        for rec in records:
            a1, a2, phased_flag = rec["gts"][s_idx]
            if not phased_flag:
                phased_ok = False
                break
            if a1 == -1 or a2 == -1:
                phased_ok = False
                break
            h1.append('1' if a1 >= 1 else '0')
            h2.append('1' if a2 >= 1 else '0')
        if not phased_ok:
            per_sample[sid] = ('NA','NA')
        else:
            per_sample[sid] = ("".join(h1), "".join(h2))
    return per_sample

# ------------------------------
# Pivot-based ancestor/descendant inference
# ------------------------------

def majority_bit_at_site(site_idx, per_sample_patterns):
    zero = 0
    one = 0
    for p1, p2 in per_sample_patterns.values():
        for pat in (p1, p2):
            if pat == 'NA':
                continue
            if site_idx >= len(pat):
                continue
            b = pat[site_idx]
            if b == '1':
                one += 1
            elif b == '0':
                zero += 1
    return '1' if one >= zero else '0'


def infer_desc_bits(pivot_idx, pivot_desc_bit, per_sample_patterns, n_sites):
    desc = []
    for j in range(n_sites):
        eq = 0
        tot = 0
        for p1, p2 in per_sample_patterns.values():
            for pat in (p1, p2):
                if pat == 'NA':
                    continue
                if pivot_idx >= len(pat) or j >= len(pat):
                    continue
                pb = pat[pivot_idx]
                sb = pat[j]
                if pb == pivot_desc_bit:
                    tot += 1
                    if sb == pivot_desc_bit:
                        eq += 1
        if tot == 0:
            desc.append(pivot_desc_bit)
        else:
            if eq >= (tot - eq):
                desc.append(pivot_desc_bit)
            else:
                desc.append('0' if pivot_desc_bit == '1' else '1')
    return desc


def count_derived_from_desc_bits(desc_bits, per_sample_patterns):
    n_sites = len(desc_bits)
    counts = [0]*n_sites
    for j in range(n_sites):
        d = desc_bits[j]
        c = 0
        for p1, p2 in per_sample_patterns.values():
            for pat in (p1, p2):
                if pat == 'NA':
                    continue
                if j >= len(pat):
                    continue
                if pat[j] == d:
                    c += 1
        counts[j] = c
    return counts

# ------------------------------
# Variant IDs & window
# ------------------------------

def build_variant_ids_and_window(varinfo_list):
    if not varinfo_list:
        return "", ""
    variant_ids = []
    chrom = varinfo_list[0].chrom
    min_pos = varinfo_list[0].pos
    max_pos = varinfo_list[0].pos
    for v in varinfo_list:
        vid = v.id if v.id else f"{v.chrom}:{v.pos}:{v.ref}:{v.alt}"
        variant_ids.append(vid)
        min_pos = min(min_pos, v.pos)
        max_pos = max(max_pos, v.pos)
    window = f"{chrom}:{min_pos}-{max_pos}"
    return ";".join(variant_ids), window

# ------------------------------
# Haplotype naming (first owner in input order)
# ------------------------------

def assign_hap_names(per_sample_patterns, samples):
    """
    Map unique allele patterns to names based on the first owner (input order).
    Names: "<sample_id>_0" for phase-0, "<sample_id>_1" for phase-1.
    Later occurrences of the same pattern reuse the earlier name.
    Returns:
      pattern_to_name: OrderedDict (in insertion order) mapping pattern -> name
      counts: Counter(patterns)
      per_sample_names: dict[sample_id] -> (name_for_phase0, name_for_phase1) or ('NA','NA')
    """
    pattern_to_name = OrderedDict()
    counts = Counter()
    per_sample_names = {}

    # Build names in input order
    for sid in samples:
        p1, p2 = per_sample_patterns.get(sid, ('NA','NA'))
        if p1 == 'NA' or p2 == 'NA':
            per_sample_names[sid] = ('NA','NA')
            continue
        # phase-0
        if p1 not in pattern_to_name:
            pattern_to_name[p1] = f"{sid}_0"
        # phase-1
        if p2 not in pattern_to_name:
            pattern_to_name[p2] = f"{sid}_1"
        # counts
        counts.update([p1, p2])
        # per-sample name tuple preserves phase order
        per_sample_names[sid] = (pattern_to_name[p1], pattern_to_name[p2])

    # For samples after initial assignment (if any NA were resolved later), ensure names resolved
    for sid in samples:
        p1, p2 = per_sample_patterns.get(sid, ('NA','NA'))
        if p1 == 'NA' or p2 == 'NA':
            continue
        # counts were updated above; ensure per_sample_names present
        if sid not in per_sample_names:
            per_sample_names[sid] = (pattern_to_name[p1], pattern_to_name[p2])

    return pattern_to_name, counts, per_sample_names

# ------------------------------
# Outputs
# ------------------------------

def write_hapfreq_tsv(out_prefix, pattern_to_name, counts, variant_ids, window, missing_policy):
    path = f"{out_prefix}.hapfreq.tsv"
    total_haps = sum(counts.values()) if counts else 0
    with open(path, 'w', encoding='utf-8') as fw:
        fw.write("hap_id\tcount\tfreq\tallele_pattern\tvariant_ids\twindow\tmissing_policy\n")
        # pattern_to_name preserves first-owner insertion order
        for pat, name in pattern_to_name.items():
            c = counts.get(pat, 0)
            freq = (c/total_haps) if total_haps>0 else 0.0
            fw.write(f"{name}\t{c}\t{freq:.6f}\t{pat}\t{variant_ids}\t{window}\t{missing_policy}\n")


def write_hapcount_by_sample_tsv(out_prefix, per_sample_names):
    path = f"{out_prefix}.hapcount_by_sample.tsv"
    with open(path, 'w', encoding='utf-8') as fw:
        fw.write("sample_id\thap1_id\thap2_id\n")
        for sid in per_sample_names:
            h1, h2 = per_sample_names[sid]
            fw.write(f"{sid}\t{h1}\t{h2}\n")

# ------------------------------
# Wavelet tree over unique haplotypes
# ------------------------------

def build_wavelet_tree_newick(pattern_to_name, site_order_len):
    """
    Build a wavelet-tree-like Newick over **unique** haplotypes.
    Leaves are haplotype names (first-owner), each corresponding to a unique pattern.
    """
    name_to_pat = {name: pat for pat, name in pattern_to_name.items()}
    leaves = list(name_to_pat.keys())
    if len(leaves)==0:
        return "();"
    if len(leaves)==1:
        return f"({leaves[0]}:0.0);"
    total = site_order_len if site_order_len>0 else 1
    def rec(names, depth):
        if len(names)==1:
            return names[0], 0.0
        if depth>=total:
            return "("+",".join(names)+")", 0.0
        zeros, ones = [], []
        for nm in names:
            pat = name_to_pat[nm]
            bit = pat[depth] if depth < len(pat) else '0'
            (zeros if bit=='0' else ones).append(nm)
        if len(zeros)==0 or len(ones)==0:
            return rec(names, depth+1)
        lstr, lh = rec(zeros, depth+1)
        rstr, rh = rec(ones, depth+1)
        bl = 1.0/total
        return f"({lstr}:{max(bl-lh,0):.6f},{rstr}:{max(bl-rh,0):.6f})", bl
    newick, _h = rec(leaves, 0)
    return newick+";"

# ------------------------------
# Filtered sites report
# ------------------------------

def write_filtered_sites_report(out_prefix, varinfo, derived_counts_all, excluded_idx, min_derived):
    path = f"{out_prefix}.filtered_sites.tsv"
    with open(path, 'w', encoding='utf-8') as fw:
        fw.write("site_index\tchrom\tpos\tref\talt\tderived_count\treason\n")
        for i in excluded_idx:
            si = varinfo[i]
            c = derived_counts_all[i]
            if c == 0:
                reason = "derived==0 (uninformative)"
            elif c < max(1, min_derived):
                reason = f"derived<{min_derived}"
            else:
                reason = "excluded"
            fw.write(f"{i}\t{si.chrom}\t{si.pos}\t{si.ref}\t{si.alt}\t{c}\t{reason}\n")

# ------------------------------
# Main
# ------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="infile", required=True, help="Input bgzip-compressed VCF (.vcf.gz) with tabix index (.tbi)")
    ap.add_argument("--region", default=None, help="Region in 'chr:start-end' format (optional)")
    ap.add_argument("--out-prefix", required=True, help="Output prefix (e.g., 'frequency')")
    ap.add_argument("--exclude-singletons", action="store_true", help="Exclude sites with derived count == 1 (singleton). Equivalent to --min-derived 2")
    ap.add_argument("--min-derived", type=int, default=None, help="Minimum derived allele count to retain (>=1 enforced). Use 2 to exclude singletons.")
    ap.add_argument("--report-filtered-sites", action="store_true", default=True, help="Write excluded sites to {out}.filtered_sites.tsv")
    ap.add_argument("--require-phased", action="store_true", help="If set, any unphased sample is assigned NA/NA in hapcount_by_sample.tsv")
    ap.add_argument("--missing-policy", default="ignore", choices=["ignore","code-3","pairwise"], help="Label recorded in hapfreq.tsv")
    # Sample filtering options
    ap.add_argument("--fam-include", default=None, help="PLINK .fam file; keep only IIDs listed")
    ap.add_argument("--fam-exclude", default=None, help="PLINK .fam file; exclude IIDs listed")
    ap.add_argument("--samples-include", default=None, help="Text file of sample IDs to keep (one per line)")
    ap.add_argument("--samples-exclude", default=None, help="Text file of sample IDs to drop (one per line)")
    # Timing option
    ap.add_argument("--time-report", action="store_true", help="Write step-wise timings to {out}.timings.tsv")

    args = ap.parse_args()

    if VCF is None:
        raise RuntimeError("cyvcf2 is not available. Please install it (pip install cyvcf2).")

    timings = {}
    t0 = time.perf_counter()

    # 1) Load VCF
    samples, records = load_vcfgz(args.infile, region=args.region)
    if len(records)==0:
        print("No variants in the specified input/region.", file=sys.stderr)
        return
    timings['load_vcf'] = time.perf_counter() - t0

    # 2) Sample filtering
    t1 = time.perf_counter()
    include_ids = set()
    exclude_ids = set()
    if args.fam_include:
        include_ids |= read_fam_ids(args.fam_include)
    if args.fam_exclude:
        exclude_ids |= read_fam_ids(args.fam_exclude)
    if args.samples_include:
        include_ids |= read_simple_id_list(args.samples_include)
    if args.samples_exclude:
        exclude_ids |= read_simple_id_list(args.samples_exclude)
    if include_ids or exclude_ids:
        samples, records = filter_samples(samples, records, include_ids=(include_ids or None), exclude_ids=(exclude_ids or None))
        if len(samples)==0:
            print("After filtering, no samples remain.", file=sys.stderr)
            return
    timings['filter_samples'] = time.perf_counter() - t1

    # 3) Build haplotypes across ALL sites (original order)
    t2 = time.perf_counter()
    per_sample_patterns_all = reconstruct_two_haps_all_sites(samples, records)
    timings['reconstruct_all_sites'] = time.perf_counter() - t2

    # 4) Heterozygosity -> pivot-1
    t3 = time.perf_counter()
    het_counts = compute_site_heterozygosity(records)
    pivot1_idx = max(range(len(records)), key=lambda i: (het_counts[i], -i))
    pivot1_desc_bit = majority_bit_at_site(pivot1_idx, per_sample_patterns_all)
    timings['pivot1_selection'] = time.perf_counter() - t3

    # 5) Pass-1 inference
    t4 = time.perf_counter()
    n_sites = len(records)
    desc_bits_p1 = infer_desc_bits(pivot1_idx, pivot1_desc_bit, per_sample_patterns_all, n_sites)
    derived_counts_p1 = count_derived_from_desc_bits(desc_bits_p1, per_sample_patterns_all)
    timings['pass1_inference'] = time.perf_counter() - t4

    # 6) Pivot-2 -> Pass-2
    t5 = time.perf_counter()
    pivot2_idx = max(range(n_sites), key=lambda i: (derived_counts_p1[i], -i))
    pivot2_desc_bit = desc_bits_p1[pivot2_idx]
    desc_bits_final = infer_desc_bits(pivot2_idx, pivot2_desc_bit, per_sample_patterns_all, n_sites)
    derived_counts_final = count_derived_from_desc_bits(desc_bits_final, per_sample_patterns_all)
    timings['pass2_inference'] = time.perf_counter() - t5

    # 7) Filtering & ordering
    t6 = time.perf_counter()
    min_derived = 2 if args.exclude_singletons else 1
    if args.min_derived is not None:
        min_derived = max(1, args.min_derived)
    keep_idx = [i for i, c in enumerate(derived_counts_final) if c >= max(1, min_derived)]
    excluded_idx = [i for i in range(n_sites) if i not in keep_idx]

    if args.report_filtered_sites:
        varinfo_all = [SiteInfo(chrom=r["chrom"], pos=r["pos"], ref=r["ref"], alt=r["alt"], id=r["id"]) for r in records]
        write_filtered_sites_report(args.out_prefix, varinfo_all, derived_counts_final, excluded_idx, min_derived)

    derived_kept = [derived_counts_final[i] for i in keep_idx]
    order_local = sorted(range(len(keep_idx)), key=lambda k: derived_kept[k])
    kept_order = [keep_idx[k] for k in order_local]
    timings['filter_sort_sites'] = time.perf_counter() - t6

    # 8) Build variant_ids & window
    t7 = time.perf_counter()
    varinfo_all = [SiteInfo(chrom=r["chrom"], pos=r["pos"], ref=r["ref"], alt=r["alt"], id=r["id"]) for r in records]
    varinfo_s = [varinfo_all[i] for i in kept_order]
    variant_ids, window = build_variant_ids_and_window(varinfo_s)
    timings['variant_window'] = time.perf_counter() - t7

    # 9) Reconstruct per-sample haps on kept+sorted order
    t8 = time.perf_counter()
    per_sample_patterns = {}
    for sid, (p1, p2) in per_sample_patterns_all.items():
        if p1=='NA' or p2=='NA':
            per_sample_patterns[sid] = ('NA','NA')
        else:
            try:
                pat1 = ''.join(p1[i] for i in kept_order)
                pat2 = ''.join(p2[i] for i in kept_order)
                per_sample_patterns[sid] = (pat1, pat2)
            except Exception:
                per_sample_patterns[sid] = ('NA','NA')
    timings['reconstruct_kept_sorted'] = time.perf_counter() - t8

    # 10) Assign names based on first owner
    t9 = time.perf_counter()
    pattern_to_name, counts, per_sample_names = assign_hap_names(per_sample_patterns, samples)
    timings['assign_names'] = time.perf_counter() - t9

    # 11) Write outputs
    t10 = time.perf_counter()
    write_hapfreq_tsv(args.out_prefix, pattern_to_name, counts, variant_ids, window, args.missing_policy)
    write_hapcount_by_sample_tsv(args.out_prefix, per_sample_names)
    timings['write_outputs'] = time.perf_counter() - t10

    # 12) Build unique-haplotype wavelet tree
    t11 = time.perf_counter()
    newick = build_wavelet_tree_newick(pattern_to_name, site_order_len=len(kept_order))
    with open(f"{args.out_prefix}.newick", 'w', encoding='utf-8') as fw:
        fw.write(newick+"\n")
    timings['build_tree'] = time.perf_counter() - t11

    # 13) Timings report
    if args.time_report:
        path = f"{args.out_prefix}.timings.tsv"
        with open(path, 'w', encoding='utf-8') as fw:
            fw.write("step\tseconds\n")
            for k in ['load_vcf','filter_samples','reconstruct_all_sites','pivot1_selection','pass1_inference','pass2_inference','filter_sort_sites','variant_window','reconstruct_kept_sorted','assign_names','write_outputs','build_tree']:
                if k in timings:
                    fw.write(f"{k}\t{timings[k]:.6f}\n")

if __name__ == "__main__":
    import time
    main()

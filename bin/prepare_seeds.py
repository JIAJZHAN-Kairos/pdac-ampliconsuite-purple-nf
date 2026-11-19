#!/usr/bin/env python3
"""Convert 1-based inclusive PURPLE segments to 0-based half-open seed BEDs."""
import csv
import sys
from collections import defaultdict
from decimal import Decimal
from pathlib import Path
from zipfile import ZipFile

CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
MAX_FOCAL_BP = 30_000_000


def bed_mask(path):
    intervals = defaultdict(list)
    for line in Path(path).read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            chrom, start, end, *_ = line.split()
            intervals[chrom].append((int(start), int(end)))
    return {chrom: sorted(blocks) for chrom, blocks in intervals.items()}


def unmasked(start, end, blocks):
    for left, right in blocks:
        if left >= end:
            break
        if right > start:
            if left > start:
                yield start, left
            start = right
    if start < end:
        yield start, end


def prepare(sample, cnv, fit, mask):
    ploidy = Decimal(fit["ploidy"])
    purity = Decimal(fit["purity"])
    if not ploidy.is_finite() or ploidy <= 0 or not 0 < purity <= 1:
        raise ValueError(f"{sample}: invalid PURPLE fit")
    calls = []
    seeds = []
    candidates = {"gain": [], "amplification": []}
    for row in cnv:
        chrom = row["chromosome"]
        if not chrom.startswith("chr"):
            chrom = "chr" + chrom
        if chrom not in CHROMS:
            continue
        start, end = int(row["start"]) - 1, int(row["end"])
        cn = Decimal(row["copyNumber"])
        if not (0 <= start < end and cn.is_finite()):
            raise ValueError(f"{sample}: invalid PURPLE segment {row}")
        if cn >= 3 * ploidy:
            candidates["amplification"].append((chrom, start, end, cn))
        elif cn >= Decimal("1.5") * ploidy:
            candidates["gain"].append((chrom, start, end, cn))
    # Gains and amplifications are grouped separately so that gain shoulders
    # do not make an adjacent focal amplification count as broad.
    for call, segments in candidates.items():
        segments.sort(key=lambda row: (CHROMS.index(row[0]), row[1]))
        groups = []
        for segment in segments:
            if groups and groups[-1][-1][0] == segment[0] and groups[-1][-1][2] == segment[1]:
                groups[-1].append(segment)
            else:
                groups.append([segment])
        for group in groups:
            broad = group[-1][2] - group[0][1] > MAX_FOCAL_BP
            for chrom, start, end, cn in group:
                pieces = [] if broad else list(unmasked(start, end, mask.get(chrom, [])))
                length = sum(right - left for left, right in pieces)
                reason = ("broad" if broad else "masked" if not pieces else
                          "partially_masked" if length < end - start else "seeded")
                calls.append((chrom, start, end, cn, cn / ploidy, call, reason))
                seeds.extend((chrom, left, right, cn) for left, right in pieces)
    calls.sort(key=lambda row: (CHROMS.index(row[0]), row[1]))
    seeds.sort(key=lambda row: (CHROMS.index(row[0]), row[1]))
    return calls, seeds


def tsv(lines):
    return csv.DictReader((line.decode() for line in lines), delimiter="\t")


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as output:
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def main():
    manifest, bundle, mask_path = sys.argv[1:]
    mask = bed_mask(mask_path)
    Path("ready").mkdir()
    Path("seed_audit").mkdir()
    summary = []
    with ZipFile(bundle) as archive, open(manifest, newline="") as sheet:
        for row in csv.DictReader(sheet):
            sample = row["sample"]
            with archive.open(f"{sample}.purple.purity.tsv") as fit_file:
                fit, = tsv(fit_file)
            with archive.open(f"{sample}.purple.cnv.somatic.tsv") as cnv_file:
                calls, seeds = prepare(sample, tsv(cnv_file), fit, mask)
            key = f"{sample}_AA_CNV_SEEDS.bed"
            seed_text = "".join("\t".join(map(str, line)) + "\n" for line in seeds)
            Path("seed_audit", key).write_text(seed_text)
            write_tsv(Path("seed_audit", f"{sample}_cn_calls.tsv"),
                      ["chrom", "start0", "end", "tumor_cn", "cn_over_ploidy",
                       "call", "seed_status"], calls)
            if seeds:
                Path("ready", key).write_text(seed_text)
            summary.append((sample, fit["purity"], fit["ploidy"], fit.get("status", ""), row.get("purple_qc", ""),
                            sum(call[5] == "gain" for call in calls),
                            sum(call[5] == "amplification" for call in calls), len(seeds),
                            "seeded" if seeds else "no_eligible_seed",
                            row["kras_qc"], row["kras_vaf"], row["kras_depth"], row["kras_variant"]))
    write_tsv("seed_audit/samples.tsv",
              ["sample", "purity", "ploidy", "purple_fit_status", "purple_qc", "gain_segments",
               "amplification_segments", "seed_intervals", "seed_status",
               "kras_qc", "kras_vaf", "kras_depth", "kras_variant"], summary)


if __name__ == "__main__":
    main()

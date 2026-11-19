#!/usr/bin/env python3
"""Verify AA/AC output and report the ecContext label for each ecDNA feature."""
import csv
import re
import sys
from pathlib import Path

ORIGIN = {
    "Simple circular simple background": "simple_excision_like",
    "Two-foldback": "simple_excision_like",
    "Simple circular complex background": "chromothripsis_associated",
    "BFB-like": "chromothripsis_associated",
    "Heavily rearranged unichromosomal": "chromothripsis_associated",
    "Heavily rearranged multichromosomal": "chromothripsis_associated",
    "Unknown": "unresolved",
}


def main():
    sample = sys.argv[1]
    results = Path("results")
    aa = results / f"{sample}_AA_results"
    ac = results / f"{sample}_classification"
    amplicons = int(re.search(r"^#Amplicons = (\d+)$",
                              (aa / f"{sample}_summary.txt").read_text(), re.M)[1])
    graphs = {path.name.removesuffix("_graph.txt") for path in aa.glob("*_amplicon*_graph.txt")}
    cycles = {path.name.removesuffix("_cycles.txt") for path in aa.glob("*_amplicon*_cycles.txt")
              if not path.name.endswith("_annotated_cycles.txt")}
    if graphs != cycles or len(graphs) != amplicons:
        raise ValueError(f"{sample}: incomplete AA graph/cycle output ({amplicons})")
    context = ac / f"{sample}_ecDNA_context_calls.tsv"
    counts = ac / f"{sample}_ecDNA_counts.tsv"
    calls = []
    if amplicons:
        if not (ac / f"{sample}_amplicon_classification_profiles.tsv").is_file():
            raise FileNotFoundError(f"{sample}: AC classification profile missing")
        with context.open() as handle:
            calls = list(csv.reader(handle, delimiter="\t"))
        with counts.open() as handle:
            raw_counts = dict(csv.reader(handle, delimiter="\t"))
        if len(calls) != int(raw_counts[sample]) or len({row[0] for row in calls}) != len(calls):
            raise ValueError(f"{sample}: ecContext and ecDNA counts disagree")
    with (results / f"{sample}_origin.tsv").open("w") as handle:
        handle.write("sample\tecDNA_feature\tecContext\torigin_group\n")
        for feature, label in calls:
            handle.write(f"{sample}\t{feature}\t{label}\t{ORIGIN.get(label, 'unresolved')}\n")
    with (results / f"{sample}_status.tsv").open("w") as handle:
        handle.write("sample\tstatus\taa_amplicons\tecdna_features\n")
        handle.write(f"{sample}\t{'classified' if amplicons else 'no_AA_amplicons'}\t"
                     f"{amplicons}\t{len(calls)}\n")


if __name__ == "__main__":
    main()

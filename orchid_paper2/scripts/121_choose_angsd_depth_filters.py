#!/usr/bin/env python3
"""Choose group-specific upper depth filters from ANGSD depth histograms."""
from __future__ import annotations
import argparse
from pathlib import Path


def weighted_quantile(counts, probability):
    target = probability * sum(counts)
    cumulative = 0
    for depth, count in enumerate(counts):
        cumulative += count
        if cumulative >= target:
            return depth
    return len(counts) - 1


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--work", required=True, type=Path)
    parser.add_argument("--quantile", type=float, default=0.995)
    args = parser.parse_args()
    output = args.work / "inputs" / "depth_filters.tsv"
    with output.open("w") as handle:
        handle.write("group\tmin_individuals\tretained_depth_profile_sites\tmedian_total_depth\tq995_total_depth\tmax_total_depth_filter\n")
        for group, minimum in (("PAM", 62), ("PAP", 55), ("PSA", 5)):
            path = args.work / "depth" / f"{group}.depthGlobal"
            counts = [int(float(x)) for x in path.read_text().split()]
            total = sum(counts)
            median = weighted_quantile(counts, 0.5)
            q995 = weighted_quantile(counts, args.quantile)
            # Never choose a threshold lower than two reads per required
            # individual; this guard matters only for pathological profiles.
            maximum = max(q995, 2 * minimum)
            handle.write(f"{group}\t{minimum}\t{total}\t{median}\t{q995}\t{maximum}\n")


if __name__ == "__main__":
    main()


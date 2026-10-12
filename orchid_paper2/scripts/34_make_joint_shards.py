#!/usr/bin/env python3
"""Split the callable BED into ordered, non-overlapping, whole-interval shards."""

from pathlib import Path
import argparse


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bed", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--n-shards", type=int, default=24)
    args = ap.parse_args()

    rows = []
    for line in args.bed.read_text().splitlines():
        if not line.strip():
            continue
        f = line.split("\t")
        if len(f) < 3:
            raise ValueError(f"Malformed BED line: {line}")
        start, end = int(f[1]), int(f[2])
        if end <= start:
            raise ValueError(f"Non-positive BED interval: {line}")
        rows.append((line, end - start))

    total = sum(n for _, n in rows)
    if not rows or args.n_shards < 1 or args.n_shards > len(rows):
        raise ValueError("Invalid number of records or shards")

    args.outdir.mkdir(parents=True, exist_ok=True)
    shards = [[] for _ in range(args.n_shards)]
    shard_bp = [0] * args.n_shards
    shard = 0
    consumed = 0
    for i, (line, length) in enumerate(rows):
        remaining_rows = len(rows) - i
        remaining_shards = args.n_shards - shard
        target = (total - consumed) / remaining_shards
        if (shard < args.n_shards - 1 and shards[shard]
                and shard_bp[shard] + length > target
                and remaining_rows >= remaining_shards):
            consumed += shard_bp[shard]
            shard += 1
        shards[shard].append(line)
        shard_bp[shard] += length

    summary = ["shard\tintervals\tcallable_bp\tbed"]
    for i, lines in enumerate(shards, 1):
        path = args.outdir / f"shard_{i:03d}.bed"
        path.write_text("\n".join(lines) + "\n")
        summary.append(f"{i}\t{len(lines)}\t{shard_bp[i-1]}\t{path}")
    (args.outdir / "shard_summary.tsv").write_text("\n".join(summary) + "\n")

    reconstructed = [line for lines in shards for line in lines]
    if reconstructed != [line for line, _ in rows]:
        raise RuntimeError("Shard concatenation does not reproduce input BED")
    if sum(shard_bp) != total:
        raise RuntimeError("Shard callable-base total differs from input BED")
    print(f"Wrote {len(shards)} shards: {len(rows)} intervals, {total} callable bp")


if __name__ == "__main__":
    main()

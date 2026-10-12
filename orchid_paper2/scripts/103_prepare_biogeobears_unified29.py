#!/usr/bin/env python3
from pathlib import Path
import csv

root = Path(__file__).resolve().parents[2]
out = root / "orchid_paper2/results/biogeobears_unified29/inputs"
out.mkdir(parents=True, exist_ok=True)

assign = root / "orchid_paper2/17_herodotools/inputs/bioregions/tip_to_thesis_bioregion.tsv"
mcc = root / "orchid_paper2/16_rad_beast/final_60m_summary/zhang_final_mcc_29tip_herodotools.newick"
post = root / "orchid_paper2/16_rad_beast/final_60m_summary/zhang_100_posterior_29tip.newick"

regions = [
    ("A", "taiwan", "Taiwan"),
    ("B", "north_luzon", "Northern Luzon-Batanes-Babuyan"),
    ("C", "south_luzon", "Southern Luzon-Bicol"),
    ("D", "greater_negros_panay", "Greater Negros-Panay"),
    ("E", "mindoro", "Mindoro"),
    ("F", "palawan", "Palawan"),
    ("G", "mindanao", "Mindanao"),
    ("H", "borneo", "Borneo"),
    ("I", "java", "Java"),
    ("J", "sulawesi_moluccas_lesser_sunda", "Sulawesi-Moluccas-Lesser Sunda Islands"),
    ("K", "new_guinea_queensland", "New Guinea-Queensland"),
]
idx = {rid: i for i, (_, rid, _) in enumerate(regions)}
rows = list(csv.DictReader(assign.open(), delimiter="\t"))
unknown = sorted({r["region_id"] for r in rows} - set(idx))
if unknown:
    raise SystemExit(f"Unrecognized regions: {unknown}")

with (out / "geography.data").open("w") as h:
    h.write(f"{len(rows)} {len(regions)} ({' '.join(x[0] for x in regions)})\n")
    for r in rows:
        bits = ["0"] * len(regions)
        bits[idx[r["region_id"]]] = "1"
        h.write(f"{r['tip']:<12}{''.join(bits)}\n")

with (out / "area_order.tsv").open("w") as h:
    h.write("code\tregion_id\tregion_name\n")
    for code, rid, name in regions:
        h.write(f"{code}\t{rid}\t{name}\n")
(out / "max_range_size.txt").write_text("4\n")
(out / "mcc_myr.newick").write_text(mcc.read_text().strip() + "\n")

trees = [x.strip() for x in post.read_text().splitlines() if x.strip()]
if len(trees) != 100:
    raise SystemExit(f"Expected 100 posterior trees, found {len(trees)}")
for i, tree in enumerate(trees, 1):
    (out / f"posterior_{i:03d}_myr.newick").write_text(tree + "\n")
print(f"Prepared {len(rows)} tips, {len(regions)} regions, MCC + {len(trees)} posterior trees")

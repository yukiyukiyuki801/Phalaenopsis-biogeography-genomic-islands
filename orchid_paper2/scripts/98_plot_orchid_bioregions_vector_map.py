#!/usr/bin/env python3
"""Plot the orchid bioregions directly from vector polygons, not rasters."""

from pathlib import Path
import csv

from osgeo import ogr
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch, Patch


ROOT = Path(__file__).resolve().parents[1]
GPKG = ROOT / "19_biogeobears/bioregions/biogeographic_regions.gpkg"
OUT = ROOT / "19_biogeobears/bioregion_plots"
OUT.mkdir(parents=True, exist_ok=True)

# Codes and colors match the single-region BioGeoBEARS legend. Region D has no
# sampled population tip and is intentionally drawn without a colored fill.
REGIONS = {
    "taiwan": ("A", "Taiwan", "#E41A1C", True),
    "north_luzon": ("B", "Northern Luzon-Batanes-Babuyan", "#FF7F00", True),
    "south_luzon": ("C", "Southern Luzon-Bicol", "#F0E442", True),
    "greater_negros_panay": ("D", "Greater Negros-Panay (unsampled)", "#FFFFFF", False),
    "mindoro": ("E", "Mindoro", "#4DAF4A", True),
    "palawan": ("F", "Palawan", "#00BFC4", True),
    "mindanao": ("G", "Mindanao", "#56B4E9", True),
    "borneo": ("H", "Borneo", "#0057D9", True),
    "java": ("I", "Java", "#6A00FF", True),
    "sulawesi_moluccas_lesser_sunda": (
        "J", "Sulawesi-Moluccas-Lesser Sunda Islands", "#F000C8", True
    ),
    "new_guinea_queensland": ("K", "New Guinea-Queensland", "#FF3B9D", True),
}

LABEL_XY = {
    "taiwan": (121.0, 23.4),
    "north_luzon": (120.8, 18.1),
    "south_luzon": (123.3, 13.9),
    "greater_negros_panay": (123.0, 10.4),
    "mindoro": (120.6, 12.5),
    "palawan": (118.2, 9.4),
    "mindanao": (125.4, 7.4),
    "borneo": (114.7, 0.7),
    "java": (110.4, -7.3),
    "sulawesi_moluccas_lesser_sunda": (122.0, -2.3),
    "new_guinea_queensland": (140.8, -5.0),
}


def polygon_path(poly):
    """Convert an OGR Polygon, including interior rings, to a matplotlib Path."""
    vertices, codes = [], []
    for ring_i in range(poly.GetGeometryCount()):
        ring = poly.GetGeometryRef(ring_i)
        points = [(ring.GetX(i), ring.GetY(i)) for i in range(ring.GetPointCount())]
        if not points:
            continue
        if points[0] != points[-1]:
            points.append(points[0])
        vertices.extend(points)
        codes.extend([MplPath.MOVETO] + [MplPath.LINETO] * (len(points) - 2) + [MplPath.CLOSEPOLY])
    return MplPath(vertices, codes)


def add_geometry(ax, geom, facecolor):
    gtype = ogr.GT_Flatten(geom.GetGeometryType())
    if gtype == ogr.wkbPolygon:
        ax.add_patch(PathPatch(polygon_path(geom), facecolor=facecolor,
                               edgecolor="#333333", linewidth=0.42, zorder=2))
    elif gtype == ogr.wkbMultiPolygon:
        for i in range(geom.GetGeometryCount()):
            add_geometry(ax, geom.GetGeometryRef(i), facecolor)
    else:
        raise ValueError(f"Unexpected geometry type: {geom.GetGeometryName()}")


ds = ogr.Open(str(GPKG))
if ds is None:
    raise FileNotFoundError(GPKG)
layer = ds.GetLayerByName("bioregions")
seen = set()

fig, ax = plt.subplots(figsize=(12.2, 8.2))
fig.subplots_adjust(left=0.075, right=0.69, bottom=0.10, top=0.91)
ax.set_facecolor("#EAF4FA")

for feature in layer:
    region_id = feature.GetField("region_id")
    if region_id not in REGIONS:
        raise KeyError(f"No code/color specification for {region_id}")
    code, name, color, sampled = REGIONS[region_id]
    add_geometry(ax, feature.GetGeometryRef(), color if sampled else "#FFFFFF")
    seen.add(region_id)

missing = set(REGIONS) - seen
if missing:
    raise RuntimeError(f"Regions missing from GeoPackage: {sorted(missing)}")

# Region codes are printed on representative islands; the complete code-name
# mapping remains in the integrated legend.
for region_id, (x, y) in LABEL_XY.items():
    code, _, _, sampled = REGIONS[region_id]
    ax.text(x, y, code, ha="center", va="center", fontsize=9.2, fontweight="bold",
            color="#111111", zorder=5,
            bbox=dict(boxstyle="circle,pad=0.16", facecolor="white",
                      edgecolor="#333333", linewidth=0.55, alpha=0.90))

ax.set_xlim(104, 158)
ax.set_ylim(-21, 27)
ax.set_aspect("equal", adjustable="box")
ax.set_xlabel("Longitude", fontsize=11)
ax.set_ylabel("Latitude", fontsize=11)
ax.tick_params(labelsize=9.5)
ax.set_title("Biogeographic regions used for the 29-population analysis",
             fontsize=14, fontweight="bold", pad=10)

handles = []
for region_id, (code, name, color, sampled) in REGIONS.items():
    handles.append(Patch(facecolor=color if sampled else "white", edgecolor="#333333",
                         linewidth=0.7, label=f"{code} — {name}"))
legend = ax.legend(handles=handles, title="Region code and name", loc="center left",
                   bbox_to_anchor=(1.015, 0.5), borderaxespad=0,
                   frameon=True, facecolor="white", edgecolor="#777777",
                   fontsize=9.2, title_fontsize=10.5, handlelength=1.6,
                   labelspacing=0.75, borderpad=0.9)
legend.get_title().set_fontweight("bold")

stem = OUT / "biogeographic_regions_10_sampled_vector_map"
fig.savefig(stem.with_suffix(".pdf"))
fig.savefig(stem.with_suffix(".png"), dpi=300)
plt.close(fig)
ds = None

with open(OUT / "biogeographic_regions_10_sampled_vector_legend.tsv", "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["code", "region_id", "region_name", "color", "sampled"])
    for region_id, (code, name, color, sampled) in REGIONS.items():
        writer.writerow([code, region_id, name, color, str(sampled).lower()])

print(f"Wrote {stem}.pdf and {stem}.png")

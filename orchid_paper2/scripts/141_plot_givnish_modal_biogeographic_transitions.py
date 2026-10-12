#!/usr/bin/env python3
"""Schematic maps of independent modal ancestral ranges, NOT BSM event estimates."""
import ast
import argparse
import csv
import json
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch, Patch, FancyArrowPatch
from matplotlib.lines import Line2D
from osgeo import ogr

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--calibration', choices=['givnish', 'zhang'], default='givnish')
parser.add_argument('--black-arrows', action='store_true')
parser.add_argument('--psa-retention-loop', action='store_true',
                    help='Show the PSA daughter retaining its parental range; not a sympatry assignment.')
args = parser.parse_args()
BASE = ROOT / f'19_biogeobears/unconstrained_29tip/{args.calibration}/max_range_3'
OUT = BASE / 'modal_range_transition_maps'
OUT.mkdir(exist_ok=True)
# Read constants without executing the original plotting script or overwriting it.
constants = {}
for statement in ast.parse((ROOT / 'scripts/98_plot_orchid_bioregions_vector_map.py').read_text()).body:
    if isinstance(statement, ast.Assign):
        for target in statement.targets:
            if isinstance(target, ast.Name) and target.id in ('REGIONS', 'LABEL_XY'):
                constants[target.id] = ast.literal_eval(statement.value)
REGIONS, LABEL_XY = constants['REGIONS'], constants['LABEL_XY']
XY = {v[0]: LABEL_XY[k] for k, v in REGIONS.items()}
BINS = [(1.5, 3 if args.calibration == 'zhang' else 2.5), (1, 1.5), (.5, 1), (0, .5)]
SPECIES_COLORS = {'PAM': '#D55E00', 'PAP': '#009E73', 'PSA': '#7B2CBF',
                  'shared ancestor': '#555555'}
with (ROOT / '17_herodotools/inputs/bioregions/tip_to_thesis_bioregion.tsv').open() as handle:
    SPECIES = {r['tip']: ('PSA' if r['tip'] == 'MDA' else 'PAM' if r['lineage'] == 'PAM' else 'PAP')
               for r in csv.DictReader(handle, delimiter='\t')}


class Node:
    def __init__(self, name='', length=0, children=()):
        self.name, self.length, self.children = name, length, list(children)


def parse_tree(text):
    """Parse the project's plain Newick, retaining every tip and branch length."""
    text, pos = text.strip(), 0
    def visit():
        nonlocal pos
        children = []
        if text[pos] == '(':
            pos += 1
            while True:
                children.append(visit())
                separator = text[pos]
                pos += 1
                if separator == ')':
                    break
                assert separator == ','
        start = pos
        while pos < len(text) and text[pos] not in ':,);':
            pos += 1
        name = text[start:pos]
        length = 0
        if pos < len(text) and text[pos] == ':':
            pos += 1
            start = pos
            while pos < len(text) and text[pos] not in ',);':
                pos += 1
            length = float(text[start:pos])
        return Node(name, length, children)
    root = visit()
    assert text[pos:] == ';'
    return root


def write_table(name, rows):
    with (OUT / name).open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


root = parse_tree((BASE / 'inputs/mcc_myr.newick').read_text())
nodes = []
def annotate(node, depth=0):
    node.depth = depth
    nodes.append(node)
    for child in node.children:
        annotate(child, depth + child.length)
    node.tips = frozenset().union(*(c.tips for c in node.children)) if node.children else frozenset([node.name])
annotate(root)
height = max(n.depth for n in nodes if not n.children)
modal = {}
with (BASE / 'all_internal_nodes_top3_ancestral_ranges.tsv').open() as handle:
    for row in csv.DictReader(handle, delimiter='\t'):
        if row['rank'] == '1':
            modal[frozenset(row['descendant_tips'].split(','))] = row
geography = {}
lines = (BASE / 'inputs/geography.data').read_text().splitlines()
codes = lines[0].split('(')[1].split(')')[0].split()
for line in lines[1:]:
    if line.strip():
        tip, values = line.split()
        geography[tip] = frozenset(c for c, value in zip(codes, values) if value == '1')
assert len(geography) == 29 and len(modal) == 28
for node in nodes:
    node.age = height - node.depth if node.children else 0
    if node.children:
        row = modal[node.tips]
        node.id, node.state, node.prob = row['node'], frozenset(row['state_code'].split('+')), float(row['probability'])
    else:
        node.id, node.state, node.prob = node.name, geography[node.name], 1.

audit, arrows, separators, retention_loops = [], [], [], []
def state_text(state):
    return '+'.join(sorted(state))
for node in nodes:
    if not node.children:
        continue
    assert len(node.children) == 2
    left, right = node.children
    p, a, b = node.state, left.state, right.state
    if a == p and b == p:
        kind = 'unchanged'
    elif a and b and not (a & b) and a | b == p:
        kind = 'vicariance-like partition'
    elif (a == p and not b & p) or (b == p and not a & p):
        kind = 'founder-like change'
    elif (a == p and b < p) or (b == p and a < p):
        kind = 'subset inheritance'
    elif a <= p and b <= p:
        kind = 'contraction/loss-like'
    elif a > p or b > p:
        kind = 'expansion-like change'
    else:
        kind = 'unresolved modal transition'
    bin_id = next(i for i, (low, high) in enumerate(BINS) if low <= node.age < high)
    row = dict(node=node.id, age_Ma=node.age, time_bin=f'{BINS[bin_id][0]:g}–{BINS[bin_id][1]:g}',
               parent_range=state_text(p), parent_probability=node.prob,
               daughter1=left.id, daughter1_range=state_text(a), daughter1_probability=left.prob,
               daughter2=right.id, daughter2_range=state_text(b), daughter2_probability=right.prob,
               classification=kind, minimum_modal_probability=min(node.prob, left.prob, right.prob),
               descendant_tips=','.join(sorted(node.tips)))
    audit.append(row)
    if args.psa_retention_loop:
        for child in node.children:
            if {SPECIES[t] for t in child.tips} == {'PSA'} and child.state == p:
                for region in sorted(p):
                    retention_loops.append(dict(row, bin_id=bin_id, region=region,
                        retained_daughter=child.id, loop_species='PSA',
                        interpretation='PSA daughter retains parental range; not a sympatry event assignment'))
    if kind == 'vicariance-like partition':
        for source in sorted(a):
            for target in sorted(b):
                separators.append(dict(row, bin_id=bin_id, source=source, target=target))
    else:
        for target in sorted((a | b) - p):
            colonizing = [c for c in node.children if target in c.state]
            colonizing_tips = frozenset().union(*(c.tips for c in colonizing))
            species_set = {SPECIES[t] for t in colonizing_tips}
            species = next(iter(species_set)) if len(species_set) == 1 else 'shared ancestor'
            for source in sorted(p):
                arrows.append(dict(row, bin_id=bin_id, source=source, target=target,
                                   arrow_species=species, descendant_species=','.join(sorted(species_set)),
                                   colonizing_tips=','.join(sorted(colonizing_tips))))
write_table('all_28_node_modal_transitions.tsv', audit)
if arrows:
    write_table('mapped_arrow_transitions.tsv', arrows)
summary = [dict(classification=k, n_nodes=v) for k, v in sorted(Counter(r['classification'] for r in audit).items())]
write_table('transition_counts.tsv', summary)
if args.psa_retention_loop:
    assert retention_loops, 'No PSA range-retaining daughter found; do not fabricate a loop'
    write_table('psa_range_retention_loop_audit.tsv', retention_loops)
cumulative = set(root.state)
active_by_bin, region_audit = [], []
for i, (low, high) in enumerate(BINS):
    previous = set(cumulative)
    for row in audit:
        if low <= row['age_Ma'] < high:
            for field in ('parent_range', 'daughter1_range', 'daughter2_range'):
                cumulative.update(row[field].split('+'))
    active_by_bin.append(set(cumulative))
    region_audit.append(dict(time_bin=f'{low:g}–{high:g}',
                             colored_regions=state_text(cumulative),
                             newly_colored_regions=state_text(cumulative - previous),
                             root_region=state_text(root.state)))
write_table('cumulative_region_coloring.tsv', region_audit)


def add_polygon(ax, rings, facecolor, edgecolor, lw, zorder):
    vertices, codes = [], []
    for ring in rings:
        points = [(p[0], p[1]) for p in ring]
        if len(points) < 3:
            continue
        if points[-1] != points[0]:
            points.append(points[0])
        vertices.extend(points)
        codes.extend([MplPath.MOVETO] + [MplPath.LINETO] * (len(points)-2) + [MplPath.CLOSEPOLY])
    if vertices:
        patch = PathPatch(MplPath(vertices, codes), facecolor=facecolor, edgecolor=edgecolor,
                          linewidth=lw, zorder=zorder, transform=ax.transData)
        # Extent is fixed explicitly below. Avoid add_patch's expensive automatic
        # bounds calculation over millions of coastline segments; geometry is unchanged.
        ax.add_artist(patch)


ds = ogr.Open(str(ROOT / '17_herodotools/inputs/bioregions/biogeographic_regions.gpkg'))
layer = ds.GetLayerByName('bioregions')
regions = [(feature.GetField('region_id'), json.loads(feature.GetGeometryRef().ExportToJson())) for feature in layer]
assert {r[0] for r in regions} == set(REGIONS)
background = []
for path in sorted((ROOT / '18_contemporary_migration/inputs/map_boundaries_geojson').glob('*.json')):
    background.extend(f['geometry'] for f in json.loads(path.read_text())['features'])
assert background, 'No country boundary features loaded'
def geometry(ax, geom, color, edge, lw, order):
    polygons = geom['coordinates'] if geom['type'] == 'MultiPolygon' else [geom['coordinates']]
    # One compound path retains ALL islands/rings without thousands of patches.
    add_polygon(ax, [ring for polygon in polygons for ring in polygon], color, edge, lw, order)

plt.rcParams.update({'font.size': 12, 'pdf.fonttype': 42})
fig, axes = plt.subplots(2, 2, figsize=(19, 14))
fig.subplots_adjust(left=.045, right=.77, bottom=.08, top=.92, wspace=.13, hspace=.18)
for i, ax in enumerate(axes.flat):
    print(f'Rendering period {i+1}/4', flush=True)
    ax.set_facecolor('#EAF4FA')
    for geom in background:
        geometry(ax, geom, '#FFFFFF', '#777777', .3, 1)
    for region_id, geom in regions:
        color = REGIONS[region_id][2] if REGIONS[region_id][0] in active_by_bin[i] else '#FFFFFF'
        geometry(ax, geom, color, '#666666', .32, 2)
    # Group repeats within species only; different species retain separate curves.
    grouped = {}
    for row in arrows:
        if row['bin_id'] != i:
            continue
        style = 'founder' if row['classification'] == 'founder-like change' else 'unresolved'
        key = row['source'], row['target'], style, row['arrow_species']
        grouped.setdefault(key, []).append(row)
    for (source, target, style, species), rows in grouped.items():
        x, y = XY[source], XY[target]
        rad = .25 if source < target else -.25
        # Reciprocal arrows are routed on opposite sides of their chord.
        if any(k[0] == target and k[1] == source for k in grouped):
            rad = .28
        parallel = sorted({k[3] for k in grouped if k[:2] == (source, target)})
        if len(parallel) > 1:
            rad += .16 * (parallel.index(species) - (len(parallel)-1)/2)
        ax.add_patch(FancyArrowPatch(x, y, connectionstyle=f'arc3,rad={rad}',
            arrowstyle='-|>', mutation_scale=17, linewidth=1.65,
            linestyle='-' if style == 'founder' else '--',
            color='#151515' if args.black_arrows else SPECIES_COLORS[species],
            shrinkA=8, shrinkB=8, zorder=4))
    for row in retention_loops:
        if row['bin_id'] != i:
            continue
        x, y = XY[row['region']]
        # Open loop east of Mindanao, with both ends pointing toward G and
        # a clear gap around its label. Route the body offshore, not across
        # the Philippine islands. This encodes retention, not travel.
        path = MplPath([(x + .8, y + .55),
                        (x + 7.0, y + 5.0),
                        (x + 7.0, y - 5.0),
                        (x + .8, y - .55)],
                       [MplPath.MOVETO] + [MplPath.CURVE4] * 3)
        ax.add_patch(FancyArrowPatch(path=path, arrowstyle='-|>',
            mutation_scale=17, linewidth=1.8, color=SPECIES_COLORS['PSA'], zorder=6))
    for row in separators:
        if row['bin_id'] == i:
            x, y = XY[row['source']], XY[row['target']]
            mx, my = (x[0]+y[0])/2, (x[1]+y[1])/2
            dx, dy = y[0]-x[0], y[1]-x[1]
            norm = (dx*dx+dy*dy)**.5
            ax.plot([mx-dy/norm*2, mx+dy/norm*2], [my+dx/norm*2, my-dx/norm*2], '--', color='black', lw=2)
    for region_id, (x, y) in LABEL_XY.items():
        if REGIONS[region_id][0] not in active_by_bin[i]:
            continue
        ax.text(x, y, REGIONS[region_id][0], ha='center', va='center', fontsize=10,
                weight='bold', zorder=5,
                bbox=dict(boxstyle='circle,pad=.15', facecolor='white', edgecolor='#555555', lw=.5))
    low, high = BINS[i]
    ax.set_title(f'{high:g}–{low:g} Ma', fontsize=17, weight='bold')
    ax.set_xlim(104, 158)
    ax.set_ylim(-21, 27)
    ax.set_aspect('equal')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.tick_params(labelsize=10)

handles = [Patch(facecolor=v[2], edgecolor='#555555', label=f'{v[0]}  {v[1]}') for v in REGIONS.values()]
legend = fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(.78, .88),
                    title='Biogeographic regions', frameon=True, fancybox=False,
                    fontsize=10, title_fontsize=12, labelspacing=1.0)
event_handles = ([Line2D([0],[0],color='#151515',lw=1.65,label='Founder-like change')]
                 if args.black_arrows else
                 [Line2D([0],[0],color=color,lw=2,label=species) for species,color in SPECIES_COLORS.items()])
if any(r['classification'] != 'founder-like change' for r in arrows):
    event_handles.append(Line2D([0],[0],color='#151515',lw=1.65,ls='--',label='Expansion / unresolved change'))
fig.legend(handles=event_handles,
           loc='upper left', bbox_to_anchor=(.78,.43), title='Transitions' if args.black_arrows else 'Arrow lineage', frameon=True, fancybox=False, fontsize=10)
lineage_note = ('Equal-width black arrows; no count labels.' if args.black_arrows else
                'Color: colonizing daughter lineage\nGray: mixed-species ancestral lineage')
if args.psa_retention_loop:
    lineage_note += '\nPurple loop: PSA retains ancestral range;\nnot dispersal or inferred sympatry.'
vicariance_note = ('No modal-state vicariance partitions.' if not separators else
                   'Dashed separators: vicariance-like splits.')
model = (BASE / 'selected_mcc_model.txt').read_text().strip()
note = ('Color accumulates from oldest to youngest.\nWhite land: not yet reached / unassigned.\n\n'
        'Arrow direction: ancestral → new region\n' + lineage_note + '\n\n' + vicariance_note +
        '\n\nLargest complete-range pie slice per node;\nindependent marginal states are not\n'
        'a jointly reconstructed event history.\n\nPanel age = ancestral branching age;\n'
        f'actual event timing is not estimated.\nModel: {model}; root age: {height:.2f} Ma\n'
        'See audit tables for probabilities.')
fig.text(.79,.29, note,
         fontsize=10, va='top', linespacing=1.5)
fig.suptitle(f'{args.calibration.capitalize()}: schematic ancestral-range transitions', fontsize=21, weight='bold', y=.975)
suffix = 'monochrome' if args.black_arrows else 'species'
if args.psa_retention_loop:
    suffix += '_psa_retention_loop'
stem = OUT / f'{args.calibration}_modal_ancestral_range_transitions_time_panels_{suffix}'
print('Saving PDF and PNG', flush=True)
fig.savefig(stem.with_suffix('.pdf'), bbox_inches='tight', pad_inches=.22)
fig.savefig(stem.with_suffix('.png'), dpi=240, bbox_inches='tight', pad_inches=.22)
plt.close(fig)
print(json.dumps({'counts': summary, 'arrows_before_aggregation':len(arrows), 'vicariance_separators':len(separators), 'output':str(stem)}, indent=2))

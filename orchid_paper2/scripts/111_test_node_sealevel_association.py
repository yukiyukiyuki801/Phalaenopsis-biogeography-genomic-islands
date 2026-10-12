#!/usr/bin/env python3
"""Test whether posterior divergence times are enriched in sea-level/climate phases."""

from pathlib import Path
from collections import defaultdict
import io

import numpy as np
import pandas as pd
from Bio import Phylo
from scipy.signal import find_peaks
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[2]
TREE_DIR = ROOT / "orchid_paper2/16_rad_beast/final_unconstrained_60m_summary"
SEA_FILE = ROOT / "orchid_paper2/results/biogeobears/sea_level/deboer2014-noaa.txt"
OUT = ROOT / "orchid_paper2/20_divergence_sealevel"
FIG = OUT / "figures"
TAB = OUT / "tables"
for directory in (OUT, FIG, TAB):
    directory.mkdir(parents=True, exist_ok=True)

PAM = set("AMB PAL PNG PON QL SAB SK SUL SUR TUA WJ".split())
PHASES = ["Glacial lowstand", "Deglaciation/rising", "Interglacial highstand",
          "Glaciation/falling"]
PHASE_COLORS = ["#2166AC", "#67A9CF", "#D6604D", "#F4A582"]
ISO_CLASSES = ["Glacial", "Intermediate", "Interglacial"]
N_NULL = 49999
SEED = 20261002


def load_climate():
    raw = pd.read_csv(SEA_FILE, sep="\t", comment="#")
    raw = raw.loc[:, ["age_calkaBP", "d18Obenth", "sealev"]].copy()
    raw["age_kyr"] = raw.age_calkaBP
    raw = raw.sort_values("age_kyr")
    # One-kyr resolution is sufficient relative to RAD node-age uncertainty.
    grid = np.arange(0, 5301, 1.0)
    sea = np.interp(grid, raw.age_kyr, raw.sealev)
    d18 = np.interp(grid, raw.age_kyr, raw.d18Obenth)
    climate = pd.DataFrame({"age_kyr": grid, "sealev": sea, "d18Obenth": d18})
    # Five-kyr smoothing suppresses numerical noise but retains orbital cycles.
    climate["sealev_smooth"] = climate.sealev.rolling(5, center=True,
        min_periods=1).mean()
    climate["d18O_smooth"] = climate.d18Obenth.rolling(5, center=True,
        min_periods=1).mean()
    # Detect individual highstand-to-highstand cycles. A 15-kyr minimum
    # separation and 0.5-m prominence retain the approximately 41-kyr early
    # Pleistocene cycles without treating numerical noise as a new cycle.
    sea_values = climate.sealev_smooth.to_numpy()
    d18_values = climate.d18O_smooth.to_numpy()
    peaks, _ = find_peaks(sea_values, distance=15, prominence=.5)
    # The present is an observed interglacial boundary rather than an
    # interior peak; include both record endpoints to retain boundary cycles.
    peaks = np.unique(np.r_[0, peaks, len(climate)-1])
    phase = np.full(len(climate), -1, dtype=int)
    iso = np.full(len(climate), -1, dtype=int)
    cycle_id = np.full(len(climate), -1, dtype=int)
    sea_cycle_position = np.full(len(climate), np.nan)
    d18_cycle_position = np.full(len(climate), np.nan)
    for cid, (left, right) in enumerate(zip(peaks[:-1], peaks[1:]), 1):
        trough = left + int(np.argmin(sea_values[left:right+1]))
        low_value = sea_values[trough]
        # Normalize each limb against its own bounding highstand so cycles
        # with different baselines and amplitudes are directly comparable.
        for start, stop, high_idx, middle_code in (
                (left, trough, left, 1),       # rising toward the present
                (trough, right, right, 3)):    # falling toward the present
            ii = np.arange(start, stop+1)
            amplitude = sea_values[high_idx] - low_value
            if amplitude <= 0:
                continue
            position = np.clip((sea_values[ii]-low_value)/amplitude, 0, 1)
            codes = np.where(position <= .25, 0,
                    np.where(position >= .75, 2, middle_code))
            phase[ii] = codes
            sea_cycle_position[ii] = position
            cycle_id[ii] = cid
        # Independent d18O check, normalized within the same sea-level cycle.
        ii = np.arange(left, right+1)
        lo, hi = d18_values[ii].min(), d18_values[ii].max()
        if hi > lo:
            position = (d18_values[ii]-lo)/(hi-lo)
            iso[ii] = np.where(position >= .75, 0,
                      np.where(position <= .25, 2, 1))
            d18_cycle_position[ii] = position
    climate["sea_phase_code"] = phase.astype(int)
    climate["isotope_class_code"] = iso.astype(int)
    climate["isotope_class"] = [ISO_CLASSES[x] if x >= 0 else "Outside complete cycle" for x in iso]
    climate["sea_phase"] = [PHASES[x] if x >= 0 else "Outside complete cycle" for x in phase]
    climate["cycle_id"] = cycle_id
    climate["sea_cycle_position"] = sea_cycle_position
    climate["d18O_cycle_position"] = d18_cycle_position
    climate.to_csv(TAB / "deboer2014_1kyr_classified_climate.tsv", sep="\t", index=False)
    return climate


def read_trees(scenario):
    path = TREE_DIR / f"{scenario}_unconstrained_100_posterior_29tip.newick"
    return list(Phylo.parse(path, "newick"))


def tree_nodes(tree):
    depths = tree.depths()
    height = max(depths[t] for t in tree.get_terminals())
    rows = []
    for node in tree.get_nonterminals(order="preorder"):
        tips = frozenset(t.name for t in node.get_terminals())
        age_ma = height - depths[node]
        if tips <= PAM:
            lineage = "PAM-internal"
        elif tips.isdisjoint(PAM):
            lineage = "PAP/PSA-internal"
        else:
            lineage = "Mixed-lineage"
        rows.append({"age_ma": age_ma, "age_kyr": age_ma*1000,
                     "tips": tips, "n_descendants": len(tips), "lineage": lineage})
    return rows


def climate_at_age(climate, age_kyr):
    idx = int(np.clip(np.rint(age_kyr), 0, len(climate)-1))
    row = climate.iloc[idx]
    return idx, int(row.sea_phase_code), int(row.isotope_class_code), row.sealev


def dataset_filter(rows, group):
    if group == "All nodes":
        return rows
    if group == "Last 1 Ma":
        return [x for x in rows if x["age_ma"] <= 1]
    if group == "Last 0.5 Ma":
        return [x for x in rows if x["age_ma"] <= .5]
    return [x for x in rows if x["lineage"] == group]


def summarize_codes(per_tree_rows, climate, scenario, groups, variable, labels, n_codes):
    rng = np.random.default_rng(SEED + (0 if scenario == "zhang" else 10000)
                                + (0 if variable == "sea_phase_code" else 20000))
    code_track = climate[variable].to_numpy(dtype=int)
    summaries, tree_summaries = [], []
    for group in groups:
        obs_matrix, exp_matrix = [], []
        null_sum = np.zeros((N_NULL, n_codes), dtype=float)
        usable = 0
        # A common circular-shift fraction is applied to every posterior tree
        # in a null replicate. Posterior trees describe uncertainty in one
        # history and are not 100 independent biological replicates; using an
        # independent shift for each tree would make the null artificially
        # narrow by averaging away the climate alignment.
        shift_fraction = rng.random(N_NULL)
        for tree_i, rows in enumerate(per_tree_rows, 1):
            # Do not assign nodes older than the 5.3-Ma proxy record to its
            # oldest boundary value. They are outside the testable interval.
            climate_max_kyr = float(climate.age_kyr.iloc[-1])
            selected = [x for x in dataset_filter(rows, group)
                        if x["age_kyr"] <= climate_max_kyr]
            if not selected:
                continue
            idx = np.array([climate_at_age(climate, x["age_kyr"])[0] for x in selected])
            root_idx = max(1, int(round(max(x["age_kyr"] for x in rows))))
            root_idx = min(root_idx, len(code_track)-1)
            valid_track = np.where(code_track >= 0)[0]
            start_idx = int(valid_track[0]); root_idx = min(root_idx, int(valid_track[-1]))
            idx = idx[(idx >= start_idx) & (idx <= root_idx) & (code_track[idx] >= 0)]
            if not len(idx):
                continue
            domain = code_track[start_idx:root_idx+1]
            relative_idx = idx-start_idx
            obs = np.bincount(code_track[idx], minlength=n_codes) / len(idx)
            duration = np.bincount(domain, minlength=n_codes) / len(domain)
            obs_matrix.append(obs); exp_matrix.append(duration); usable += 1
            shifts = np.floor(shift_fraction*len(domain)).astype(int)
            shifted_idx = (relative_idx[None, :] + shifts[:, None]) % len(domain)
            shifted_codes = domain[shifted_idx]
            for code in range(n_codes):
                null_sum[:, code] += np.mean(shifted_codes == code, axis=1)
            for code, label in enumerate(labels):
                tree_summaries.append({"scenario":scenario,"group":group,"tree":tree_i,
                    "category":label,"observed_fraction":obs[code],
                    "duration_fraction":duration[code],
                    "duration_adjusted_enrichment":obs[code]/duration[code] if duration[code] else np.nan,
                    "n_nodes":len(idx),"root_age_ma":root_idx/1000})
        obs_matrix = np.asarray(obs_matrix); exp_matrix = np.asarray(exp_matrix)
        null = null_sum / usable
        observed_stat = obs_matrix.mean(axis=0)
        for code, label in enumerate(labels):
            null_mean = null[:,code].mean()
            summaries.append({"scenario":scenario,"group":group,"category":label,
                "n_posterior_trees":usable,
                "median_nodes_per_tree":np.median([x["n_nodes"] for x in tree_summaries
                    if x["scenario"]==scenario and x["group"]==group and x["category"]==label]),
                "observed_mean_fraction":observed_stat[code],
                "observed_sd_fraction":np.std(obs_matrix[:,code],ddof=1),
                "observed_median_fraction":np.median(obs_matrix[:,code]),
                "observed_q025_fraction":np.quantile(obs_matrix[:,code],.025),
                "observed_q975_fraction":np.quantile(obs_matrix[:,code],.975),
                "expected_duration_fraction":exp_matrix[:,code].mean(),
                "circular_shift_null_mean_fraction":null_mean,
                "circular_shift_null_sd_fraction":np.std(null[:,code],ddof=1),
                "enrichment_ratio":observed_stat[code]/null_mean if null_mean else np.nan,
                "null_q025_fraction":np.quantile(null[:,code],.025),
                "null_q975_fraction":np.quantile(null[:,code],.975),
                "p_enriched_one_sided":(1+np.sum(null[:,code]>=observed_stat[code]))/(N_NULL+1),
                "p_depleted_one_sided":(1+np.sum(null[:,code]<=observed_stat[code]))/(N_NULL+1)})
    return pd.DataFrame(summaries), pd.DataFrame(tree_summaries)


def mcc_node_phase_probabilities(scenario, per_tree_rows, climate):
    mcc = Phylo.read(TREE_DIR/f"{scenario}_unconstrained_final_mcc_29tip.newick", "newick")
    out = []
    for node_i, node in enumerate(mcc.get_nonterminals(order="preorder"),1):
        tips = frozenset(t.name for t in node.get_terminals())
        climate_max_kyr = float(climate.age_kyr.iloc[-1])
        matched = [x for rows in per_tree_rows for x in rows
                   if x["tips"] == tips and x["age_kyr"] <= climate_max_kyr]
        phase = [climate_at_age(climate,x["age_kyr"])[1] for x in matched]
        iso = [climate_at_age(climate,x["age_kyr"])[2] for x in matched]
        phase = [x for x in phase if x >= 0]
        iso = [x for x in iso if x >= 0]
        base = {"scenario":scenario,"node":f"N{node_i}","n_descendants":len(tips),
            "descendants":",".join(sorted(tips)),"posterior_tree_occurrences":len(matched),
            "clade_frequency_in_100_trees":len(matched)/100,
            "median_age_ma":np.median([x["age_ma"] for x in matched]) if matched else np.nan}
        for i,label in enumerate(PHASES): base[f"P_{label}"] = np.mean(np.array(phase)==i) if phase else np.nan
        for i,label in enumerate(ISO_CLASSES): base[f"P_isotope_{label}"] = np.mean(np.array(iso)==i) if iso else np.nan
        out.append(base)
    return pd.DataFrame(out)


def plot_observed_expected(summary):
    data = summary[summary.group=="All nodes"].copy()
    fig, axes = plt.subplots(2,1,figsize=(11,8.5),sharex=True,sharey=True)
    x = np.arange(len(PHASES)); width=.34
    for ax, scenario in zip(axes,("zhang","givnish")):
        z=data[data.scenario==scenario].set_index("category").loc[PHASES]
        obs=z.observed_mean_fraction.to_numpy()*100; exp=z.circular_shift_null_mean_fraction.to_numpy()*100
        obs_sd=z.observed_sd_fraction.to_numpy()*100
        exp_sd=z.circular_shift_null_sd_fraction.to_numpy()*100
        ax.bar(x-width/2,obs,width,color=PHASE_COLORS,label="Observed nodes",edgecolor="none")
        ax.errorbar(x-width/2,obs,yerr=obs_sd,fmt="none",ecolor="black",capsize=3,lw=.8)
        ax.bar(x+width/2,exp,width,color="white",edgecolor="#444",hatch="//",label="Circular-shift expectation")
        ax.errorbar(x+width/2,exp,yerr=exp_sd,fmt="none",ecolor="#333",capsize=3,lw=.9)
        ax.set_title(f"{scenario.capitalize()} calibration",fontweight="bold")
        ax.set_ylabel("Nodes (%)")
        ax.spines[["top","right"]].set_visible(False)
        ax.legend(frameon=False,ncol=2,loc="upper right")
    axes[-1].set_xticks(x,PHASES,rotation=22,ha="right")
    fig.suptitle("Divergence-node occupancy of sea-level phases\n100 posterior trees per calibration",fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIG/"all_nodes_observed_vs_expected_sealevel_phases.pdf")
    fig.savefig(FIG/"all_nodes_observed_vs_expected_sealevel_phases.png",dpi=300)
    plt.close(fig)


def plot_enrichment(summary):
    data=summary.copy(); groups=list(data.group.drop_duplicates())
    fig,axes=plt.subplots(1,2,figsize=(14,8),sharey=True)
    for ax,scenario in zip(axes,("zhang","givnish")):
        z=data[data.scenario==scenario]
        positions=[]; labels=[]
        for gi,group in enumerate(groups):
            for pi,phase in enumerate(PHASES):
                row=z[(z.group==group)&(z.category==phase)].iloc[0]
                yy=gi*5+pi; positions.append(yy); labels.append(f"{group} — {phase}")
                ax.plot(row.enrichment_ratio,yy,"o",color=PHASE_COLORS[pi],ms=6)
        ax.axvline(1,color="#555",ls="--",lw=1)
        ax.set_title(f"{scenario.capitalize()} calibration",fontweight="bold")
        ax.set_xlabel("Enrichment ratio (observed / shifted expectation)")
        ax.set_yticks(positions)
        if ax is axes[0]:
            ax.set_yticklabels(labels,fontsize=8)
        else:
            ax.tick_params(axis="y",labelleft=False)
        ax.grid(axis="x",color="#ddd",lw=.5)
        ax.spines[["top","right"]].set_visible(False)
    fig.suptitle("Sea-level phase enrichment across divergence-node subsets",fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIG/"sealevel_phase_enrichment_all_subsets.pdf",bbox_inches="tight")
    fig.savefig(FIG/"sealevel_phase_enrichment_all_subsets.png",dpi=300,bbox_inches="tight")
    plt.close(fig)


def plot_climate_curve(climate, all_nodes):
    fig,axes=plt.subplots(2,1,figsize=(14,8),sharex=True,gridspec_kw={"height_ratios":[2,1]})
    age=climate.age_kyr/1000
    axes[0].plot(age,climate.sealev_smooth,color="#1769AA",lw=.65)
    axes[0].axhline(0,color="#777",ls="--",lw=.7)
    axes[0].set_ylabel("Sea level relative to present (m)")
    axes[0].spines[["top","right"]].set_visible(False)
    for si,scenario in enumerate(("zhang","givnish")):
        ages=[x["age_ma"] for rows in all_nodes[scenario] for x in rows]
        axes[1].hist(ages,bins=np.arange(0,5.31,.1),histtype="step",lw=1.4,
                     color=("#0072B2","#D55E00")[si],label=scenario.capitalize())
    axes[1].set_ylabel("Posterior node-age count")
    axes[1].set_xlabel("Millions of years before present")
    axes[1].invert_xaxis();axes[1].legend(frameon=False);axes[1].spines[["top","right"]].set_visible(False)
    fig.suptitle("Posterior divergence times and the de Boer et al. sea-level reconstruction",fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIG/"posterior_node_ages_and_sealevel_curve.pdf")
    fig.savefig(FIG/"posterior_node_ages_and_sealevel_curve.png",dpi=300)
    plt.close(fig)


def plot_isotope_observed_expected(summary):
    data=summary[summary.group=="All nodes"].copy()
    colors=["#2166AC","#BDBDBD","#D6604D"]
    primary=data[data.category=="Glacial"].set_index("scenario").loc[["zhang","givnish"]]
    holm=holm_adjust(primary.p_enriched_one_sided.to_numpy())
    holm_by_scenario=dict(zip(("zhang","givnish"),holm))
    fig,axes=plt.subplots(2,1,figsize=(9.5,8),sharex=True,sharey=True)
    x=np.arange(len(ISO_CLASSES));width=.34
    for ax,scenario in zip(axes,("zhang","givnish")):
        z=data[data.scenario==scenario].set_index("category").loc[ISO_CLASSES]
        obs=z.observed_mean_fraction.to_numpy()*100
        exp=z.circular_shift_null_mean_fraction.to_numpy()*100
        obs_sd=z.observed_sd_fraction.to_numpy()*100
        exp_sd=z.circular_shift_null_sd_fraction.to_numpy()*100
        ax.bar(x-width/2,obs,width,color=colors,edgecolor="none",label="Observed nodes")
        ax.errorbar(x-width/2,obs,yerr=obs_sd,fmt="none",ecolor="black",capsize=3,lw=.8)
        ax.bar(x+width/2,exp,width,color="white",edgecolor="#444",hatch="//",
               label="Circular-shift expectation")
        ax.errorbar(x+width/2,exp,yerr=exp_sd,fmt="none",
                    ecolor="#333",capsize=3,lw=.9)
        p=float(z.loc["Glacial","p_enriched_one_sided"])
        y=max(obs[0]+obs_sd[0],exp[0]+exp_sd[0])+2.0
        ax.text(0,y,f"one-sided P = {p:.4f}\nHolm P = {holm_by_scenario[scenario]:.4f}",
                ha="center",va="bottom",fontsize=8.5)
        ax.set_title(f"{scenario.capitalize()} calibration",fontweight="bold")
        ax.set_ylabel("Nodes (%)");ax.spines[["top","right"]].set_visible(False)
        ax.legend(frameon=False,ncol=2,loc="upper right")
        ax.set_ylim(0,62)
    axes[-1].set_xticks(x,ISO_CLASSES)
    fig.suptitle("Divergence-node occupancy of benthic δ18O climate classes\n"
                 "Mean ± SD across posterior trees and circular-shift null replicates",
                 fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIG/"all_nodes_observed_vs_expected_isotope_classes.pdf")
    fig.savefig(FIG/"all_nodes_observed_vs_expected_isotope_classes.png",dpi=300)
    plt.close(fig)


def add_bh_qvalues(table):
    table=table.copy()
    table["q_enriched_BH_within_scenario"]=np.nan
    table["q_depleted_BH_within_scenario"]=np.nan
    for scenario in table.scenario.unique():
        ii=np.where(table.scenario.to_numpy()==scenario)[0]
        for pcol,qcol in (("p_enriched_one_sided","q_enriched_BH_within_scenario"),
                          ("p_depleted_one_sided","q_depleted_BH_within_scenario")):
            p=table.iloc[ii][pcol].to_numpy();order=np.argsort(p);ranked=p[order]
            q=np.minimum.accumulate((ranked*len(ranked)/np.arange(1,len(ranked)+1))[::-1])[::-1]
            q=np.minimum(q,1);restored=np.empty_like(q);restored[order]=q
            table.loc[table.index[ii],qcol]=restored
    return table


def holm_adjust(p):
    p=np.asarray(p,dtype=float);m=len(p);order=np.argsort(p);ranked=p[order]
    adjusted=np.maximum.accumulate(ranked*(m-np.arange(m)))
    restored=np.empty_like(adjusted);restored[order]=np.minimum(adjusted,1)
    return restored


climate = load_climate()
groups = ["All nodes","Last 1 Ma","Last 0.5 Ma","PAM-internal","PAP/PSA-internal","Mixed-lineage"]
all_nodes, phase_summaries, phase_tree_summaries, iso_summaries, mcc_tables = {}, [], [], [], []
for scenario in ("zhang","givnish"):
    per_tree = [tree_nodes(t) for t in read_trees(scenario)]
    all_nodes[scenario] = per_tree
    summary, tree_summary = summarize_codes(per_tree,climate,scenario,groups,
        "sea_phase_code",PHASES,len(PHASES))
    phase_summaries.append(summary); phase_tree_summaries.append(tree_summary)
    iso_summary,_ = summarize_codes(per_tree,climate,scenario,groups,
        "isotope_class_code",ISO_CLASSES,len(ISO_CLASSES))
    iso_summaries.append(iso_summary)
    mcc_tables.append(mcc_node_phase_probabilities(scenario,per_tree,climate))

phase_summary=add_bh_qvalues(pd.concat(phase_summaries,ignore_index=True))
phase_tree_summary=pd.concat(phase_tree_summaries,ignore_index=True)
iso_summary=add_bh_qvalues(pd.concat(iso_summaries,ignore_index=True))
mcc_summary=pd.concat(mcc_tables,ignore_index=True)
phase_summary.to_csv(TAB/"sealevel_phase_enrichment_summary.tsv",sep="\t",index=False)
phase_tree_summary.to_csv(TAB/"sealevel_phase_ratios_by_posterior_tree.tsv",sep="\t",index=False)
iso_summary.to_csv(TAB/"isotope_glacial_class_enrichment_summary.tsv",sep="\t",index=False)
mcc_summary.to_csv(TAB/"mcc_nodes_phase_probabilities.tsv",sep="\t",index=False)
primary=pd.concat([
    phase_summary[(phase_summary.group=="All nodes")&(phase_summary.category=="Glacial lowstand")],
    iso_summary[(iso_summary.group=="All nodes")&(iso_summary.category=="Glacial")]
],ignore_index=True)
primary.to_csv(TAB/"primary_lowstand_glacial_hypothesis.tsv",sep="\t",index=False)
primary_iso=iso_summary[(iso_summary.group=="All nodes")&
                        (iso_summary.category=="Glacial")].copy()
primary_iso["p_holm_two_calibrations"]=holm_adjust(
    primary_iso.p_enriched_one_sided.to_numpy())
# An intersection-union test requires enrichment under both alternative clock
# calibrations; its valid P value is the larger of the two component P values.
primary_iso["p_calibration_robust_intersection_union"]=float(
    primary_iso.p_enriched_one_sided.max())
primary_iso.to_csv(TAB/"primary_isotope_glacial_circular_shift_test.tsv",
                   sep="\t",index=False)
plot_observed_expected(phase_summary)
plot_enrichment(phase_summary)
plot_climate_curve(climate,all_nodes)
plot_isotope_observed_expected(iso_summary)

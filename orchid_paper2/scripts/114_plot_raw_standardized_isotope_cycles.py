#!/usr/bin/env python3
"""Plot raw and within-cycle standardized benthic d18O classifications."""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D

ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/"orchid_paper2/20_divergence_sealevel"
data=pd.read_csv(BASE/"tables/deboer2014_1kyr_classified_climate.tsv",sep="\t")

labels=["Glacial","Intermediate","Interglacial"]
colors=np.array(["#2166AC","#BDBDBD","#D6604D"])
x=data.age_kyr.to_numpy()/1000
raw=data.d18O_smooth.to_numpy()
standard=data.d18O_cycle_position.to_numpy()
codes=data.isotope_class_code.to_numpy(dtype=int)


def colored_line(ax,y,lw=.9):
    points=np.column_stack([x,y])
    segments=np.stack([points[:-1],points[1:]],axis=1)
    lc=LineCollection(segments,colors=colors[codes[:-1]],linewidths=lw,
                      capstyle="butt",rasterized=False)
    ax.add_collection(lc);ax.autoscale_view()


fig,axes=plt.subplots(2,1,figsize=(14,8.2),sharex=True,
    gridspec_kw={"height_ratios":[1.25,1]})

colored_line(axes[0],raw,1.0)
# Every short dashed segment is the 25% or 75% threshold calculated from that
# cycle's own raw minimum and maximum; they need not share one absolute value.
for cid,z in data[data.cycle_id>=0].groupby("cycle_id",sort=True):
    lo=float(z.d18O_smooth.min());hi=float(z.d18O_smooth.max())
    xmin=z.age_kyr.min()/1000;xmax=z.age_kyr.max()/1000
    axes[0].hlines([lo+.25*(hi-lo),lo+.75*(hi-lo)],xmin,xmax,
                   colors=["#D6604D","#2166AC"],linestyles="--",
                   linewidth=.42,alpha=.55)
axes[0].set_ylabel("Raw benthic δ18O (‰)")
axes[0].set_title("Raw proxy values with cycle-specific classification thresholds",
                  fontweight="bold")

colored_line(axes[1],standard,1.0)
axes[1].axhline(.25,color="#D6604D",ls="--",lw=1.1,
               label="Interglacial/intermediate threshold (0.25)")
axes[1].axhline(.75,color="#2166AC",ls="--",lw=1.1,
               label="Intermediate/glacial threshold (0.75)")
axes[1].set_ylim(-.04,1.04)
axes[1].set_yticks([0,.25,.5,.75,1])
axes[1].set_ylabel("Within-cycle standardized δ18O")
axes[1].set_xlabel("Millions of years before present")
axes[1].set_title("Each cycle rescaled from its own minimum (0) to maximum (1)",
                  fontweight="bold")

phase_legend=[Line2D([0],[0],color=colors[i],lw=3,label=labels[i]) for i in range(3)]
threshold_legend=[Line2D([0],[0],color="#D6604D",ls="--",lw=1.1,label="0.25 threshold"),
                  Line2D([0],[0],color="#2166AC",ls="--",lw=1.1,label="0.75 threshold")]
axes[0].legend(handles=phase_legend+threshold_legend,ncol=5,frameon=False,
               loc="upper center",bbox_to_anchor=(.5,1.02),fontsize=9)
axes[1].legend(handles=phase_legend+threshold_legend,ncol=5,frameon=False,
               loc="upper center",bbox_to_anchor=(.5,1.02),fontsize=9)
for ax in axes:
    ax.spines[["top","right"]].set_visible(False)
    ax.grid(axis="x",color="#e6e6e6",lw=.45)
axes[1].invert_xaxis()
fig.suptitle("Cycle-specific standardization of the de Boer et al. benthic δ18O record",
             fontsize=16,fontweight="bold")
fig.tight_layout(rect=[0,0,1,.965])
fig.savefig(BASE/"figures/raw_vs_cycle_standardized_isotope.pdf")
fig.savefig(BASE/"figures/raw_vs_cycle_standardized_isotope.png",dpi=300)
plt.close(fig)

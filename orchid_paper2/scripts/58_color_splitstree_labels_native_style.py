#!/usr/bin/env python3
"""Change only taxon-label fill colors in native SplitsTree SVG sources."""
from pathlib import Path
import argparse,re,pandas as pd

COL={"PAM":"rgb(0,114,178)","PAP":"rgb(213,94,0)","PSA":"rgb(0,158,115)"}

def recolor(src,dst,label_species):
    text=Path(src).read_text()
    def repl(m):
        tag,label=m.group(1),m.group(2)
        species=label_species.get(label)
        if species:
            tag=re.sub(r"fill:rgb\([^;]+\)",f"fill:{COL[species]}",tag)
        return tag+label+"</text>"
    Path(dst).write_text(re.sub(r"(<text\b[^>]*>)([^<]+)</text>",repl,text))

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--base",default="orchid_paper2");ap.add_argument("--tmp",required=True)
    a=ap.parse_args();root=Path(a.base)/"results"/"neighbornet";tmp=Path(a.tmp);tmp.mkdir(parents=True,exist_ok=True)
    im=pd.read_csv(root/"individual_metadata.tsv",sep="\t")
    recolor(root/"splitstree_170_individuals_neighbornet.svg",tmp/"individuals.svg",dict(zip(im.id,im.species)))

if __name__=="__main__":main()

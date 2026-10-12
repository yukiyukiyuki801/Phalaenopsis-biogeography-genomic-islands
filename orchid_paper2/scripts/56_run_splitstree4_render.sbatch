#!/bin/bash
#SBATCH --job-name=splitstree4_nn
#SBATCH --partition=day
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=orchid_paper2/logs/splitstree4_nn.%j.out
#SBATCH --error=orchid_paper2/logs/splitstree4_nn.%j.err
set -euo pipefail
JAVA=orchid_paper2/software/SplitsTree4/jre/bin/java
JARS='orchid_paper2/software/SplitsTree4/jars/*'
OUT=orchid_paper2/results/neighbornet
for LEVEL in 170_individuals; do
 INPUT="$OUT/splitstree_${LEVEL}_distances.nex"
 CMD="$OUT/splitstree_${LEVEL}.commands"
 xvfb-run -a "$JAVA" -Djava.util.prefs.userRoot=/tmp/orchid_splitstree_prefs \
  -cp "$JARS" splitstree4.main.SplitsTree -g -s -v -c "$CMD"
done

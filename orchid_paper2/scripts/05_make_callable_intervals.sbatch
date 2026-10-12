#!/usr/bin/env bash
#SBATCH --job-name=orchid_callable
#SBATCH --partition=day
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --output=orchid_paper2/logs/callable_%j.out
#SBATCH --error=orchid_paper2/logs/callable_%j.err

set -euo pipefail
module load SAMtools/1.21-GCC-13.3.0

manifest=${MANIFEST:-orchid_paper2/00_metadata/run_manifest.tsv}
bamdir=orchid_paper2/04_bam
outdir=orchid_paper2/05_all_sites_vcf/callable
mkdir -p "$outdir"

tail -n +2 "$manifest" | awk -v d="$bamdir" '{print d "/" $2 ".sorted.bam"}' \
  > "$outdir/bam.list"

n=$(wc -l < "$outdir/bam.list")
min_dp=${MIN_DP:-3}
min_fraction=${MIN_SAMPLE_FRACTION:-0.30}
min_samples=$(awk -v n="$n" -v f="$min_fraction" 'BEGIN{x=int(n*f); if(x<n*f)x++; if(x<1)x=1; print x}')

echo "samples=$n min_dp=$min_dp min_samples=$min_samples" \
  > "$outdir/callable_settings.txt"

# Stream positions covered in any BAM and retain bases callable in the required
# number of samples. Consecutive passing bases are collapsed to BED intervals.
samtools depth -@ "$SLURM_CPUS_PER_TASK" -q 20 -Q 20 -d 0 \
  -f "$outdir/bam.list" \
  | awk -v dp="$min_dp" -v need="$min_samples" '
      function flush(){if(active) print chr, start-1, last}
      {
        pass=0;
        for(i=3;i<=NF;i++) if($i>=dp) pass++;
        keep=(pass>=need);
        if(keep){
          if(active && $1==chr && $2==last+1){last=$2}
          else {flush(); chr=$1; start=$2; last=$2; active=1}
        } else {flush(); active=0}
      }
      END{flush()}
    ' OFS='\t' > "$outdir/cohort_callable.bed"

awk '{bp += $3-$2} END{print "callable_bp\t" bp}' \
  "$outdir/cohort_callable.bed" > "$outdir/cohort_callable.summary.tsv"


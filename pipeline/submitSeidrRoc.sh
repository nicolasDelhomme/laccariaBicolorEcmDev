#!/bin/bash

# Load the tools
module load bioinfo-tools seidr-devel

# process the argument
indir=$(realpath ../data/seidr/backbone)
out=$(realpath ../data/seidr/roc)
gs=$(realpath ../goldStandard/Populus-tremula_KEGG-based-positive-gold-standard.tsv)
sh=$(realpath ../UPSCb-common/pipeline/runSeidrRoc.sh)
proj=u2019007
mail=nicolas.delhomme@umu.se

if [ ! -d $out ]; then
  mkdir -p $out
fi

# find the network files
for f in $(find $indir -name "*.sf"); do
#f="$indir/aggregated.sf"
  fnam=$(basename ${f/.sf/})
  
  # identify networks
  read -r -a algos <<< $(seidr view -H $f | grep "\[A\]" | cut -d" " -f3 | xargs)
  
  # for every algos
  for i in $(seq 0 $(expr ${#algos[@]} - 1)); do
     sbatch -o $out/${fnam}_roc.out -e $out/${fnam}_roc.err \
     -A $proj --mail-user=$mail $sh $f $gs $(expr $i + 1) \
     $out/${fnam}_${algos[$i]}_roc.tsv
  done
done

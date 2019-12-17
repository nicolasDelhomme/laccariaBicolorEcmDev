#!/bin/bash

# Load the tools
module load bioinfo-tools seidr-devel

# arguments
mail=nicolas.delhomme@umu.se
network=$(realpath ../data/seidr/aggregate/aggregated.sf)
out=$(realpath ../data/seidr/threshold)
sh=$(realpath ../UPSCb-common/pipeline/runSeidrThreshold.sh)
proj=u2019007

if [ ! -d $out ]; then
  mkdir -p $out
fi

sbatch -A $proj -o $out/threshold.out -e $out/threshold.err -J lbc-threshold --mail-user=$mail \
  $sh $network $out/threshold.txt

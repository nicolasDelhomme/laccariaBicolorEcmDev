#!/bin/bash -l

# define the thresholds
thresholds=( 2.33 2.05 1.88 1.75 1.64 1.55 1.48 1.41 1.34 1.28 )

# Load the tools
module load bioinfo-tools seidr-devel

# process the argument
network=$(realpath ../data/seidr/aggregate/aggregated.sf)
out=$(realpath ../data/seidr/backbone)
sh=$(realpath ../UPSCb-common/pipeline/runSeidrBackbone.sh)
proj=u2019007

if [ ! -d $out ]; then
  mkdir -p $out
fi

# submit
for i in {0..9}; do
  j=$(expr $i + 1)
  
  sbatch -A $proj -o $out/backbone-${j}-percent.out \
  -e $out/backbone-${j}-percent.err -J cold-backbone-${j} \
  $sh $network ${thresholds[$i]} $out/backbone-${j}-percent.sf
done

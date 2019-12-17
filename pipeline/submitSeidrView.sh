#!/bin/bash -l
module load bioinfo-tools seidr-devel

mail=nicolas.delhomme@umu.se
proj=u2019007

in=$(realpath ../data/seidr/aggregate/aggregated.sf)
out=$(realpath ../data/seidr/backbone/hard-threshold-dot3.sf)
sh=$(realpath runSeidrView.sh)
th=0.3

sbatch -A $proj --mail-user=$mail -o ${out/.sf/.out} -e ${out/.sf/.err} \
-J seidr-view $sh $in $out $th

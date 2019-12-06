#!/bin/bash -l

## be verbose and print
set -eux

proj=u2019007
mail=nicolas.delhomme@umu.se

## source functions
source ../UPSCb-common/src/bash/functions.sh

ref=$(realpath ../fungi/indices/salmon/Aspen-Laccaria-transcripts_salmon-version-dot-14-dot1.inx)
bind=/mnt:/mnt
img=/mnt/picea/projects/singularity/salmon-0.14.1.simg

## January
in=$(realpath ../data/trimmomatic)
out=$(realpath ../data/Salmon/Aspen-Laccaria)

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

## for every file
for f in $(find $in -name "*_trimmomatic_1.fq.gz"); do
  fnam=$(basename ${f/_1.fq.gz/})
  
  ## execute
 sbatch -A $proj --mail-user=$mail \
  -e $out/$fnam.err -o $out/$fnam.out -J salmon.$fnam \
  ../UPSCb-common/pipeline/runSalmon.sh -b $bind \
  -i $img $ref $f $in/${fnam}_2.fq.gz $out

done

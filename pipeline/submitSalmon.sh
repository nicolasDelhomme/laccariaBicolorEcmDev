#!/bin/bash -l

## be verbose and print
set -eux

proj=snic2019-8-101
mail=nicolas.delhomme@umu.se

## define a function
usage () {
    echo "The UPSCb env. var. needs to be set to your Git UPSCb checkout directory."
}

## source functions
source $UPSCb/src/bash/functions.sh

## process the argument
ref=/proj/uppstore2019033/reference/Potri03/indices/salmon/Ptrichocarpa_v3.0_210_transcript_salmon-v14dot1.inx
bind=/proj/uppstore2019033:/proj/uppstore2019033
img=/proj/uppstore2019033/singularity/salmon-0.14.1.simg

## check vars
if [ -z $UPSCb ]; then
    abort "The UPSCb var needs to be set."
fi

## January
in=/proj/uppstore2019033/results/January/trimmomatic
out=/proj/uppstore2019033/results/January/Salmon/Potri03

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

## for every file
cd $out
for f in $(find $in -name "*_trimmomatic_1.fq.gz"); do
  fnam=$(basename ${f/_1.fq.gz/})
  
  ## execute
 sbatch -A $proj --mail-user=$mail \
  -e $out/$fnam.err -o $out/$fnam.out -J salmon.$fnam \
  $UPSCb/pipeline/runSalmon.sh -b $bind \
  -i $img $ref $f $in/${fnam}_2.fq.gz $out

done

## March
in=/proj/uppstore2019033/results/March/trimmomatic
out=/proj/uppstore2019033/results/March/Salmon/Potri03

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

## for every file
cd $out
for f in $(find $in -name "*_trimmomatic_1.fq.gz"); do
  fnam=$(basename ${f/_1.fq.gz/})
  
  ## execute
  sbatch -A $proj --mail-user=$mail \
  -e $out/$fnam.err -o $out/$fnam.out -J salmon.$fnam \
  $UPSCb/pipeline/runSalmon.sh -b $bind \
  -i $img $ref $f $in/${fnam}_2.fq.gz $out

done






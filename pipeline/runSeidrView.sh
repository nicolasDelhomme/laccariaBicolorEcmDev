#!/bin/bash -l
#SBATCH -p core -n 1
#SBATCH --mem=6G
#SBATCH -t 1-00:00:00
#SBATCH --mail-type=ALL

USAGETXT=\
"
  Usage: $0 <seidr file> <out file> <threshold>

"

source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

isExec seidr

if [ $# -ne 3 ]; then
  abort "This script expects 3 arguments"
fi

if [ ! -f $1 ]; then
  abort "The first argument needs to be an existing file"
fi

if [ ! -d $(dirname $2) ]; then
  abort "The second argument directory needs to exist"
fi

#if [ $3 -gt 1 ] || [ $3 -lt 0 ] ; then
#  abort "The third argument needs to be in [0-1]"
#fi

# run
seidr view -t $3 -b $1 -o $2


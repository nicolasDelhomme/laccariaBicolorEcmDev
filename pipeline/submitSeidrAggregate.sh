#!/bin/bash

# project vars
account=SNIC2019-3-207
mail=nicolas.delhomme@umu.se
in=/home/d/delhomme/pfs/seidr/aspen-laccaria/sf
out=/home/d/delhomme/pfs/seidr/aspen-laccaria/aggregate

# tools
source ../UPSCb-common/src/bash/functions.sh
source /pfs/nobackup/home/b/bastian/seidr/build/sourcefile

# dir
if [ ! -d $out ]; then
  mkdir -p $out
fi

# submit
sbatch -A $account --mail-user=$mail -e $out/aggregate.err -o $out/aggregate.err \
-J lbed-aggregate ../UPSCb-common/pipeline/runSeidrAggregate.sh $out $in/*.sf 

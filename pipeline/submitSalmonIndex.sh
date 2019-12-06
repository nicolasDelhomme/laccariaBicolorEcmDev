#!/bin/bash

# fail on ERROR
set -eux

# load helpers
source ../UPSCb-common/src/bash/functions.sh

# vars
mail=nicolas.delhomme@umu.se
speciesList=(Arabidopsis-thaliana Aspen-Laccaria Birch-fungi Betula-pendula Picea-abies Populus-tremula)
OPTIONS="-p" # default to perfectHash creation

# usage
USAGETXT=\
"
  $0 <species>
  
  Valid species:
  ${speciesList[@]}
"

# test
if [ "$#" -ne "1" ]; then
  abort "This script expects 1 argument"
fi

case "$1" in
	    "Arabidopsis-thaliana")
	      in=/mnt/picea/storage/reference/Arabidopsis-thaliana/ARAPORT11/fasta/Araport11_all.201606.cdna.fasta.gz
	      out=/mnt/picea/storage/reference/Arabidopsis-thaliana/ARAPORT11/indices/salmon/Araport11_all.201606.inx
	    ;;
	    "Aspen-Laccaria")
	      in=/mnt/picea/storage/reference/Laccaria-bicolor/Lacbi2/fasta/Aspen_Laccaria_transcripts.fasta.gz
	      out=/mnt/picea/storage/reference/Laccaria-bicolor/Lacbi2/indices/salmon/Aspen-Laccaria-transcripts_salmon-version-dot-14-dot1.inx
	    ;;
	    "Betula-pendula")
	      in=/mnt/picea/storage/reference/Betula-pendula/v1.4/fasta/mRNA.fa.gz
	      out=/mnt/picea/storage/reference/Betula-pendula/v1.4/indices/salmon/Betula-pendula-mRNA_salmon-version-dot-14-dot1.inx
	    ;;
	    "Birch-fungi")
	      in=/mnt/picea/storage/reference/Phytophthora-cactorum/fasta/Birch-Fungi-transcript.fa.gz
	      out=/mnt/picea/storage/reference/Phytophthora-cactorum/indices/salmon/Birch-Fungi-transcript_salmon-version-dot-14-dot1.inx
	    ;;	    
	    "Populus-tremula")
	      in=/mnt/picea/storage/reference/Populus-tremula/v1.1/fasta/Potra01-mRNA.fa.gz
	      out=/mnt/picea/storage/reference/Populus-tremula/v1.1/indices/salmon/Potra01-mRNA.inx
	    ;;
	    "Picea-abies")
	      in=/mnt/picea/storage/reference/Picea-abies/v1.0/fasta/GenePrediction/phased/Pabies1.0-all.phase.gff3.CDS.fa
	      out=/mnt/picea/storage/reference/Picea-abies/v1.0/indices/salmon/Pabies1.0-all.phased.inx
	    ;;
		  \?) ## unknown flag
		  abort "Unknow species";;
esac

outdir=`dirname $out`

sbatch -e $outdir/salmonIndex.err -o $outdir/salmonIndex.out \
-A facility --mail-user $mail ../UPSCb-common/pipeline/runSalmonIndex.sh $OPTIONS $in $out


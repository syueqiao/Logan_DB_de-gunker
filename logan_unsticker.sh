#!/bin/bash
#script to run overlapping windows of input genomes against input DB
#accepts URLs of genomes in 2bit and zipped fasta format
#Dependency: gunzip, seqkit, diamond
#Output directories and names are ugly (for now) but it works

#usage
function usage {
  echo ""
  echo "Usage: ./logan_DB_de-gunker.sh https://ftp.edu/genome.fna.gz protein_DB.fa"
  echo ""
  echo "    Input parameters"
  echo "    -b    Input genome is in 2bit format"
  echo ""
  echo "    Output : genome.fna.pro"
  exit 1
}

#setup
BIT='false'
GNM=''
DB=''

while getopts ":b:h:" arg; do
  case "$arg" in
    b)
      BIT='true'
      ;;
     h)
      usage
      ;; 
    \?) #unrecognized option - show help
      echo "Input parameter not recognized"
      usage
      ;;
  esac
done
shift $((OPTIND-2))

#set CPU Parameters
THREADS='8'

#input genomes for scanning and parse names
GNM=$1
name="${GNM##*/}"
basename="${name%.*}"

#input query protein db
DB=$2


#create diamond db with sequences above
diamond makedb --in ${DB} -d ${DB}.dmnd

#diamond search
if [ $BIT = 'true' ]; then

twoBitToFa "${GNM}" -udcDir=. stdout | seqkit sliding -s 800 -W 1000 - | diamond blastx \
	  -q - \
	  -d ${DB}.dmnd \
	  --masking 0 \
	  --ultra-sensitive \
      --seed-cut 0.9 \
	  -p $THREADS -k1 \
	  -f 6 qseqid  qstart qend qlen qstrand \
	       sseqid  sstart send slen \
	       pident evalue cigar \
	       qseq_translated full_qseq full_qseq_mate \
	  > $basename.pro

else

wget -O - "${GNM}" | gunzip | seqkit sliding -s 800 -W 1000 - | diamond blastx \
	  -q - \
	  -d ${DB}.dmnd \
	  --masking 0 \
	  --ultra-sensitive \
      --seed-cut 0.9 \
	  -p $THREADS -k1 \
	  -f 6 qseqid  qstart qend qlen qstrand \
	       sseqid  sstart send slen \
	       pident evalue cigar \
	       qseq_translated full_qseq full_qseq_mate \
	  > $basename.pro
     
fi

#use seed-cut "Cutoff for masking low complexity seeds in bits/letter. The defaults are 0.9 for the fast mode, 0.8 for the default mode and 1.0 otherwise. (Supported since v2.0.12)"
#seed-cut parameter needs to be further investigated, but seems to filter out some false positives

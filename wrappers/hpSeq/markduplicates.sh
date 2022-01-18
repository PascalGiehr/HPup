#$ -cwd
#$ -j y
#$ -l mem_free=32G
#$ -V
#$ -S /bin/bash

name="markduplicates.sh"

#TODO: get parameters from configuration file
samtools=/TL/deep-share/archive00/software/bin/samtools

printHelp(){
 echo -e ""
 echo -e "Usage: $name <options>"
 echo -e ""
 echo -e " Mandatory:"
 echo -e "  -i read in, comma separated list of bam files"
 echo -e "  -o read out as text file"
 echo -e ""
 echo -e " Optional"
}


while getopts ":i:o:h" opt
do
 case "$opt" in
  h) printHelp; exit 0 ;;
  i) bamFile="${OPTARG}" ;;
  o) output="${OPTARG}" ;;
 esac
done


if [ -z "$bamFile" ] || [ -z "$output" ]
 then
  echo ""
  echo "ERROR ($name): All mandatory options must be set"
  printHelp
 exit 1
fi

#create tmp directory and files to replace multiple append to the same file
TMPWD=`mktemp --tmpdir=$TMPDIR -d duplicates.XXXXXXXXXX`

NDUP=$TMPWD/info1.txt
NWOB=$TMPWD/info2.txt
STARTPOS=$TMPWD/info3.txt
NDUPWOB=$TMPWD/info4.txt

#cut -f4|sort|uniq -c |mawk '{TOT+=($1-1)} END {print "number of duplicates:\t" TOT}'
#samtools view -F4 merge_replace.bam |cut -f2|sort|uniq -c |awk -vOFS='\t' '$2=="1024" {C+=$1} $2=="1040" {T+=$1} END {print "number of duplicates:\t" C+T}'
# merge all bam files and read only mapped read ( not flag four)
# calculate number of duplicates by using start position
# filter out reads with added wobble information and write statistics
$samtools merge - $(echo $bamFile | tr ',' ' ') \
|$samtools view -F4 - \
|tee >(cut -f2|sort|uniq -c |awk -vOFS='\t' '$2=="1024" {N+=$1} $2=="1040" {D+=$1} END {print "number of duplicates:\t" N+D}' > $NDUP) \
|grep -E "ZW:Z:[ACGT]{6}" \
|tee >(echo -e "number of reads with wobble:\t$(wc -l)" > $NWOB) \
     >(cut -f4|sort|uniq -c|echo -e "number of startposition groups with wobble:\t$(wc -l)" > $STARTPOS) \
|mawk -vOFS='\t' '$2==0 {print $0}'|cut -f4|sort|uniq -c |mawk '{TOT+=($1-1)} END {print "number of duplicates with wobble:\t" TOT}' > $NDUPWOB

cat $NDUP $NWOB $STARTPOS $NDUPWOB > $output

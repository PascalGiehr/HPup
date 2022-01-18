#$ -cwd
#$ -j y
#$ -l mem_free=32G
#$ -V
#$ -S /bin/bash

name="wobblemapRead.sh"


#TODO: take parameters from config file...
awk=/usr/bin/awk
java=/usr/lib/jvm/java-8-oracle/jre/bin/java
picardjar=/DEEP_fhgfs/projects/karln/hpseq/pipeline/local_pipelines/java/jar/picard.jar
samtools=/TL/deep-share/archive00/software/bin/samtools

#TODO: add linker sequence as input

printHelp(){
 echo -e ""
 echo -e "Usage: $name <options>"
 echo -e ""
 echo -e " Mandatory:"
 echo -e "  -i text file with found wobble sequence"
 echo -e "  -j bam file"
 echo -e "  -o result1 as text file"
 echo -e ""
 echo -e " Optional"
}


while getopts ":i:j:o:h" opt
 do
 case "$opt" in
  h) printHelp; exit 0 ;;
  i) wobble="${OPTARG}" ;;
  j) bamFile="${OPTARG}" ;;
  o) output="${OPTARG}" ;;
 esac
done


if [ -z "$wobble" ] || [ -z "$bamFile" ] || [ -z "$output" ]
 then
  echo ""
  echo "ERROR ($name): All mandatory options must be set"
  printHelp
 exit 1
fi

# add wobble flag and wobble sequence to bam file if the wobble sequence was found in the read
$samtools view $bamFile \
| awk -vOFS="\t" 'NR==FNR {a[$1]=$2;next} NF {print $0, (($1 in a)?"ZW:Z:"a[$1]:"")}' $wobble - \
| cat <($samtools view -H $bamFile) - \
| $samtools view -b - \
| $samtools sort -T $bamFile.tmp - > $output


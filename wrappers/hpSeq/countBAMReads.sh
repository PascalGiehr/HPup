#$ -cwd
#$ -j y
#$ -l mem_free=8G
#$ -V
#$ -S /bin/bash

name="countBAMReads.sh"

#TODO: get parameters from configuration file
samtools=/TL/deep-share/archive00/software/bin/samtools
awk=/usr/bin/awk

printHelp(){
 echo -e ""
 echo -e "Usage: $name <options>"
 echo -e ""
 echo -e " Mandatory:"
 echo -e "  -i comma separated list of bam files"
 echo -e "  -o output as text file"
 echo -e "  -n descriptive name"
 echo -e ""
 echo -e " Optional"
}


while getopts ":i:o:n:h" opt
do
 case "$opt" in
  h) printHelp; exit 0 ;;
  i) bamFiles="${OPTARG}" ;;
  o) output="${OPTARG}" ;;
  n) description="${OPTARG}" ;;
 esac
done


if [ -z "$bamFiles" ] || [ -z "$output" ] || [ -z "$description" ]
 then
  echo ""
  echo "ERROR ($name): All mandatory options must be set"
  printHelp
 exit 1
fi

count = 0
filter="-F260"  # 4 - unmapped read and 256 - not primary alignment

# replace coma with space
for x in ${bamFiles//,/ }; do
# count number of reads, whose are not unmapped and not primary alignment 
 count=$((count + $(samtools view $filter -c $x) ))
done
# write count in outputfile
echo  -e "$description\t$count" > $output

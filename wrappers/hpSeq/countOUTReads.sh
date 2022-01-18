#$ -cwd
#$ -j y
#$ -l mem_free=8G
#$ -V
#$ -S /bin/bash


name="countOUTPUTReads.sh"

printHelp(){
 echo -e ""
 echo -e "Usage: $name <options>"
 echo -e ""
 echo -e " Mandatory:"
 echo -e "  -i input text file from pileup"
 echo -e "  -o output statistics as text file"
 echo -e "  -n descriptive name"
 echo -e ""
 echo -e " Optional"
}


while getopts ":i:o:n:h" opt
do
 case "$opt" in
  h) printHelp; exit 0 ;;
  i) input="${OPTARG}" ;;
  o) output="${OPTARG}" ;;
  n) description="${OPTARG}" ;;
 esac
done


if [ -z "$input" ] || [ -z "$output" ] || [ -z "$description" ]
 then
  echo ""
  echo "ERROR ($name): All mandatory options must be set"
  printHelp
 exit 1
fi

#create tmp directory and files to replace multiple append to the same file
TMPWD=`mktemp --tmpdir=$TMPDIR -d conout.XXXXXXXXXX`

METCOV=$TMPWD/info1.txt
UNMCOV=$TMPWD/info2.txt
SITES=$TMPWD/info3.txt
# chr  position  base   total   Meth.   Unmeth.  type     F       P       M       U
# 1       3238380 G       1       0       1       GC      0       0       0       1
# calculate average conversion of methylated CGs and unmethylated GCs
# count number of CGs,GCs,...
zcat $input \
| tee >(mawk -vOFS='\t' '$7=="CG" {TOT+=$4;METHCG+=$5} END {print "average_methylation_of_CGs", (METHCG/TOT)}' > $METCOV) \
| tee >(mawk -vOFS='\t' '$7=="GC" {TOT+=$4;METHGC+=$5} END {print "average_methylation_of_GCs", (METHGC/TOT)}' > $UNMCOV) \
| cut -f7 \
| sort \
| uniq -c \
| awk -v OFS='\t' -vlabel="$description" '
    {
    print label"_"$2, $1
    }
   ' \
> $SITES

cat $METCOV $UNMCOV $SITES > $output

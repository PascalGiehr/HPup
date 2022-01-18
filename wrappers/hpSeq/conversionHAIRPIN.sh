#$ -cwd
#$ -j y
#$ -l mem_free=32G
#$ -V
#$ -S /bin/bash

name="conversion2HAIRPIN.sh"

python3=/opt/bioinfo/bin/python3

# TODO: add linker sequence as input

printHelp(){
 echo -e ""
 echo -e "Usage: $name <options>"
 echo -e ""
 echo -e " Mandatory:"
 echo -e "  -i compressed fastq file"
 echo -e "  -x conversion file within number of hp's and conversion rates"
 echo -e "  -y header of read with found wobble sequence"
 echo -e ""
 echo -e " Optional"
 echo -e "  -c config file"
 echo -e "  -I path of InstallFolder"
}


while getopts ":i:x:y:c:I:h" opt
do
 case "$opt" in
  h) printHelp; exit 0 ;;
  i) fastqFile="${OPTARG}" ;;
  x) conversion="${OPTARG}" ;;
  y) wobble="${OPTARG}" ;;
  c) configFile="${OPTARG}" ;;
  I) installFolder="${OPTARG}" ;;
 esac
done


if [ -z "$fastqFile" ] || [ -z "$conversion" ]|| [ -z "$wobble" ] || [ -z "$configFile" ] || [ -z "$installFolder" ]
 then
  echo ""
  echo "ERROR ($name): All mandatory options must be set"
  printHelp
 exit 1
fi


#create tmp directory and files to replace multiple append to the same file
TMPWD=`mktemp --tmpdir=$TMPDIR -d conversion.XXXXXXXXXX`

NLINKER=$TMPWD/info1.txt
CONCYT=$TMPWD/info2.txt
BSCON=$TMPWD/info3.txt

# TODO: make undependent on linker sequence
# IUPAC nucleotide code Y= C or T and D = A or G or T
linkerGrep="$(python3 $installFolder/wrappers/getArgumentfromConfig.py $configFile Analysis linkerGrep)"
wobbleSed="$(python3 $installFolder/wrappers/getArgumentfromConfig.py $configFile Analysis linkerSed)"
id="$(python3 $installFolder/wrappers/getArgumentfromConfig.py $configFile Analysis fastqid)"

#id="@HISEQ"
# concat all fastq files and use mawk to get only the header and sequence of each read
# get reads with searched linker sequence by using grep and the given linker sequence
# get header and sequence of wobble position (see wobble output) and get number of linker (see conversion output)
# use sed to filter out the first and last group of the linker sequence
# $() execute
# create statistics, e.g. count found number of C's and T's to calculate the bisulfite conversion rate
#rm $conversion

zcat $(echo $fastqFile | tr ',' ' ') \
| mawk -vOFS="\t" -v RS="$id" -vFS="\n" 'NR>1{print RS$1,$2}' \
| grep -E "$linkerGrep" \
| tee >(sed -e 's/@\([^ \t]*\).*'$wobbleSed'.*/\1\t\3\4/' -e 's/\([^\n]\)/\1/g' > $wobble) \
      >(wc -l|awk -v OFS='\t' '{print "numberLinker", $1}' > $NLINKER) \
| sed -e 's/.*'$wobbleSed'.*/\1\4/' -e 's/\([^\n]\)/\1\n/g' \
| sort \
| uniq -c \
| tee >(awk -v OFS='\t' '$2=="C" {C+=$1} $2=="T" {T+=$1} END {print "numberUnconvertedCytosine", C"\n"  "numberAllCytosine", (C+T) }' >$CONCYT) \
| awk -v OFS='\t' '
    $2=="C" {C+=$1}
    $2=="T" {T+=$1}
    END {print "bisulfiteConversionRate", 1-(C/T)}
' \
> $BSCON
#echo $NLINKER $CONCYT $BSCON 
cat $NLINKER $CONCYT $BSCON > $conversion
rm -r $TMPWD 

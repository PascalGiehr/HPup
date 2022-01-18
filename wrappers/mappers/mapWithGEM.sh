#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l mem_free=64G
#$ -V

name="mapWithGEM.array.sh"
TMPDIR=/DEEP_fhgfs/tmp
readGroupString="ID=library_1_1,LB=library_1,SM=library"
uniqueHitsOnly=0

#TODO: needs an update if there should be a general config file for software

export PATH=/TL/deep-share/archive00/software/bin:$PATH
#LONGLINECUTOFF=500000
LONGLINECUTOFF=100000

printHelp(){
 echo -e ""
 echo -e "Usage: $name <options>"
 echo -e ""
 echo -e " Mandatory:"
 echo -e "  -a file1"
 echo -e "  -o outputfile"
 echo -e "  -r REFERENCE"
 echo -e ""
 echo -e " Optional"
 echo -e "  -b file2"
 echo -e "  -g readGroup (default: $readGroupString)"
 echo -e "  -u \t\tif only unique hits are to be considered"
 echo -e "  -t TMPDIR\tdefault: $TMPDIR"
}

while getopts ":a:b:o:r:g:t:uh" opt
do
 case "$opt" in
  h) printHelp; exit 0 ;;
  a) file1="$OPTARG" ;;
  b) file2="$OPTARG" ;;
  o) outputFile="$OPTARG" ;;
  r) reference="$OPTARG" ;;
  g) readGroupString="$OPTARG" ;;
  u) uniqueHitsOnly=1 ;;
  t) TMPDIR="$OPTARG" ;;
  *) printHelp; exit 1 ;;
 esac
done

if [ -z "$file1" ] || [ -z "$outputFile" ] || [ -z "$reference" ]
then
 echo ""
 echo "ERROR ($name): All mandatory options must be set"
 printHelp
 exit 1
fi


# qsub -t 1-`find $inputRoot -name "*_R1_*q.gz" |wc -l` -o log.txt mapWitGEM.array.sh $inputRoot $outputRoot $reference

echo `date` "LOGG "$OUTPUTNAME": START"
 
TMPWD=`mktemp --tmpdir=$TMPDIR -d mapGEM.XXXXXXXXXX`

set -o pipefail



#TODO: gem uses the whole header name of a sequence in an index as the target name in the output
#TODO: >1 dna  becomes "1 dna" in the bam file

samtoolsFilter="samtools view -F256"


if [ ! -f "$file2" ]
then
 if [ "$uniqueHitsOnly" -eq 1 ]
 then
  extraOptions="$extraOptions --unique-mapping"
#  samtoolsFilter="samtools view -F256"
 fi

 
 #TODO: gem sets the mapped in proper pair flag for single end reads, picard tools does not like this 
 ER1=$TMPWD/read1.fq

 gunzip -c $file1 > $ER1 || exit 1 
 
 touch $TMPWD/long.gi

 gem-mapper -I $reference -q offset-33 $extraOptions -T 1 -i $ER1 |
 mawk -vFS='\t' -vOFS='\t' -vCUTOFF=$LONGLINECUTOFF 'length($0)>CUTOFF {$4=0;$5="-";print $1 >"'$TMPWD/long.gi'"} {print}' |
 gem-2-sam --read-group $readGroupString -l -I $reference -q offset-33 -T 1|
 awk -vOFS='\t' -vFS='\t' '$1!~/^@/ {$2=$2-and($2,66)} {print}' |
 $samtoolsFilter -bS - | 
 samtools sort -m 5000000000 -O bam -T $TMPWD/bamSort - > $outputFile || exit 3

 

 echo `date` "LOGG "$OUTPUTNAME": Removed `cat $TMPWD/long.gi |wc -l` pairs with extremely many hits"
else

 if [ "$uniqueHitsOnly" -eq 1 ]
 then
  extraOptions="$extraOptions --unique-pairing"
#  samtoolsFilter="samtools view -F256"
 fi
 ER1=$TMPWD/read1.fq
 ER2=$TMPWD/read2.fq

 #TODO: change long lines to missing alignments by setting $4 to 0 and $5 to -. Then pipe directly into gem-2-sam. Both here and for paired end. IS THIS THE SAME FOR PAIRED END READS???

 gunzip -c $file1 > $ER1 || exit 4 
 gunzip -c $file2 > $ER2 || exit 5

 gem-mapper -I $reference -q offset-33 -p -b $extraOptions -T 1 -1 $ER1 -2 $ER2 |tee $TMPWD/out.gem |mawk -vCUTOFF=$LONGLINECUTOFF 'length($0)>CUTOFF {print $1}' > $TMPWD/long.gi || exit 6 

 echo `date` "LOGG "$OUTPUTNAME": Removed `cat $TMPWD/long.gi |wc -l` pairs with extremely many hits"

 grep -w -vFf $TMPWD/long.gi $TMPWD/out.gem | gem-2-sam --read-group $readGroupString -l -I $reference -q offset-33 -T 1|$samtoolsFilter -bS - | samtools sort -m 5000000000 -O bam -T $TMPWD/bamSort - > $outputFile || exit 7 

fi

rm -rf $TMPWD
echo `date` "LOGG "$OUTPUTNAME": DONE"

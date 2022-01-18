#$ -cwd
#$ -j y
#$ -l mem_free=8G
#$ -V
#$ -S /bin/bash

name="amplicons.sh"

gsnap=/TL/deep-share/archive00/software/bin/gsnap
samtools=/TL/deep-share/archive00/software/bin/samtools

tmpWd=/DEEP_fhgfs/projects/karln/hpseq/pipeline/data
padSize=100
fastaFile=$tmpWd/Spikes.ref.upper.fa

printHelp(){
echo -e ""
echo -e "Usage: $name <options>"
echo -e ""
echo -e " Mandatory:"
echo -e "  -i read1 in, comma separated list of files"
echo -e "  -o read2 out as text file"
echo -e ""
echo -e " Optional"
}


while getopts ":i:o:h" opt
do
case "$opt" in
h) printHelp; exit 0 ;;
i) file1="${OPTARG}" ;;
o) file1out="${OPTARG}" ;;
r) fastaFile="${OPTARG}" ;;
esac
done


if [ -z "$file1" ] || [ -z "$file1out" ]
then
echo ""
echo "ERROR ($name): All mandatory options must be set"
printHelp
exit 1
fi

mkdir -p $tmpWd/db
set -o pipefail

$gsnap --nthreads=8 --npaths=1 --db=repeats --dir=$tmpWd/db --sam-use-0M --format=sam --gunzip --mode=cmet-stranded $file1 |$samtools view -u -F 260 - | $samtools sort - > $file1out

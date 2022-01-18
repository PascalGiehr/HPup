#$ -cwd
#$ -j y
#$ -l mem_free=32G
#$ -V
#$ -S /bin/bash

name="replaceSeqmarkDup.sh"

#TODO: get parameters from configuration file
jarFolder=/DEEP_fhgfs/projects/karln/hpseq/pipeline/local_pipelines/java/jar
java=/usr/lib/jvm/java-8-oracle/jre/bin/java
samtools=/TL/deep-share/archive00/software/bin/samtools

#TODO: add linker sequence as input

printHelp(){
echo -e ""
echo -e "Usage: $name <options>"
echo -e ""
echo -e " Mandatory:"
echo -e "  -i fastq file"
echo -e "  -j bam file"
echo -e "  -o result1 as bam file"
echo -e ""
echo -e " Optional"
}


while getopts ":i:j:o:h" opt
do
case "$opt" in
h) printHelp; exit 0 ;;
i) fastqFile="${OPTARG}" ;;
j) bamFile="${OPTARG}" ;;
o) file1out="${OPTARG}" ;;
esac
done


if [ -z "$fastqFile" ] || [ -z "$bamFile" ] || [ -z "$file1out" ]
then
echo ""
echo "ERROR ($name): All mandatory options must be set"
printHelp
exit 1
fi

# replace the sequence from the bam files with the sequence from the fastq files to have informations about the methylation
$java -cp .:$jarFolder/htsjdk-1.130.jar:$jarFolder/bamUtils.jar tools.bam.bamUtils hpReplaceSeq $fastqFile $bamFile \
| cat <($samtools view -H $bamFile) - \
| $samtools view -b - \
| $samtools sort -T $bamFile.tmp - > $file1out

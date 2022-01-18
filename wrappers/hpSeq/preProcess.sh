#$ -cwd
#$ -j y
#$ -l mem_free=32G
#$ -V
#$ -S /bin/bash

name="preProcess.sh"
TMPDIR=/DEEP_fhgfs/tmp

#TODO: get parameters from configuration file

flash=/TL/deep-share/archive00/software/bin/flash
cutadapt=/TL/deep-share/archive00/software/bin/cutadapt
java=/usr/lib/jvm/java-8-oracle/jre/bin/java
fastqUtils=/DEEP_fhgfs/projects/karln/hpseq/pipeline/local_pipelines/java/jar/fastqUtils.jar

export PYTHONPATH=$PYTHONPATH:/TL/deep-share/archive00/software/lib/python2.7/site-packages

nrCpu=3
hpf="GGGTTT[AGT][AGT][AGT]T[AGT][AGT][AGT]AGGTTT"


printHelp(){
 echo -e ""
 echo -e "Usage: $name <options>"
 echo -e ""
 echo -e " Mandatory:"
 echo -e "  -a read1 in"
 echo -e "  -b read2 in"
 echo -e "  -x read1 out gzipped"
 echo -e "  -y read2 out gzipped"
 echo -e "  -o OUTPUTDIR"
 echo -e ""
 echo -e " Optional"
 echo -e "  -L exp\tlinker sequence regular expression (default: $hpf)"
 echo -e "  -t TMPDIR\tdefault: $TMPDIR"
}


while getopts ":a:b:x:y::L:l:t:h" opt
do
 case "$opt" in
  h) printHelp; exit 0 ;;
  a) file1="${OPTARG}" ;;
  b) file2="${OPTARG}" ;;
  x) file1out="${OPTARG}" ;;
  y) file2out="${OPTARG}" ;;
  L) hpf="${OPTARG}" ;;
  t) TMPDIR="$OPTARG" ;;
  *) printHelp; exit 1 ;;
 esac
done

if [ -z "$file1" ] || [ -z "$file2" ] || [ -z "$file1out" ] || [ -z "$file2out" ] || [ -z "$hpf" ]
then
 echo ""
 echo "ERROR ($name): All mandatory options must be set"
 printHelp
 exit 1
fi

echo `date` "LOGG "$name": START"

#setup work directory
tmpWd=$( mktemp --tmpdir=$TMPDIR -d ${name%.sh}.XXXXXXXXXX ) || exit 13
#inputFolder=$tmpWd/00_input
joinFolder=$tmpWd/05_join
splitFolder=$tmpWd/10_split
mergeFolder=$tmpWd/15_merge

mkdir -p $joinFolder $splitFolder $mergeFolder

#generate N string
hpfN=$(echo $hpf |sed -e 's/\[[^]]*\]/N/g')
hpfLength=${#hpfN}
hprN=$(echo $hpfN |tr 'ACGTN' 'TGCAN' |rev)



#get reads and extract them




#TODO: This is now generating a lot of files... and temporary files are compressed... not nice
#TODO: filter out short reads


$cutadapt -n 7 -a $hpfN $file1 2> /dev/null > $mergeFolder/seq_R1.fq ||exit 5
$cutadapt -n 7 -a $hprN $file2 2> /dev/null > $mergeFolder/seq_R2.fq || exit 10

$java -Xmx1G -jar $fastqUtils removeWindowPaired $mergeFolder/seq_R1.fq $mergeFolder/seq_R2.fq 0 20 $file1out $file2out ||exit 15

rm -r $tmpWd

echo `date` "LOGG DONE: " $name

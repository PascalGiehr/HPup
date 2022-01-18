#$ -cwd
#$ -j y
#$ -l mem_free=32G
#$ -V
#$ -S /bin/bash

name="preProcess.sh"
TMPDIR=/DEEP_fhgfs/tmp

flash=/TL/deep-share/archive00/software/bin/flash
cutadapt=/TL/deep-share/archive00/software/bin/cutadapt
java=/TL/deep-share/archive00/software/packages/java/jre1.7.0_60/bin/java
fastqUtils=/TL/deep-share/archive00/software/lib/jars/fastqUtils.jar

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


#merge reads
$flash -q -M 1000  -d $joinFolder -t $nrCpu $file1 $file2  || exit 1

#split on linker
#TODO: generate stats here....

#if the files are empty
touch $splitFolder/out.extendedFrags_R1.fq $splitFolder/out.extendedFrags_R2.fq

cat $joinFolder/out.extendedFrags.fastq | awk --re-interval 'BEGIN{
 RS = "'`cat $joinFolder/out.extendedFrags.fastq | head -1 |cut -f 1 -d ':'`'";
 hpl = "[AGCTN]{30,600}'$hpf'";
 fname="'$splitFolder'/out.extendedFrags";
 ll = '$hpfLength';
}
NR>1 && match($3,hpl) {
 rs=RSTART; re=RLENGTH;
 printf "%s%s %s\n%s\n%s\n%s\n",RS,$1,$2,substr($3,rs,re-ll),$4,substr($5,rs,re-ll) >  fname"_R1.fq";
 gsub("1:", "2:", $2)
 printf "%s%s %s\n%s\n%s\n%s\n",RS,$1,$2,substr($3,re,length($3)),$4,substr($5,re,length($5)) > fname"_R2.fq";
}' || exit 2

#TODO: This is now generating a lot of files... and temporary files are compressed... not nice
#TODO: filter out short reads
cat $splitFolder/out.extendedFrags_R1.fq <($cutadapt -a $hpfN  $joinFolder/out.notCombined_1.fastq 2> /dev/null) > $mergeFolder/seq_R1.fq ||exit 5

cat <($java -Xmx1G -jar $fastqUtils reverseComplement $splitFolder/out.extendedFrags_R2.fq) <($cutadapt -a $hprN $joinFolder/out.notCombined_2.fastq 2> /dev/null) > $mergeFolder/seq_R2.fq || exit 10

$java -Xmx1G -jar $fastqUtils removeWindowPaired $mergeFolder/seq_R1.fq $mergeFolder/seq_R2.fq 0 20 $file1out $file2out ||exit 15

#rm -r $tmpWd

echo `date` "LOGG DONE: " $name

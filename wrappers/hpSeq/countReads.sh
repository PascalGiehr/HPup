#$ -cwd
#$ -j y
#$ -l mem_free=8G
#$ -V
#$ -S /bin/bash

name="countReads.sh"

#zcat {inputfile} | echo $((`wc -l`/4)) > {outputfile}

printHelp(){
 echo -e ""
 echo -e "Usage: $name <options>"
 echo -e ""
 echo -e " Mandatory:"
 echo -e "  -i read1 in,comma separated list of fastq file"
 echo -e "  -o read2 out as text file"
 echo -e "  -n descriptive name"
 echo -e ""
 echo -e " Optional"
}


while getopts ":i:o:n:h" opt
do
 case "$opt" in
  h) printHelp; exit 0 ;;
  i) fastq="${OPTARG}" ;;
  o) output="${OPTARG}" ;;
  n) description="${OPTARG}" ;;
 esac
done


if [ -z "$fastq" ] || [ -z "$output" ] || [ -z "$description" ]
 then
  echo ""
  echo "ERROR ($name): All mandatory options must be set"
  printHelp
  exit 1
fi

# merge all fastq files and count lines, all four lines are one read
zcat $(echo $fastq | tr ',' ' ') | echo -e "$description\t$(($(wc -l)/4))" > $output


#$ -cwd
#$ -j y
#$ -l mem_free=8G
#$ -V
#$ -S /bin/bash

name="mergeLink.sh"

printHelp(){
echo -e ""
echo -e "Usage: $name <options>"
echo -e ""
echo -e " Mandatory:"
echo -e "  -i read in"
echo -e "  -o read out"
echo -e ""
echo -e " Optional"
}


while getopts ":i:o:h" opt
do
case "$opt" in
h) printHelp; exit 0 ;;
i) input="${OPTARG}" ;;
o) output="${OPTARG}" ;;
esac
done

if [ -z "$input" ] || [ -z "$output" ]
then
echo ""
echo "ERROR ($name): All mandatory options must be set"
printHelp
exit 1
fi

# TODO use awk if first six keys are equal sum up counts and delete one line
# merge all input text files which are in a comma separated list
sort -k1,1 -k2,2n -k3,3n -k4,4n -k5,5 -k6,6 $(echo $input |tr ',' ' ') \
| awk -vOFS='\t' '
    $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6!=prev
    {if(NR>1)print prev,count;
    prev=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6;
    count=0;}
    {count+=$7}
    END{print prev,count}
' \
| awk -vOFS='\t' '
    {print $1,$2,$3,$4,($5$6),$7}
' \
> $output

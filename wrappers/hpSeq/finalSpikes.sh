#$ -cwd
#$ -j y
#$ -l mem_free=32G
#$ -V
#$ -S /bin/bash

name="finalSpikes.sh"

mawk=/usr/bin/mawk
bedtools=/TL/deep-share/archive00/software/bin/bedtools
python3=/usr/bin/python3

printHelp(){
echo -e ""
echo -e "Usage: $name <options>"
echo -e ""
echo -e " Mandatory:"
echo -e "  -a forward merged txt file"
echo -e "  -b reverse merged txt file"
echo -e "  -o read out as text file"
echo -e ""
echo -e " Optional"
echo -e "  -c config file"
echo -e "  -I path of InstallFolder"
}

while getopts ":a:b:o:c:I:h" opt
do
 case "$opt" in
  h) printHelp; exit 0 ;;
  a) forward="${OPTARG}" ;;
  b) reverse="${OPTARG}" ;;
  o) file1out="${OPTARG}" ;;
  c) configFile="${OPTARG}" ;;
  I) installFolder="${OPTARG}" ;;
 esac
done

if [ -z "$forward" ] || [ -z "$reverse" ] || [ -z "$file1out" ] || [ -z "$configFile" ] || [ -z "$installFolder" ]
then
echo ""
echo "ERROR ($name): All mandatory options must be set"
printHelp
exit 1
fi

# first mawk creates a file in bed format including all necessary properties of spikes
# than the result will be intersected with the padded version of our own created spike file
# spike in        start   end    base    total  current   spike in       start    end  meth_state strand
# Q1hmC_upper     107     108     C       9       9       Q1hmC_upper     107     108     C       +
# Q1hmC_upper     109     110     C       9       1       Q1hmC_upper     109     110     5-mC    +
# Q1hmC_upper     112     113     C       9       9       Q1hmC_upper     112     113     C       +
# Q1hmC_upper     117     118     C       9       0       Q1hmC_upper     117     118     5-mC    +
# Q1hmC_upper     104     105     G       8       1       Q1hmC_upper     104     105     5-mG    -
# Q1hmC_upper     110     111     G       8       8       Q1hmC_upper     110     111     G       -
# calculate average conversion and frequency
annotationFile="$(python3 $installFolder/wrappers/getArgumentfromConfig.py $configFile System annotation)"
padSize="$(python3 $installFolder/wrappers/getArgumentfromConfig.py $configFile Analysis padsize)"


cat <(
      mawk -v OFS='\t' '
          $3=="C" && $4!=0 {print $1,$2,($2+1),$3,$4,gsub("[Tt]","x",$5)}
        ' $forward
    ) <(
        mawk -v OFS='\t' '
           $3=="G" && $4!=0 {print $1,$2,($2+1),$3,$4,gsub("[Aa]","x",$5)}
         ' $reverse
    ) \
    | $bedtools intersect -wb -a - -b <(
     awk -vpad=$padSize -vOFS='\t' '{$2+=pad; $3+=pad; print}' $annotationFile
    ) \
    | mawk -vOFS='\t' '
      {
        tot[$10]+=$6
        cov[$10]+=$5
        count[$10]++
      }
      END{
        for(s in tot) {
          print "spikein_"s"_numbofsites",count[s]"\n" "spikein_"s"_avgcov",cov[s]/count[s]"\n" "spikein_"s"_avgfreq",tot[s]/cov[s]
         }
      }
    ' \
> $file1out


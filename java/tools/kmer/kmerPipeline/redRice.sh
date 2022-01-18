echo `date` "Starting ..."

FASTQA1=${1}
FASTQA2=${2}
FASTQB1=${3}
FASTQB2=${4}

OUTPREFIXA=${5}
OUTPREFIXB=${6}


CORES=20 #number of cores to use for blast searches
KMERSIZE=31
MAXCUTOFF=100 # maximum value for a cutoff from the histogram
LMERSIZE=25 # l-mer for seed generation
MINSIZE=21 # minimum number of kmers to make up a seed (minlength= MINSIZE+KMERSIZE-1)
HIGHBOUNDARY=10000 # kmers with a total coverage of this many kmers are discarded
TOLERANCE=10 #tolerance in seed pairing

JAVA=java
SCRIPTS=/projects/dep_coupland/grp_nordstrom/projects/kmerMutants/scripts

BLASTUTILS=~/jars/blastUtils.jar
KMERUTILS=~/jars/kmerUtils.jar
FASTAUTILS=~/jars/fastaUtil.jar
FASTQUTILS=~/jars/fastqUtils.jar
CSVUTILS=~/jars/csvUtils.jar
SHOREMAPANNOTATE=~/shoreMap/NEW/SHOREmap_annotate.pl

BLASTDBGENOME="/projects/dep_coupland/grp_nordstrom/data/Rice/japonica/shore_rice/Rice.fa.shore"
GENOMEGFF="/projects/dep_coupland/grp_nordstrom/data/Rice/japonica/Rice.gff3"


# 2.1. find cutoff values
#    < kmer dist
#    > cutoff values

#A
CUTOFFA=$(cat ${OUTPREFIXA}.kmers.hist.csv | awk -v MAXCUTOFF=${MAXCUTOFF} -f ${SCRIPTS}/histCutoffPipe.awk)

#B
CUTOFFB=$(cat ${OUTPREFIXB}.kmers.hist.csv | awk -v MAXCUTOFF=${MAXCUTOFF} -f ${SCRIPTS}/histCutoffPipe.awk)

echo `date` "Extracting reads ..."

bsub -K -g /kn/kmerCount "bash ${SCRIPTS}/mergeCountFilesMarkFilter.sh ${SCRIPTS} ${OUTPREFIXA}.kmers.gz ${OUTPREFIXB}.kmers.gz ${CUTOFFA} ${CUTOFFA} ${CUTOFFB} ${CUTOFFB} ${OUTPREFIXA}_unique.kmerDiff ${OUTPREFIXB}_unique.kmerDiff"

# 4. generate seeds and remove highly repetitive sequences
#  < unique kmers, lmer-length (% of kmer-length??), minSize, high boundary
#  > seeds
#

echo `date` "Generating seeds A ..."

#A
bsub -K -R "rusage[mem=96000]" -g /kn "${JAVA} -Xmx92G -Xms92G -XX:-UseGCOverheadLimit -jar ${KMERUTILS} generateSeeds3 ${OUTPREFIXA}_unique.kmerDiff ${OUTPREFIXB}_unique.kmerDiff 0 ${MINSIZE} ${LMERSIZE} ${OUTPREFIXA}_seedsRep > ${OUTPREFIXA}_seedsRep_endCondition.csv"
bsub -K -g /kn "cat ${OUTPREFIXA}_seedsRep_long.txt | awk -f ${SCRIPTS}/longToInfo.awk > ${OUTPREFIXA}_seedsRep_longTable.csv"
cat ${OUTPREFIXA}_seedsRep_longTable.csv| awk -v hb=${HIGHBOUNDARY} '$3<hb && $4<hb {print $1}' > ${OUTPREFIXA}_seeds.gi
${JAVA} -jar ${FASTAUTILS} sub ${OUTPREFIXA}_seedsRep.fa ${OUTPREFIXA}_seeds.gi > ${OUTPREFIXA}_seeds.fa

echo `date` "Generating seeds B ..."

#B
bsub -K -R "rusage[mem=96000]" -g /kn "${JAVA} -Xmx92G -Xms92G -XX:-UseGCOverheadLimit -jar ${KMERUTILS} generateSeeds3 ${OUTPREFIXB}_unique.kmerDiff ${OUTPREFIXA}_unique.kmerDiff 0 ${MINSIZE} ${LMERSIZE} ${OUTPREFIXB}_seedsRep > ${OUTPREFIXB}_seedsRep_endCondition.csv"
bsub -K -g /kn "cat ${OUTPREFIXB}_seedsRep_long.txt | awk -f ${SCRIPTS}/longToInfo.awk > ${OUTPREFIXB}_seedsRep_longTable.csv"
cat ${OUTPREFIXB}_seedsRep_longTable.csv| awk -v hb=${HIGHBOUNDARY} '$3<hb && $4<hb {print $1}' > ${OUTPREFIXB}_seeds.gi
${JAVA} -jar ${FASTAUTILS} sub ${OUTPREFIXB}_seedsRep.fa ${OUTPREFIXB}_seeds.gi > ${OUTPREFIXB}_seeds.fa

if [ $( grep -c '>' ${OUTPREFIXA}_seeds.fa ) -lt $( grep -c '>' ${OUTPREFIXB}_seeds.fa ) ]
then
 #start with A
 PREFIX1=$OUTPREFIXA
 PREFIX2=$OUTPREFIXB
else
 #start with B
 PREFIX1=$OUTPREFIXB
 PREFIX2=$OUTPREFIXA
fi

echo `date` "Cleaning small set ..."

makeblastdb -in ${PREFIX1}_seeds.fa -dbtype nucl
bsub -K -R "rusage[cpu=${CORES}]" -g /kn/blast "blastn -num_threads ${CORES} -db ${PREFIX1}_seeds.fa -query ${PREFIX1}_seeds.fa -dust no -outfmt '6 std qlen slen'  |${JAVA} -jar ${KMERUTILS} blastRevCompFilterPipe ${PREFIX1}_seeds.fa > ${PREFIX1}_seeds_filtered.fa"

echo `date` "Reducing large set ..."

makeblastdb -in ${PREFIX2}_seeds.fa -dbtype nucl
#prefilter second data set
bsub -K -R "rusage[cpu=${CORES}]" -g /kn/blast "blastn -num_threads ${CORES} -db ${PREFIX2}_seeds.fa -query ${PREFIX1}_seeds_filtered.fa -task blastn -outfmt 6| cut -f 2| sort -g |uniq > ${PREFIX2}_seeds_toUse.gi"
${JAVA} -jar ${FASTAUTILS} sub ${PREFIX2}_seeds.fa ${PREFIX2}_seeds_toUse.gi > ${PREFIX2}_seeds_toUse.fa

echo `date` "Cleaning reduced large set ..."

makeblastdb -in ${PREFIX2}_seeds_toUse.fa -dbtype nucl
bsub -K -R "rusage[cpu=${CORES}]" -g /kn/blast "blastn -num_threads ${CORES} -db ${PREFIX2}_seeds_toUse.fa -query ${PREFIX2}_seeds_toUse.fa -dust no -outfmt '6 std qlen slen'  |${JAVA} -jar ${KMERUTILS} blastRevCompFilterPipe ${PREFIX2}_seeds_toUse.fa > ${PREFIX2}_seeds_filtered.fa"

# 5. contrast seeds (BLASTN -task blastn) remove sequences with multiple hits
#

echo `date` "Contrasting seeds A ..."
 
#A
makeblastdb -in ${OUTPREFIXB}_seeds_filtered.fa -dbtype nucl
bsub -K -R "rusage[cpu=${CORES}]" -g /kn/blast "blastn -num_threads ${CORES} -task blastn -dust no -outfmt '6 std qseq sseq' -db ${OUTPREFIXB}_seeds_filtered.fa -query ${OUTPREFIXA}_seeds_filtered.fa > ${OUTPREFIXA}_seeds_filtered_inOther_blast.csv"

echo `date` "Contrasting seeds B ..."

#B
makeblastdb -in ${OUTPREFIXA}_seeds_filtered.fa -dbtype nucl
bsub -K -R "rusage[cpu=${CORES}]" -g /kn/blast "blastn -num_threads ${CORES} -task blastn -dust no -outfmt '6 std qseq sseq' -db ${OUTPREFIXA}_seeds_filtered.fa -query ${OUTPREFIXB}_seeds_filtered.fa > ${OUTPREFIXB}_seeds_filtered_inOther_blast.csv"

# 5.1. SNPs
#    < seeds, tolerance
#    > candidate seed pairs, annotation file
#

#A

cat ${OUTPREFIXA}_seeds_filtered_inOther_blast.csv | ${JAVA} -jar ${KMERUTILS} blastSNPFilterPipe ${KMERSIZE} ${TOLERANCE} > ${OUTPREFIXA}_SNP_blastRep.csv

#B

cat ${OUTPREFIXB}_seeds_filtered_inOther_blast.csv | ${JAVA} -jar ${KMERUTILS} blastSNPFilterPipe ${KMERSIZE} ${TOLERANCE} > ${OUTPREFIXB}_SNP_blastRep.csv

${JAVA} -jar ${KMERUTILS} blastRemoveRepetitive ${OUTPREFIXA}_SNP_blastRep.csv ${OUTPREFIXB}_SNP_blastRep.csv ${OUTPREFIXA}_SNP_blast.csv ${OUTPREFIXB}_SNP_blast.csv


# 5.2. InDels
#    < seeds, ???
#    > candidate seed pairs???, annotation file
#

# 5.3. Merge
#    < annotations files
#    > merged annotations file
#

#A
cat ${OUTPREFIXA}_SNP_blast.csv > ${OUTPREFIXA}_annotated_blast.csv

#B
cat ${OUTPREFIXB}_SNP_blast.csv > ${OUTPREFIXB}_annotated_blast.csv

# 5.4. Extract fasta file with all seeds
#    < seeds fasta file, all paired seeds
#    > candidate seeds
#

#A
cat ${OUTPREFIXA}_SNP_blast.csv |cut -f 1 |sort -g|uniq > ${OUTPREFIXA}_candidates.gi
${JAVA} -jar ${FASTAUTILS} sub ${OUTPREFIXA}_seeds.fa ${OUTPREFIXA}_candidates.gi > ${OUTPREFIXA}_candidates.fa

#B
cat ${OUTPREFIXB}_SNP_blast.csv |cut -f 1 |sort -g|uniq > ${OUTPREFIXB}_candidates.gi
${JAVA} -jar ${FASTAUTILS} sub ${OUTPREFIXB}_seeds.fa ${OUTPREFIXB}_candidates.gi > ${OUTPREFIXB}_candidates.fa


# 6. extend candidate seeds
#

# 6.1. kmerize candidate seeds
#    < candidate seeds
#    > candidate seeds kmers
#

#A
${JAVA} -jar ${FASTAUTILS} kmerize ${OUTPREFIXA}_candidates.fa ${KMERSIZE} |sort|uniq >  ${OUTPREFIXA}_candidates.kmers

#B
${JAVA} -jar ${FASTAUTILS} kmerize ${OUTPREFIXB}_candidates.fa ${KMERSIZE} |sort|uniq >  ${OUTPREFIXB}_candidates.kmers

# 6.2. extract reads
#    < fastq files, candidate seeds kmers
#    > subseted read files
#

echo `date` "Extracting reads A ..."

#A
bsub -K -g /kn/kmerCount "${JAVA} -jar ${FASTQUTILS} kmerCoverage_extractCoveredZ ${FASTQA1} ${FASTQA2} ${OUTPREFIXA}_candidates.kmers 0 ${OUTPREFIXA}_candidates_reads"

echo `date` "Extracting reads B ..."

#B
bsub -K -g /kn/kmerCount "${JAVA} -jar ${FASTQUTILS} kmerCoverage_extractCoveredZ ${FASTQB1} ${FASTQB2} ${OUTPREFIXB}_candidates.kmers 0 ${OUTPREFIXB}_candidates_reads"

# 6.3. extend each seed
#    < candidate seed, subseted read files
#    > extended seed
#

echo `date` "Assembly A ..."

#A
mkdir ${OUTPREFIXA}_candidates
${JAVA} -jar ${FASTAUTILS} split ${OUTPREFIXA}_candidates.fa ${OUTPREFIXA}_candidates/seq 1
for s in `ls ${OUTPREFIXA}_candidates/seq*.fa |grep -v extended`
do
 bash ${SCRIPTS}/assemble.sh ${s} ${OUTPREFIXA}_candidates_reads_1.fastq ${OUTPREFIXA}_candidates_reads_2.fastq
done
cat ${OUTPREFIXA}_candidates/seq*extended.fa > ${OUTPREFIXA}_candidates_ext.fa

${JAVA} -jar ${FASTAUTILS} lengths ${OUTPREFIXA}_candidates_ext.fa > ${OUTPREFIXA}_candidates_ext.lengths

echo `date` "Assembly B ..."

#B
mkdir ${OUTPREFIXB}_candidates
${JAVA} -jar ${FASTAUTILS} split ${OUTPREFIXB}_candidates.fa ${OUTPREFIXB}_candidates/seq 1
for s in `ls ${OUTPREFIXB}_candidates/seq*.fa |grep -v extended`
do
 bash ${SCRIPTS}/assemble.sh ${s} ${OUTPREFIXB}_candidates_reads_1.fastq ${OUTPREFIXB}_candidates_reads_2.fastq
done
cat ${OUTPREFIXB}_candidates/seq*extended.fa > ${OUTPREFIXB}_candidates_ext.fa

${JAVA} -jar ${FASTAUTILS} lengths ${OUTPREFIXB}_candidates_ext.fa > ${OUTPREFIXB}_candidates_ext.lengths

# 7. annotate position in extended seed
#  < seeds, extended seed, seed annotation file
#  > extended seed annotation file
#

echo `date` "Annotation ..."

#A
makeblastdb -in ${OUTPREFIXA}_candidates_ext.fa -dbtype nucl
bsub -K -R "rusage[cpu=${CORES}]" -g /kn/blast "blastn -num_threads ${CORES} -query ${OUTPREFIXA}_candidates.fa -db ${OUTPREFIXA}_candidates_ext.fa -dust no -outfmt '6 std qlen slen' | ${JAVA} -jar ${KMERUTILS} blastTransferAnnotationPipe ${OUTPREFIXA}_annotated_blast.csv | ${JAVA} -jar ${KMERUTILS} blastMutationEMSAnnotationPipe > ${OUTPREFIXA}_candidates_ext_annotation.csv"

#B
makeblastdb -in ${OUTPREFIXB}_candidates_ext.fa -dbtype nucl
bsub -K -R "rusage[cpu=${CORES}]" -g /kn/blast "blastn -num_threads ${CORES} -query ${OUTPREFIXB}_candidates.fa -db ${OUTPREFIXB}_candidates_ext.fa -dust no -outfmt '6 std qlen slen' | ${JAVA} -jar ${KMERUTILS} blastTransferAnnotationPipe ${OUTPREFIXB}_annotated_blast.csv | ${JAVA} -jar ${KMERUTILS} blastMutationEMSAnnotationPipe > ${OUTPREFIXB}_candidates_ext_annotation.csv"

# 8. functional annotation, by comparison to protein data set
#  < extended seeds, extended seed annotation file, protein data set
#  > extended seed functional annotation
#

#A
bsub -K -R "rusage[cpu=${CORES}]" -g /kn/blast "blastn -num_threads ${CORES} -task blastn -db ${BLASTDBGENOME} -query ${OUTPREFIXA}_candidates_ext.fa -dust no -outfmt '6 std' -evalue 0.01 -out ${OUTPREFIXA}_candidates_ext_inDBgenome_blast.csv"
${JAVA} -jar ${KMERUTILS} blastMutatedHits ${OUTPREFIXA}_candidates_ext_inDBgenome_blast.csv ${OUTPREFIXA}_candidates_ext_annotation.csv 2 | ${JAVA} -jar ${BLASTUTILS} topResultPipe 1 | ${JAVA} -jar ${KMERUTILS} blastTransferAnnotationToChrPipe ${OUTPREFIXA}_candidates_ext_annotation.csv 1 |sort -k 3,3 -g |sort -k 2,2 -s > ${OUTPREFIXA}_candidates_ext_genomeAnnotation.csv
cat ${OUTPREFIXA}_candidates_ext_genomeAnnotation.csv | awk '{print $2"_"$3"\t"$1} ' > ${OUTPREFIXA}_candidates_ext_genomeAnnotation.trans

cat ${OUTPREFIXA}_candidates_ext_genomeAnnotation.csv |cut -f 1-3 > ${OUTPREFIXA}_candidates_ext_mutated.csv

${JAVA} -jar ${KMERUTILS} getGFF3lines ${OUTPREFIXA}_candidates_ext_mutated.csv $GENOMEGFF |awk 'BEGIN{ h["CDS"]=1; h["exon"]=4; h["five_prime_UTR"]=2; h["gene"]=6; h["mRNA"]=5; h["three_prime_UTR"]=3;} { print $1"\t"h[$4]"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' |sort -g -k1,1 -k2,2 -k11,11 |sed s/';'/'\t'/| awk 'BEGIN {key=""} $1!=key {print $1"\t"$5"\t"substr($11,4,14)"\t"substr($11,4);key=$1}' > ${OUTPREFIXA}_candidates_ext_mutatedGFF.csv 

#B
bsub -K -R "rusage[cpu=${CORES}]" -g /kn/blast "blastn -num_threads ${CORES} -task blastn -db ${BLASTDBGENOME} -query ${OUTPREFIXB}_candidates_ext.fa -dust no -outfmt '6 std' -evalue 0.01 -out ${OUTPREFIXB}_candidates_ext_inDBgenome_blast.csv"
${JAVA} -jar ${KMERUTILS} blastMutatedHits ${OUTPREFIXB}_candidates_ext_inDBgenome_blast.csv ${OUTPREFIXB}_candidates_ext_annotation.csv 2 | ${JAVA} -jar ${BLASTUTILS} topResultPipe 1 | ${JAVA} -jar ${KMERUTILS} blastTransferAnnotationToChrPipe ${OUTPREFIXB}_candidates_ext_annotation.csv 1 |sort -k 3,3 -g |sort -k 2,2 -s > ${OUTPREFIXB}_candidates_ext_genomeAnnotation.csv
cat ${OUTPREFIXB}_candidates_ext_genomeAnnotation.csv | awk '{print $2"_"$3"\t"$1} ' > ${OUTPREFIXB}_candidates_ext_genomeAnnotation.trans

cat ${OUTPREFIXB}_candidates_ext_genomeAnnotation.csv |cut -f 1-3 > ${OUTPREFIXB}_candidates_ext_mutated.csv

${JAVA} -jar ${KMERUTILS} getGFF3lines ${OUTPREFIXB}_candidates_ext_mutated.csv $GENOMEGFF |awk 'BEGIN{ h["CDS"]=1; h["exon"]=4; h["five_prime_UTR"]=2; h["gene"]=6; h["mRNA"]=5; h["three_prime_UTR"]=3;} { print $1"\t"h[$4]"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' |sort -g -k1,1 -k2,2 -k11,11 |sed s/';'/'\t'/| awk 'BEGIN {key=""} $1!=key {print $1"\t"$5"\t"substr($11,4,14)"\t"substr($11,4);key=$1}' > ${OUTPREFIXB}_candidates_ext_mutatedGFF.csv

# 9. summarize and report in a human readable format
#  < 
#  > summary table

echo `date` "Merging results ..."

#A
cat ${OUTPREFIXA}_seedsRep_longTable.csv |grep -f <(cat ${OUTPREFIXA}_candidates.gi | awk '{print "^"$1"\t"}') |awk '{print $1"\t"$3"\t"$4}' > ${OUTPREFIXA}_candidates_counts.csv
${JAVA} -jar ${CSVUTILS} leftjoin ${OUTPREFIXA}_candidates.gi ${OUTPREFIXA}_candidates_ext.lengths ${OUTPREFIXA}_candidates_ext_annotation.csv ${OUTPREFIXA}_candidates_counts.csv ${OUTPREFIXA}_candidates_ext_mutated.csv ${OUTPREFIXA}_candidates_ext_mutatedGFF.csv 0 |sort -g|grep -v key > ${OUTPREFIXA}_raw
cat ${OUTPREFIXA}_raw | awk 'BEGIN{print "A.id\tA.length\tA.mutation\tA.direction\tA.countA\tA.countB\tA.chr\tA.pos\tA.type\tA.gene\tA.geneId"} {print $0}' > ${OUTPREFIXA}_summary.csv

#B
cat ${OUTPREFIXB}_seedsRep_longTable.csv |grep -f <(cat ${OUTPREFIXB}_candidates.gi | awk '{print "^"$1"\t"}') |awk '{print $1"\t"$3"\t"$4}' > ${OUTPREFIXB}_candidates_counts.csv
${JAVA} -jar ${CSVUTILS} leftjoin ${OUTPREFIXB}_candidates.gi ${OUTPREFIXB}_candidates_ext.lengths ${OUTPREFIXB}_candidates_ext_annotation.csv ${OUTPREFIXB}_candidates_counts.csv ${OUTPREFIXB}_candidates_ext_mutated.csv ${OUTPREFIXB}_candidates_ext_mutatedGFF.csv 0 |sort -g |grep -v key > ${OUTPREFIXB}_raw
cat ${OUTPREFIXB}_raw |awk 'BEGIN{print "B.id\tB.length\tB.mutation\tB.direction\tB.countA\tB.countB\tB.chr\tB.pos\tB.type\tB.gene\tB.geneId"} {print $0}' > ${OUTPREFIXB}_summary.csv

#merge
cat ${OUTPREFIXA}_annotated_blast.csv |cut -f 1-12 > ${OUTPREFIXA}_tmpBlast
cat ${OUTPREFIXB}_raw | awk '{print "<AB>\t"$0}' > ${OUTPREFIXB}_raw2
${JAVA} -jar ${CSVUTILS} leftjoin ${OUTPREFIXA}_tmpBlast ${OUTPREFIXA}_raw 0 |grep -v key > ${OUTPREFIXA}_raw2
${JAVA} -jar ${CSVUTILS} leftjoin ${OUTPREFIXA}_raw2 ${OUTPREFIXB}_raw2 1 |grep -v key|awk 'BEGIN {OFS="\t"} {tmp=$1;$1=$2;$2=tmp;print $0}'|sort -g |awk 'BEGIN{print "A.id\tB.id\tidentity\tseed_alignment_length\tmismatches\tgaps\tA.start\tA.end\tB.start\tB.end\tEvalue\tscore\tA.length\tA.mutation\tA.direction\tA.countA\tA.countB\tA.chr\tA.pos\tA.type\tA.gene\tA.geneId\t<AB>\tB.length\tB.mutation\tB.direction\tB.countA\tB.countB\tB.chr\tB.pos\tB.type\tB.gene\tB.geneId\tmirror count\tsupport count\ttolerance"} {print $0"\t"($17+$27)"\t"($16+$28)"\t"($4-length($15)-60)}' > ${OUTPREFIXA}_summaryBOTH.csv

echo `date` "Done!"








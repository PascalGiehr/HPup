# the pipeline

echo `initializing ...`

#parameters



#How to allow multiple fastq files per sample??????
#Sample A fastq files
#FASTQA1="../data/pep1-1/530_216A_ATTCCT_L001_R1_0*fastq.gz"
#FASTQA2="../data/pep1-1/530_216A_ATTCCT_L001_R2_0*fastq.gz"
#Sample B fastq files
#FASTQB1="../pep2BC1F2_396___pep1-2_397/397_ACTGAT_L002_R1_0*fastq.gz"
#FASTQB2="../pep2BC1F2_396___pep1-2_397/397_ACTGAT_L002_R2_0*fastq.gz"

#OUTPREFIXA=""
#OUTPREFIXB=""



#BLASTDB="/projects/dep_coupland/grp_nordstrom/data/Athal/TAIR10/TAIR10_cdna_20110103_representative_gene_model"

#Athaliana
#BLASTDBGENOME="/projects/dep_coupland/grp_nordstrom/data/Athal/TAIR9/shore/genome.fas.shore"
#GENOMEGFF="/projects/dep_coupland/grp_nordstrom/data/Athal/TAIR10/TAIR10_GFF3_genes.gff"

#Rice
#BLASTDBGENOME="/projects/dep_coupland/grp_nordstrom/data/Rice/japonica/shore_rice/Rice.fa.shore"
#GENOMEGFF="/projects/dep_coupland/grp_nordstrom/data/Rice/japonica/Rice.gff3"

CORES=10 #number of cores to use for blast searches
KMERSIZE=31
MAXCUTOFF=100 # maximum value for a cutoff from the histogram
LMERSIZE=25 # l-mer for seed generation
MINSIZE=21 # minimum number of kmers to make up a seed (minlength= MINSIZE+KMERSIZE-1)
HIGHBOUNDARY=10000 # kmers with a total coverage of this many kmers are discarded
TOLERANCE=10 #tolerance in seed pairing

JAVA=java
SCRIPTS=/projects/dep_coupland/grp_nordstrom/projects/kmerMutants/scripts
#BLASTUTILS=${SCRIPTS}/jars/blastUtils.jar
#KMERUTILS=${SCRIPTS}/jars/kmerUtils.jar
#FASTAUTILS=${SCRIPTS}/jars/fastaUtil.jar
#FASTQUTILS=${SCRIPTS}/jars/fastqUtils.jar
#CSVUTILS=${SCRIPTS}/jars/csvUtils.jar
#SHOREMAPANNOTATE=${SCRIPTS}/SHOREmap_annotate.pl

BLASTUTILS=~/jars/blastUtils.jar
KMERUTILS=~/jars/kmerUtils.jar
FASTAUTILS=~/jars/fastaUtil.jar
FASTQUTILS=~/jars/fastqUtils.jar
CSVUTILS=~/jars/csvUtils.jar
SHOREMAPANNOTATE=~/shoreMap/NEW/SHOREmap_annotate.pl


# 1. kmerize
#  < fastq-files, kmer-length, outprefix
#  > kmer-distribution
#


#A
#mkdir ${OUTPREFIXA}_kmers
#COUNTER=1
#for s in `ls ${FASTQA1}`
#do
# ${JAVA} -jar ${FASTQUTILS} kmerize ${s} ${KMERSIZE} |bash ${SCRIPTS}/sortCountPipe.sh | gzip -c > ${OUTPREFIXA}_kmers/${COUNTER}_1.kmers.gz
# ((COUNTER+=1))
#done

#COUNTER=1
#for s in `ls ${FASTQA2}`
#do
# ${JAVA} -jar ${FASTQUTILS} kmerize ${s} ${KMERSIZE} |bash ${SCRIPTS}/sortCountPipe.sh | gzip -c > ${OUTPREFIXA}_kmers/${COUNTER}_2.kmers.gz
# ((COUNTER+=1))
#done

# merge
#CMD="sort -k2 -m"
#COUNTER=1

#for s in `ls ${OUTPREFIXA}_kmers/*kmers.gz`
#do
# CMD=$CMD" <(zcat "${s}")"
# if [ $(($COUNTER%10)) -eq 0 ]
# then
#  CMD=$CMD" |sort -k2 -m -"
# fi
# ((COUNTER+=1))
#done

#echo $CMD > ${OUTPREFIXA}_merge.sh

#bash ${OUTPREFIXA}_merge.sh | bash ${SCRIPTS}/addPipe.sh |gzip -c > ${OUTPREFIXA}.kmers.gz

#rm -rf ${OUTPREFIXA}_kmers

#B
#mkdir ${OUTPREFIXB}_kmers
#COUNTER=1
#for s in `ls ${FASTQB1}`
#do
# ${JAVA} -jar ${FASTQUTILS} kmerize ${s} ${KMERSIZE} |bash ${SCRIPTS}/sortCountPipe.sh | gzip -c > ${OUTPREFIXB}_kmers/${COUNTER}_1.kmers.gz
# ((COUNTER+=1))
#done

#COUNTER=1
#for s in `ls ${FASTQB2}`
#do
# ${JAVA} -jar ${FASTQUTILS} kmerize ${s} ${KMERSIZE} |bash ${SCRIPTS}/sortCountPipe.sh | gzip -c > ${OUTPREFIXB}_kmers/${COUNTER}_2.kmers.gz
# ((COUNTER+=1))
#done

# merge
#CMD="sort -k2 -m"
#COUNTER=1

#for s in `ls ${OUTPREFIXB}_kmers/*kmers.gz`
#do
# CMD=$CMD" <(zcat "${s}")"
# if [ $(($COUNTER%10)) -eq 0 ]
# then
#  CMD=$CMD" |sort -k2 -m -"
# fi
# ((COUNTER+=1))
#done

#echo $CMD > ${OUTPREFIXB}_merge.sh

#bash ${OUTPREFIXB}_merge.sh | bash ${SCRIPTS}/addPipe.sh |gzip -c > ${OUTPREFIXB}.kmers.gz

#rm -rf ${OUTPREFIXB}_kmers

#JELLYFISH

#A

mkdir ${OUTPREFIXA}_kmers

echo `date` "Counting kmers A ..."

echo "jellyfish count -C -o ${OUTPREFIXA}_kmers/tmp -m ${KMERSIZE} -t ${CORES} -s 10G <( zcat ${FASTQA1} ${FASTQA2} )" > ${OUTPREFIXA}_kmerCount.sh

bsub -K -g /kn/jellyfish -R "rusage[cpu=${CORES}]" -R "rusage[mem=50000]" "bash ${OUTPREFIXA}_kmerCount.sh"

#if the hash size was too small... do a merge
COUNT=$(ls ${OUTPREFIXA}_kmers/tmp* |wc -l)

echo `date` "merging ${COUNT} datasets A ..."

if [ $COUNT -eq 1 ]
then
 mv ${OUTPREFIXA}_kmers/tmp_0 ${OUTPREFIXA}_kmers_jellyfish
else
 jellyfish merge -o ${OUTPREFIXA}_kmers_jellyfish ${OUTPREFIXA}_kmers/tmp*
fi
rm -rf ${OUTPREFIXA}_kmers

echo `date` "Calculatng hist A ..."

bsub -K -g /kn/jellyfish -R "rusage[cpu=${CORES}]" "jellyfish histo -f -o ${OUTPREFIXA}.kmers.hist.csv -t ${CORES} ${OUTPREFIXA}_kmers_jellyfish"
awk '{print $2"\t"$1}' ${OUTPREFIXA}.kmers.hist.csv > ${OUTPREFIXA}_tmp
mv ${OUTPREFIXA}_tmp ${OUTPREFIXA}.kmers.hist.csv

#B

mkdir ${OUTPREFIXB}_kmers

echo `date` "Counting kmers B ..."

echo "jellyfish count -C -o ${OUTPREFIXB}_kmers/tmp -m ${KMERSIZE} -t ${CORES} -s 10G <( zcat ${FASTQB1} ${FASTQB2} )" > ${OUTPREFIXB}_kmerCount.sh

bsub -K -g /kn/jellyfish -R "rusage[cpu=${CORES}]" -R "rusage[mem=50000]" "bash ${OUTPREFIXB}_kmerCount.sh"

#if the hash size was too small... do a merge
COUNT=$(ls ${OUTPREFIXB}_kmers/tmp* |wc -l)

echo `date` "merging ${COUNT} datasets B ..."

if [ $COUNT -eq 1 ]
then
 mv ${OUTPREFIXB}_kmers/tmp_0 ${OUTPREFIXB}_kmers_jellyfish
else
 jellyfish merge -o ${OUTPREFIXB}_kmers_jellyfish ${OUTPREFIXB}_kmers/tmp*
fi
rm -rf ${OUTPREFIXB}_kmers

echo `date` "Calculatng hist B ..."

bsub -K -g /kn/jellyfish -R "rusage[cpu=${CORES}]" "jellyfish histo -f -o ${OUTPREFIXB}.kmers.hist.csv -t ${CORES} ${OUTPREFIXB}_kmers_jellyfish"
awk '{print $2"\t"$1}' ${OUTPREFIXB}.kmers.hist.csv > ${OUTPREFIXB}_tmp
mv ${OUTPREFIXB}_tmp ${OUTPREFIXB}.kmers.hist.csv


# 2. calculate histogram
#  > cutoff values
#

#A
#zcat ${OUTPREFIXA}.kmers.gz | bash ${SCRIPTS}/histPipe.sh > ${OUTPREFIXA}.kmers.hist.csv

#B
#zcat ${OUTPREFIXB}.kmers.gz | bash ${SCRIPTS}/histPipe.sh > ${OUTPREFIXB}.kmers.hist.csv

# 2.1. find cutoff values
#    < kmer dist
#    > cutoff values

#A
CUTOFFA=$(awk -v MAXCUTOFF=${MAXCUTOFF} -f ${SCRIPTS}/histCutoffPipe.awk ${OUTPREFIXA}.kmers.hist.csv)

echo $CUTOFFA > ${OUTPREFIXA}_cutoff.csv

#B
CUTOFFB=$(awk -v MAXCUTOFF=${MAXCUTOFF} -f ${SCRIPTS}/histCutoffPipe.awk ${OUTPREFIXB}.kmers.hist.csv)

echo $CUTOFFB > ${OUTPREFIXB}_cutoff.csv

# 3. contrast kmer-distributions
#  < cutoff values
#  > unique kmers
#

#bash ${SCRIPTS}/mergeCountFilesMarkFilter.sh ${SCRIPTS} ${OUTPREFIXA}.kmers.gz ${OUTPREFIXB}.kmers.gz ${CUTOFFA} ${CUTOFFA} ${CUTOFFB} ${CUTOFFB} ${OUTPREFIXA}_unique.kmerDiff ${OUTPREFIXB}_unique.kmerDiff

#JELLYFISH
bsub -K -g /kn/kmerCount "bash ${SCRIPTS}/mergeCountFilesMarkFilterJellyfish3.sh ${SCRIPTS} ${OUTPREFIXA}_kmers_jellyfish ${OUTPREFIXB}_kmers_jellyfish ${CUTOFFA} ${CUTOFFA} ${CUTOFFB} ${CUTOFFB} ${OUTPREFIXA}_unique.kmerDiff ${OUTPREFIXB}_unique.kmerDiff"


# 4. generate seeds and remove highly repetitive sequences
#  < unique kmers, lmer-length (% of kmer-length??), minSize, high boundary
#  > seeds
#

echo `date` "Generating seeds A ..."

#A
bsub -K -R "rusage[mem=64000]" -g /kn/seed "${JAVA} -Xmx64G -Xms64G -XX:-UseGCOverheadLimit -jar ${KMERUTILS} generateSeeds3 ${OUTPREFIXA}_unique.kmerDiff ${OUTPREFIXB}_unique.kmerDiff 0 ${MINSIZE} ${LMERSIZE} ${OUTPREFIXA}_seedsRep |gzip -c > ${OUTPREFIXA}_seedsRep_endCondition.csv.gz"
awk -f ${SCRIPTS}/longToInfo.awk ${OUTPREFIXA}_seedsRep_long.txt |gzip -c> ${OUTPREFIXA}_seedsRep_longTable.csv.gz
gunzip -c ${OUTPREFIXA}_seedsRep_longTable.csv.gz| awk -v hb=${HIGHBOUNDARY} '$3<hb && $4<hb {print $1}' > ${OUTPREFIXA}_seeds.gi
${JAVA} -jar ${FASTAUTILS} sub ${OUTPREFIXA}_seedsRep.fa ${OUTPREFIXA}_seeds.gi > ${OUTPREFIXA}_seeds.fa

echo `date` "Generating seeds B ..."

#B
bsub -K -R "rusage[mem=64000]" -g /kn/seed "${JAVA} -Xmx64G -Xms64G -XX:-UseGCOverheadLimit -jar ${KMERUTILS} generateSeeds3 ${OUTPREFIXB}_unique.kmerDiff ${OUTPREFIXA}_unique.kmerDiff 0 ${MINSIZE} ${LMERSIZE} ${OUTPREFIXB}_seedsRep |gzip -c > ${OUTPREFIXB}_seedsRep_endCondition.csv.gz"
awk -f ${SCRIPTS}/longToInfo.awk ${OUTPREFIXB}_seedsRep_long.txt |gzip -c> ${OUTPREFIXB}_seedsRep_longTable.csv.gz
gunzip -c ${OUTPREFIXB}_seedsRep_longTable.csv.gz| awk -v hb=${HIGHBOUNDARY} '$3<hb && $4<hb {print $1}' > ${OUTPREFIXB}_seeds.gi
${JAVA} -jar ${FASTAUTILS} sub ${OUTPREFIXB}_seedsRep.fa ${OUTPREFIXB}_seeds.gi > ${OUTPREFIXB}_seeds.fa


# 4.1. clean reverse complement
#    < seeds
#    > filtered seeds
#

### remove multiple hitting sequences here

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
bsub -K -g /kn/blastScript bash ${SCRIPTS}/blastnWrapper.sh ${SCRIPTS} ${PREFIX1}_seeds.fa ${PREFIX1}_seeds.fa "-dust no -outfmt '6 std qlen slen'" ${PREFIX1}_seeds_inSelfBlast.csv.gz
#blastn -num_threads ${CORES} -db ${PREFIX1}_seeds.fa -query ${PREFIX1}_seeds.fa -dust no -outfmt '6 std qlen slen'  
gunzip -c ${PREFIX1}_seeds_inSelfBlast.csv.gz |${JAVA} -XX:-UseGCOverheadLimit -jar ${KMERUTILS} blastRevCompFilterPipe ${PREFIX1}_seeds.fa > ${PREFIX1}_seeds_filtered.fa

echo `date` "Reducing large set ..."
#prefilter second data set
makeblastdb -in ${PREFIX2}_seeds.fa -dbtype nucl
bsub -K -g /kn/blastScript bash ${SCRIPTS}/blastnWrapper.sh ${SCRIPTS} ${PREFIX1}_seeds_filtered.fa ${PREFIX2}_seeds.fa "-outfmt 6 -task blastn" ${PREFIX2}_seeds_toUse_blast.csv.gz
#blastn -num_threads ${CORES} -outfmt 6 -db ${PREFIX2}_seeds.fa -query ${PREFIX1}_seeds_filtered.fa -task blastn
gunzip -c ${PREFIX2}_seeds_toUse_blast.csv.gz | cut -f 2| sort -gu > ${PREFIX2}_seeds_toUse.gi
${JAVA} -jar ${FASTAUTILS} sub ${PREFIX2}_seeds.fa ${PREFIX2}_seeds_toUse.gi > ${PREFIX2}_seeds_toUse.fa

echo `date` "Cleaning reduced large set ..."

makeblastdb -in ${PREFIX2}_seeds_toUse.fa -dbtype nucl
bsub -K -g /kn/blastScript bash ${SCRIPTS}/blastnWrapper.sh ${SCRIPTS} ${PREFIX2}_seeds_toUse.fa ${PREFIX2}_seeds_toUse.fa "-dust no -outfmt '6 std qlen slen'" ${PREFIX2}_seeds_toUse_inSelfBlast.csv.gz
#blastn -num_threads ${CORES} -db ${PREFIX2}_seeds_toUse.fa -query ${PREFIX2}_seeds_toUse.fa -dust no -outfmt '6 std qlen slen' 
gunzip -c ${PREFIX2}_seeds_toUse_inSelfBlast.csv.gz |${JAVA} -XX:-UseGCOverheadLimit -jar ${KMERUTILS} blastRevCompFilterPipe ${PREFIX2}_seeds_toUse.fa > ${PREFIX2}_seeds_filtered.fa

# 5. contrast seeds (BLASTN -task blastn) remove sequences with multiple hits
#

echo `date` "Contrasting seeds A ..."
 
#A
makeblastdb -in ${OUTPREFIXB}_seeds_filtered.fa -dbtype nucl
bsub -K -g /kn/blastScript bash ${SCRIPTS}/blastnWrapper.sh ${SCRIPTS} ${OUTPREFIXA}_seeds_filtered.fa ${OUTPREFIXB}_seeds_filtered.fa "-task blastn -dust no -outfmt '6 std qseq sseq'" ${OUTPREFIXA}_seeds_filtered_inOther_blast.csv.gz

#blastn -num_threads ${CORES} -task blastn -dust no -outfmt '6 std qseq sseq' -db ${OUTPREFIXB}_seeds_filtered.fa -query ${OUTPREFIXA}_seeds_filtered.fa | gzip -c > ${OUTPREFIXA}_seeds_filtered_inOther_blast.csv.gz

echo `date` "Contrasting seeds B ..."

#B
makeblastdb -in ${OUTPREFIXA}_seeds_filtered.fa -dbtype nucl
bsub -K -g /kn/blastScript bash ${SCRIPTS}/blastnWrapper.sh ${SCRIPTS} ${OUTPREFIXB}_seeds_filtered.fa ${OUTPREFIXA}_seeds_filtered.fa "-task blastn -dust no -outfmt '6 std qseq sseq'" ${OUTPREFIXB}_seeds_filtered_inOther_blast.csv.gz

#blastn -num_threads ${CORES} -task blastn -dust no -outfmt '6 std qseq sseq' -db ${OUTPREFIXA}_seeds_filtered.fa -query ${OUTPREFIXB}_seeds_filtered.fa | gzip -c > ${OUTPREFIXB}_seeds_filtered_inOther_blast.csv.gz

# 5.1. SNPs
#    < seeds, tolerance
#    > candidate seed pairs, annotation file
#

#A

gunzip -c ${OUTPREFIXA}_seeds_filtered_inOther_blast.csv.gz | ${JAVA} -jar ${KMERUTILS} blastSNPFilterPipe ${KMERSIZE} ${TOLERANCE} > ${OUTPREFIXA}_SNP_blastRep.csv

#B

gunzip -c ${OUTPREFIXB}_seeds_filtered_inOther_blast.csv.gz | ${JAVA} -jar ${KMERUTILS} blastSNPFilterPipe ${KMERSIZE} ${TOLERANCE} > ${OUTPREFIXB}_SNP_blastRep.csv

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
cut -f 1 ${OUTPREFIXA}_SNP_blast.csv |sort -gu > ${OUTPREFIXA}_candidates.gi
${JAVA} -jar ${FASTAUTILS} sub ${OUTPREFIXA}_seeds_filtered.fa ${OUTPREFIXA}_candidates.gi > ${OUTPREFIXA}_candidates.fa

#B
cut -f 1 ${OUTPREFIXB}_SNP_blast.csv |sort -gu > ${OUTPREFIXB}_candidates.gi
${JAVA} -jar ${FASTAUTILS} sub ${OUTPREFIXB}_seeds_filtered.fa ${OUTPREFIXB}_candidates.gi > ${OUTPREFIXB}_candidates.fa


# 6. extend candidate seeds
#

# 6.1. kmerize candidate seeds
#    < candidate seeds
#    > candidate seeds kmers
#

#A
${JAVA} -jar ${FASTAUTILS} kmerize ${OUTPREFIXA}_candidates.fa ${KMERSIZE} |sort -u >  ${OUTPREFIXA}_candidates.kmers

#B
${JAVA} -jar ${FASTAUTILS} kmerize ${OUTPREFIXB}_candidates.fa ${KMERSIZE} |sort -u >  ${OUTPREFIXB}_candidates.kmers

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
 bsub -o /dev/null -e /dev/null -q short bash ${SCRIPTS}/assemble.sh ${s} ${OUTPREFIXA}_candidates_reads_1.fastq ${OUTPREFIXA}_candidates_reads_2.fastq
done
cat ${OUTPREFIXA}_candidates/seq*extended.fa > ${OUTPREFIXA}_candidates_ext.fa

${JAVA} -jar ${FASTAUTILS} lengths ${OUTPREFIXA}_candidates_ext.fa > ${OUTPREFIXA}_candidates_ext.lengths

tar czf ${OUTPREFIXA}_candidates_extData.tar.gz ${OUTPREFIXA}_candidates

echo `date` "Assembly B ..."

#B
mkdir ${OUTPREFIXB}_candidates
${JAVA} -jar ${FASTAUTILS} split ${OUTPREFIXB}_candidates.fa ${OUTPREFIXB}_candidates/seq 1
for s in `ls ${OUTPREFIXB}_candidates/seq*.fa |grep -v extended`
do
 bsub -o /dev/null -e /dev/null -q short bash ${SCRIPTS}/assemble.sh ${s} ${OUTPREFIXB}_candidates_reads_1.fastq ${OUTPREFIXB}_candidates_reads_2.fastq
done
cat ${OUTPREFIXB}_candidates/seq*extended.fa > ${OUTPREFIXB}_candidates_ext.fa

${JAVA} -jar ${FASTAUTILS} lengths ${OUTPREFIXB}_candidates_ext.fa > ${OUTPREFIXB}_candidates_ext.lengths

tar czf ${OUTPREFIXB}_candidates_extData.tar.gz ${OUTPREFIXB}_candidates

# 7. annotate position in extended seed
#  < seeds, extended seed, seed annotation file
#  > extended seed annotation file
#

echo `date` "Annotation ..."

#A
makeblastdb -in ${OUTPREFIXA}_candidates_ext.fa -dbtype nucl
blastn -query ${OUTPREFIXA}_candidates.fa -db ${OUTPREFIXA}_candidates_ext.fa -dust no -outfmt '6 std qlen slen' | ${JAVA} -jar ${KMERUTILS} blastTransferAnnotationPipe ${OUTPREFIXA}_annotated_blast.csv | ${JAVA} -jar ${KMERUTILS} blastMutationEMSAnnotationPipe > ${OUTPREFIXA}_candidates_ext_annotation.csv

#B
makeblastdb -in ${OUTPREFIXB}_candidates_ext.fa -dbtype nucl
blastn -query ${OUTPREFIXB}_candidates.fa -db ${OUTPREFIXB}_candidates_ext.fa -dust no -outfmt '6 std qlen slen' | ${JAVA} -jar ${KMERUTILS} blastTransferAnnotationPipe ${OUTPREFIXB}_annotated_blast.csv | ${JAVA} -jar ${KMERUTILS} blastMutationEMSAnnotationPipe > ${OUTPREFIXB}_candidates_ext_annotation.csv

# 8. functional annotation, by comparison to protein data set
#  < extended seeds, extended seed annotation file, protein data set
#  > extended seed functional annotation
#

#A
blastn -num_threads ${CORES} -task blastn -db ${BLASTDBGENOME} -query ${OUTPREFIXA}_candidates_ext.fa -dust no -outfmt '6 std qseq sseq' -evalue 0.01 -out ${OUTPREFIXA}_candidates_ext_inDBgenome_blast.csv
${JAVA} -jar ${KMERUTILS} blastMutatedHits ${OUTPREFIXA}_candidates_ext_inDBgenome_blast.csv ${OUTPREFIXA}_candidates_ext_annotation.csv 2 | ${JAVA} -jar ${BLASTUTILS} topResultPipe 1 | ${JAVA} -jar ${KMERUTILS} blastTransferAnnotationToChrPipe ${OUTPREFIXA}_candidates_ext_annotation.csv 1 |sort -k 3,3 -g |sort -k 2,2 -s > ${OUTPREFIXA}_candidates_ext_genomeAnnotation.csv
awk '{print $2"_"$3"\t"$1} '  ${OUTPREFIXA}_candidates_ext_genomeAnnotation.csv > ${OUTPREFIXA}_candidates_ext_genomeAnnotation.trans

#cut -f 1-3 ${OUTPREFIXA}_candidates_ext_genomeAnnotation.csv > ${OUTPREFIXA}_candidates_ext_mutated.csv

#Rice
#${JAVA} -jar ${KMERUTILS} getGFF3lines ${OUTPREFIXA}_candidates_ext_mutated.csv $GENOMEGFF |awk 'BEGIN{ h["CDS"]=1; h["exon"]=4; h["five_prime_UTR"]=2; h["gene"]=6; h["mRNA"]=5; h["three_prime_UTR"]=3;} { print $1"\t"h[$4]"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' |sort -g -k1,1 -k2,2 -k11,11 |sed s/';'/'\t'/| awk 'BEGIN {key=""} $1!=key {print $1"\t"$5"\t"substr($11,4,14)"\t"substr($11,4);key=$1}' > ${OUTPREFIXA}_candidates_ext_mutatedGFF.csv

rm -f ${OUTPREFIXA}_candidates_ext_mutated_raw.csv
for s in `${JAVA} -jar ${FASTAUTILS} giList ${BLASTDBGENOME}`
do
 perl -I ${SCRIPTS} ${SHOREMAPANNOTATE} --snp ${OUTPREFIXA}_candidates_ext_genomeAnnotation.csv --chrom ${s} --start 0 --end 99999999999 --gff ${GENOMEGFF} --genome ${BLASTDBGENOME}
 cat snp.priority.txt >> ${OUTPREFIXA}_candidates_ext_mutated_raw.csv
done
rm snp.priority.txt

cat ${OUTPREFIXA}_candidates_ext_mutated_raw.csv |awk '{if(NF<11 ||length($11)==0){$11=0} if(NF<14 ||length($14)==0){$14=0} if(NF<15 ||length($15)==0){$15=0} if(NF<16 ||length($16)==0){$16=0} if(NF<17 ||length($17)==0){$17=0} print $1"_"$2"\t"$10"\t"$11"\t"$14"\t"$15"\t"$16"\t"$17}'  > ${OUTPREFIXA}_candidates_ext_mutated_rawExtRed.csv

${JAVA} -jar ${CSVUTILS} leftjoin ${OUTPREFIXA}_candidates_ext_genomeAnnotation.trans ${OUTPREFIXA}_candidates_ext_mutated_rawExtRed.csv 0 |grep -v key |cut -f 2-8|sort -g|uniq >${OUTPREFIXA}_candidates_ext_mutated.csv

#${JAVA} -jar ${KMERUTILS} getGFF3lines ${OUTPREFIXA}_candidates_ext_mutated.csv $GENOMEGFF > ${OUTPREFIXA}_candidates_ext_mutatedGFF.csv

#cat ${OUTPREFIXA}_candidates_ext_genomeAnnotation.csv |cut -f 1-3 > ${OUTPREFIXA}_candidates_ext_mutated.csv

#blastn -num_threads ${CORES} -task blastn -db ${BLASTDB} -query ${OUTPREFIXA}_candidates_ext.fa -dust no -outfmt '6 std' -evalue 0.01 -out ${OUTPREFIXA}_candidates_ext_inDB_blast.csv
#${JAVA} -jar ${KMERUTILS} blastMutatedHits ${OUTPREFIXA}_candidates_ext_inDB_blast.csv ${OUTPREFIXA}_candidates_ext_annotation.csv 2 | ${JAVA} -jar ${BLASTUTILS} topResultPipe 1 > ${OUTPREFIXA}_candidates_ext_inDB_blast_mutated.csv
#cat ${OUTPREFIXA}_candidates_ext_inDB_blast_mutated.csv |cut -f 1,2 > ${OUTPREFIXA}_candidates_ext_mutated.csv


#${JAVA} -jar ${BLASTUTILS} topResultNonOverlapping ${OUTPREFIXA}_candidates_ext_inDB_blast.csv 100 > ${OUTPREFIXA}_candidates_ext_inDB_blast_close.csv

#B
blastn -num_threads ${CORES} -task blastn -db ${BLASTDBGENOME} -query ${OUTPREFIXB}_candidates_ext.fa -dust no -outfmt '6 std qseq sseq' -evalue 0.01 -out ${OUTPREFIXB}_candidates_ext_inDBgenome_blast.csv
${JAVA} -jar ${KMERUTILS} blastMutatedHits ${OUTPREFIXB}_candidates_ext_inDBgenome_blast.csv ${OUTPREFIXB}_candidates_ext_annotation.csv 2 | ${JAVA} -jar ${BLASTUTILS} topResultPipe 1 | ${JAVA} -jar ${KMERUTILS} blastTransferAnnotationToChrPipe ${OUTPREFIXB}_candidates_ext_annotation.csv 1 |sort -k2,2 -k3,3g > ${OUTPREFIXB}_candidates_ext_genomeAnnotation.csv
awk '{print $2"_"$3"\t"$1} ' ${OUTPREFIXB}_candidates_ext_genomeAnnotation.csv > ${OUTPREFIXB}_candidates_ext_genomeAnnotation.trans

#cut -f 1-3 ${OUTPREFIXB}_candidates_ext_genomeAnnotation.csv > ${OUTPREFIXB}_candidates_ext_mutated.csv

#Rice
#${JAVA} -jar ${KMERUTILS} getGFF3lines ${OUTPREFIXB}_candidates_ext_mutated.csv $GENOMEGFF |awk 'BEGIN{ h["CDS"]=1; h["exon"]=4; h["five_prime_UTR"]=2; h["gene"]=6; h["mRNA"]=5; h["three_prime_UTR"]=3;} { print $1"\t"h[$4]"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' |sort -g -k1,1 -k2,2 -k11,11 |sed s/';'/'\t'/| awk 'BEGIN {key=""} $1!=key {print $1"\t"$5"\t"substr($11,4,14)"\t"substr($11,4);key=$1}' > ${OUTPREFIXB}_candidates_ext_mutatedGFF.csv


rm -f ${OUTPREFIXB}_candidates_ext_mutated_raw.csv
for s in `${JAVA} -jar ${FASTAUTILS} giList ${BLASTDBGENOME}`
do
 perl -I ${SCRIPTS} ${SHOREMAPANNOTATE} --snp ${OUTPREFIXB}_candidates_ext_genomeAnnotation.csv --chrom ${s} --start 0 --end 99999999999 --gff ${GENOMEGFF} --genome ${BLASTDBGENOME}
 cat snp.priority.txt >> ${OUTPREFIXB}_candidates_ext_mutated_raw.csv
done
rm snp.priority.txt

cat ${OUTPREFIXB}_candidates_ext_mutated_raw.csv |awk '{if(NF<11 ||length($11)==0){$11=0} if(NF<14 ||length($14)==0){$14=0} if(NF<15 ||length($15)==0){$15=0} if(NF<16 ||length($16)==0){$16=0} if(NF<17 ||length($17)==0){$17=0} print $1"_"$2"\t"$10"\t"$11"\t"$14"\t"$15"\t"$16"\t"$17}'  > ${OUTPREFIXB}_candidates_ext_mutated_rawExtRed.csv

${JAVA} -jar ${CSVUTILS} leftjoin ${OUTPREFIXB}_candidates_ext_genomeAnnotation.trans ${OUTPREFIXB}_candidates_ext_mutated_rawExtRed.csv 0 |grep -v key |cut -f 2-8|sort -g|uniq >${OUTPREFIXB}_candidates_ext_mutated.csv

#${JAVA} -jar ${KMERUTILS} getGFF3lines ${OUTPREFIXB}_candidates_ext_mutated.csv $GENOMEGFF > ${OUTPREFIXB}_candidates_ext_mutatedGFF.csv

#cat ${OUTPREFIXB}_candidates_ext_genomeAnnotation.csv |cut -f 1-3 > ${OUTPREFIXB}_candidates_ext_mutated.csv


#blastn -num_threads ${CORES} -task blastn -db ${BLASTDB} -query ${OUTPREFIXB}_candidates_ext.fa -dust no -outfmt '6 std' -evalue 0.01 -out ${OUTPREFIXB}_candidates_ext_inDB_blast.csv
#${JAVA} -jar ${KMERUTILS} blastMutatedHits ${OUTPREFIXB}_candidates_ext_inDB_blast.csv ${OUTPREFIXB}_candidates_ext_annotation.csv 2 | ${JAVA} -jar ${BLASTUTILS} topResultPipe 1 > ${OUTPREFIXB}_candidates_ext_inDB_blast_mutated.csv
#cat ${OUTPREFIXB}_candidates_ext_inDB_blast_mutated.csv |cut -f 1,2 > ${OUTPREFIXB}_candidates_ext_mutated.csv


#${JAVA} -jar ${BLASTUTILS} topResultNonOverlapping ${OUTPREFIXB}_candidates_ext_inDB_blast.csv 100 > ${OUTPREFIXB}_candidates_ext_inDB_blast_close.csv

# 9. summarize and report in a human readable format
#  < 
#  > summary table

echo `date` "Merging results ..."

#A
join -1 1 -2 1 <(gunzip -c ${OUTPREFIXA}_seedsRep_longTable.csv.gz |grep -Ff <(cut -f 1 ${OUTPREFIXA}_candidates.gi ) |sort -k1,1) <(cut -f 1 ${OUTPREFIXA}_candidates.gi |sort -k1,1) |awk '{print $1"\t"$3"\t"$4}' |sort -k1,1g > ${OUTPREFIXA}_candidates_counts.csv
#gunzip -c ${OUTPREFIXA}_seedsRep_longTable.csv.gz |grep -f <(awk '{print "^"$1"\t"}' ${OUTPREFIXA}_candidates.gi) |awk '{print $1"\t"$3"\t"$4}' > ${OUTPREFIXA}_candidates_counts.csv
${JAVA} -jar ${CSVUTILS} leftjoin ${OUTPREFIXA}_candidates.gi ${OUTPREFIXA}_candidates_ext.lengths ${OUTPREFIXA}_candidates_ext_annotation.csv ${OUTPREFIXA}_candidates_counts.csv ${OUTPREFIXA}_candidates_ext_mutated.csv ${OUTPREFIXA}_candidates_ext_mutatedGFF.csv 0 |sort -g|grep -v key > ${OUTPREFIXA}_raw
#Athal
#${JAVA} -jar ${CSVUTILS} leftjoin ${OUTPREFIXA}_candidates.gi ${OUTPREFIXA}_candidates_ext.lengths ${OUTPREFIXA}_candidates_ext_annotation.csv ${OUTPREFIXA}_candidates_counts.csv ${OUTPREFIXA}_candidates_ext_mutated.csv 0 |sort -g|grep -v key > ${OUTPREFIXA}_raw
awk 'BEGIN{print "A.id\tA.length\tA.mutation\tA.direction\tA.countA\tA.countB\tA.chr\tA.pos\tA.type\tA.gene\tA.geneId"} {print $0}' ${OUTPREFIXA}_raw > ${OUTPREFIXA}_summary.csv

#B
join -1 1 -2 1 <(gunzip -c ${OUTPREFIXB}_seedsRep_longTable.csv.gz |grep -Ff <(cut -f 1 ${OUTPREFIXB}_candidates.gi ) |sort -k1,1) <(cut -f 1 ${OUTPREFIXB}_candidates.gi |sort -k1,1) |awk '{print $1"\t"$3"\t"$4}' |sort -k1,1g > ${OUTPREFIXB}_candidates_counts.csv
#gunzip -c ${OUTPREFIXB}_seedsRep_longTable.csv.gz |grep -f <(awk '{print "^"$1"\t"}' ${OUTPREFIXB}_candidates.gi) |awk '{print $1"\t"$3"\t"$4}' > ${OUTPREFIXB}_candidates_counts.csv
${JAVA} -jar ${CSVUTILS} leftjoin ${OUTPREFIXB}_candidates.gi ${OUTPREFIXB}_candidates_ext.lengths ${OUTPREFIXB}_candidates_ext_annotation.csv ${OUTPREFIXB}_candidates_counts.csv ${OUTPREFIXB}_candidates_ext_mutated.csv ${OUTPREFIXB}_candidates_ext_mutatedGFF.csv 0 |sort -g |grep -v key > ${OUTPREFIXB}_raw
#Athal
#${JAVA} -jar ${CSVUTILS} leftjoin ${OUTPREFIXB}_candidates.gi ${OUTPREFIXB}_candidates_ext.lengths ${OUTPREFIXB}_candidates_ext_annotation.csv ${OUTPREFIXB}_candidates_counts.csv ${OUTPREFIXB}_candidates_ext_mutated.csv 0 |sort -g |grep -v key > ${OUTPREFIXB}_raw
awk 'BEGIN{print "B.id\tB.length\tB.mutation\tB.direction\tB.countA\tB.countB\tB.chr\tB.pos\tB.type\tB.gene\tB.geneId"} {print $0}' ${OUTPREFIXB}_raw > ${OUTPREFIXB}_summary.csv

#merge
cut -f 1-12 ${OUTPREFIXA}_annotated_blast.csv > ${OUTPREFIXA}_tmpBlast
awk '{print "<AB>\t"$0}' ${OUTPREFIXB}_raw > ${OUTPREFIXB}_raw2
${JAVA} -jar ${CSVUTILS} leftjoin ${OUTPREFIXA}_tmpBlast ${OUTPREFIXA}_raw 0 |grep -v key > ${OUTPREFIXA}_raw2
#${JAVA} -jar ${CSVUTILS} leftjoin ${OUTPREFIXA}_raw2 ${OUTPREFIXB}_raw2 1 |grep -v key|awk 'BEGIN {OFS="\t"} {tmp=$1;$1=$2;$2=tmp;print $0}'|sort -g |awk 'BEGIN{print "A.id\tB.id\tidentity\tseed_alignment_length\tmismatches\tgaps\tA.start\tA.end\tB.start\tB.end\tEvalue\tscore\tA.length\tA.mutation\tA.direction\tA.countA\tA.countB\tA.chr\tA.pos\tA.type\tA.gene\tA.geneId\t<AB>\tB.length\tB.mutation\tB.direction\tB.countA\tB.countB\tB.chr\tB.pos\tB.type\tB.gene\tB.geneId\tmirror count\tsupport count\ttolerance"} {print $0"\t"($17+$27)"\t"($16+$28)"\t"($4-length($15)-60)}' > ${OUTPREFIXA}_summaryBOTH.csv

#Athal
${JAVA} -jar ${CSVUTILS} leftjoin ${OUTPREFIXA}_raw2 ${OUTPREFIXB}_raw2 1 |grep -v key|awk 'BEGIN {OFS="\t"} {tmp=$1;$1=$2;$2=tmp;print $0}'|sort -g |awk 'BEGIN{print "A.id\tB.id\tidentity\tseed_alignment_length\tmismatches\tgaps\tA.start\tA.end\tB.start\tB.end\tEvalue\tscore\tA.length\tA.mutation\tA.direction\tA.countA\tA.countB\tA.type\tA.gene\tA.codon\tA.syn\tA.refAa\tA.mutAa\t<AB>\tB.length\tB.mutation\tB.direction\tB.countA\tB.countB\tB.type\tB.gene\tB.codon\tB.syn\tB.refAa\tB.mutAa\tmirror count\tsupport count\ttolerance"} {print $0"\t"($17+$28)"\t"($16+$29)"\t"($4-length($15)-60)}' > ${OUTPREFIXA}_summaryBOTH.csv


echo `date` "Done!"


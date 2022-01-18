
rg="RGPL=Illumina RGPU=flowcell RGID=default_1_1 RGLB=default_1 RGSM=default"
kmerSize=11
#should be read length minus one?
padSize=100

fastaFile=data/repeats/humanRepeats.fa
bamFile=data/RRBS/51_Hf08_BmTM4_SP1_RRBS_HaeIII/51_Hf08_BmTM4_SP1_RRBS_HaeIII.MCSv3.20150825.raw.bam
outputFolder=


tmpWd=tmp2

mkdir -p $tmpWd/{db,biq}

#TODO: find out where the reads align within the file. Could they be extracted directly?

#create database
java -jar /home/karln/jars/fastaUtils.jar pad $fastaFile $padSize > $tmpWd/padded.ref.fa
gmap_build -D $tmpWd/db -d repeats $tmpWd/padded.ref.fa
cmetindex -F $tmpWd/db -d repeats

#generate kmers
java -jar /home/karln/jars/fastaUtils.jar kmerizeBis $fastaFile $kmerSize |sort -u |mawk '{printf ">seq_%s\n%s\n",NR,$0}' > $tmpWd/kmers.fa


#prefilter??
source /DEEP_fhgfs/projects/karln/software/khmer/khmerEnv/bin/activate

time java -cp .:/home/karln/jars/htsjdk-1.130.jar:/home/karln/jars/bamUtils.jar tools.bam.bamUtils toFastqRev2 $bamFile | pv -s $(samtools idxstats $bamFile |awk '{s+=$3+$4} END{print s*4}') -lcN input | python sweep-reads2.py -x 2000000 -k $kmerSize $tmpWd/kmers.fa /dev/stdin |pv -lcN filter | gsnap -n 1 -t 1 -d repeats -D $tmpWd/db --nofails --sam-use-0M -A sam --mode=cmet-stranded |pv -lcN mapped |samtools view -bS - |samtools sort - $tmpWd/out


deactivate

time java -cp .:/home/karln/jars/htsjdk-1.130.jar:/home/karln/jars/bamUtils.jar tools.bam.bamUtils bamToBiQtab $tmpWd/out.bam $tmpWd/padded.ref.fa $tmpWd/biq

mv $tmpWd/biq $outputFolder

${SCRIPTFOLDER}/wrappers/bisSNP_realignerTargetCreatorChr.sh
${SCRIPTFOLDER}/wrappers/bisSNP_realignerChr.sh
${SCRIPTFOLDER}/wrappers/bisSNP_recalibrateCountCovariates.sh
${SCRIPTFOLDER}/wrappers/bisSNP_recalibrateChr.sh
${SCRIPTFOLDER}/wrappers/bisSNP_genotypingChr.sh
${SCRIPTFOLDER}/wrappers/bisSNP_filterAndPostprocess.sh



gsnap=/TL/deep-share/archive00/software/packages/gsnap/gmap-2015-11-20/install/bin/gsnap

$gsnap --nthreads=8 --npaths=1 --merge-overlap --db=amplicons --dir=db --sam-use-0M --format=sam --mode=cmet-stranded read1_val_1.fq read2_val_2.fq |java -jar ~/jars/AddOrReplaceReadGroups.jar RGPL=Illumina RGPU=flowcell RGID=default_1_1 RGLB=default_1 RGSM=default I=/dev/stdin o=tmp.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT SORT_ORDER=coordinate


/usr/lib/jvm/java-6-openjdk-amd64/jre/bin/java -jar ~/jars/BisSNP-0.82.2.jar -R ref.fa -I mapped.sorted.bam -T BisulfiteRealignerTargetCreator -o mapped.sorted.indel_target_interval.intervals -nt 8 -S LENIENT
/usr/lib/jvm/java-6-openjdk-amd64/jre/bin/java -jar ~/jars/BisSNP-0.82.2.jar -R ref.fa -I mapped.sorted.bam -T BisulfiteIndelRealigner -targetIntervals mapped.sorted.indel_target_interval.intervals  -o mapped.sorted.realigned.bam -S LENIENT -cigar --maxReadsInMemory 1500000


/usr/lib/jvm/java-6-openjdk-amd64/jre/bin/java -jar ~/jars/BisSNP-0.82.2.jar --maximum_read_cov 100000 -R ref.fa -T BisulfiteGenotyper -I mapped.sorted.realigned.bam  -vfn1 mapped.cpg.raw.vcf -vfn2 mapped.snp.raw.vcf -S LENIENT -nt 14 -out_modes EMIT_VARIANT_AND_CYTOSINES -stand_emit_conf 0 -C CG,1 -C CA,1 -C CC,1 -C CT,1 -C CH,1 -C CAG,1 -C CHH,1 -C CHG,1 -C GC,2 -C GCH,2 -C GCG,2 -C HCG,2 -C HCH,2 -C HCA,2 -C HCC,2 -C HCT,2


perl /TL/deep-share/nobackup/deep_svn/Deep/bisulfite/trunk/bisulfitePipeline/src/scripts/third-party/bis-SNP_Utils/sortByRefAndCor.pl --k 1 --c 2  mapped.snp.raw.vcf ref.fa.fai > mapped.snp.raw.sorted.vcf
perl /TL/deep-share/nobackup/deep_svn/Deep/bisulfite/trunk/bisulfitePipeline/src/scripts/third-party/bis-SNP_Utils/sortByRefAndCor.pl --k 1 --c 2 mapped.cpg.raw.vcf ref.fa.fai > mapped.cpg.raw.sorted.vcf 


/usr/lib/jvm/java-6-openjdk-amd64/jre/bin/java -jar ~/jars/BisSNP-0.82.2.jar -maxCov 100000 -R ref.fa -T VCFpostprocess -oldVcf mapped.snp.raw.sorted.vcf -newVcf mapped.snp.filtered.vcf -snpVcf mapped.snp.raw.sorted.vcf -o mapped.snp.raw.filter.summary.txt

/usr/lib/jvm/java-6-openjdk-amd64/jre/bin/java -jar ~/jars/BisSNP-0.82.2.jar -maxCov 100000 -R ref.fa -T VCFpostprocess -oldVcf mapped.cpg.raw.sorted.vcf -newVcf mapped.cpg.filtered.vcf -snpVcf mapped.snp.raw.sorted.vcf -o mapped.cpg.raw.filter.summary.txt 


perl /TL/deep-share/nobackup/deep_svn/Deep/bisulfite/trunk/bisulfitePipeline/src/scripts/tools/vcf2bed.NOME.pl mapped.cpg.filtered.vcf GCH > mapped.cpg.filtered.GCH.bed 



self, bam_path, gpc_path, or_path, fasta, motifs, mapq, proc, folder, cutoff











mainFolder=/DEEP_fhgfs/projects/karln/hpseq/150929/out/hpRRBS_02/wd

#fastq
rawFolder=$mainFolder/00_raw
trimFolder=$mainFolder/05_trim
#sync and join mixed in one folder
joinFolder=$mainFolder/10_join

#bam
mapFolder=$mainFolder/15_map


echo readCount_raw $(cat $rawFolder/seq_R1_????.fastq.gz.readCount |mawk '{s+=$1} END{print s}')
echo readCount_trimmed $(cat $trimFolder/seq_R1_????_val_1.fq.gz.readCount |mawk '{s+=$1} END{print s}')
echo readCount_joined $(cat $joinFolder/seq_R1_????.fq.gz.readCount |mawk '{s+=$1} END{print s}')
echo readCount_unsynced $(cat $joinFolder/seq_????_R1_unsync.fq.gz.readCount |mawk '{s+=$1} END{print s}')
echo readCount_synced $(cat $joinFolder/seq_????_R1_sync.fq.gz.readCount |mawk '{s+=$1} END{print s}')



zcat $s |mawk -vOFS='\t' '

$7=="other" {
 otherMeth+=$5
 otherTot+=$5+$6
 next
}

$7=="GC" {
 gcMeth+=$8+ ($9+$10)/2
 gcTot+= $8 +$9 +$10+$11
 gcHemi+=$9+$10
 next
}

$7=="CG" {
 cgMeth+=$8+ ($9+$10)/2
 cgTot+= $8 +$9 +$10+$11
 cgHemi+=$9+$10
 cgCov+=$4
 cgN+=1
}

$7=="CG" && $5+$6>4 {
 CG5+=1
}

$7=="CG" && $5+$6>9 {
 CG10+=1
}

END {
 print "gcMeth",gcMeth/gcTot
 print "gcHemi",gcHemi/gcTot
 print "otherMeth",otherMeth/otherTot
 print "cgNcov5",CG5
 print "cgNcov10",CG10
 print "cgNtot",cgN
 print "cgAvgCov",cgCov/cgN
 print "cgMeth",cgMeth/cgTot
 print "cgHemi",cgHemi/cgTot
}
'







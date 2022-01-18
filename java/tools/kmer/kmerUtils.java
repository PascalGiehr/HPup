package tools.kmer;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.zip.GZIPInputStream;

import cern.colt.list.IntArrayList;
import cern.colt.map.OpenIntObjectHashMap;

import com.google.common.hash.BloomFilter;
import com.google.common.hash.Funnel;
import com.google.common.hash.PrimitiveSink;

import tools.aligning.LocalMutation;
import tools.aligning.LocalMutationList;
import tools.blast.blastM8Alignment;
import tools.blast.blastM8Parser;
import tools.fasta.FastaSeq;
import tools.fasta.fastaKmerParser;
import tools.fasta.fastaParser;
import tools.fastq.FastqSeq;
import tools.fastq.fastqParser;
import tools.gff.gffLine;
import tools.kmer.KmerSet_binary;
import tools.kmer.bisulfite.BisulfiteKmerSet2;
import tools.kmer.bisulfite.BisulfiteTarget;
import tools.kmer.bisulfite.BisulfiteTargetKey;
//import tools.shoreMap.Mutation;
//import tools.shoreMap.ShoreMapLine;
import tools.kmer.bisulfite.BisulfiteKmerSet;
import tools.sequences.sequenceUtils;
import utils.perfectHash.EasyPerfectMinimalHashing;

public class kmerUtils {

	private static HashMap<String, String> methodsHelp= new HashMap<String, String>();
	private final static String sep="\t";
	private final static String[] nucleotides= new String[]{"A","C","G","T","N"};
	
	public static void main(String[] args)throws Exception{
		//load helpMap
		methodsHelp.put("sortKmers", "sortKmers - sorts the kmers\n\targs = <kmer file> <kmer position> <bin size> <tmpPrefix>\n");
		methodsHelp.put("splitGraph", "splitGraph - generates a list of clusters for each read pair given k\n\targs = <fastq file 1> <fastq file 2> <kmerSize> <nrOfReads>\n");
		methodsHelp.put("sortKmersLink", "sortKmersLink - sorts the kmers and adds a flag if it is linked to the previous read\n\targs = <kmer file> <kmer position> <max jump>\n");
		methodsHelp.put("generateKmers", "generateKmers - generates all kmers of length n\n\targs = <n>\n");
		methodsHelp.put("generateSeeds", "generateSeeds - generates seeds from a kmer file\n\targs = <kmerFile> <kmer position> <min seed size> <outPrefix>\n");
		methodsHelp.put("generateSeedsNoLong", "generateSeedsNoLong - generates seeds from a kmer file\n\targs = <kmerFile> <kmer position> <min seed size> <outPrefix>\n");
		methodsHelp.put("generateSeeds2", "generateSeeds2 - generates seeds from a kmer file and the contrast file\n\targs = <kmerFile> <contrast kmerFile> <kmer position> <min seed size> <outPrefix>\n");
		methodsHelp.put("generateSeeds3", "generateSeeds3 - generates seeds from a kmer file and the contrast file. finding starting points with l-mers\n\targs = <kmerFile> <contrast kmerFile> <kmer position> <min seed size> <lmerSize> <outPrefix>\n");
//		methodsHelp.put("mapListMutationKmers", "mapListMutationKmers - takes a map.list and a list of mutations (Chr\\tposition\\tnewNuc) and prints two files, one with the kmers of size n to remove and one with the ones to add\n\targs = <map.list> <mutation file> <kmer size n> <outPrefix>\n");
		methodsHelp.put("kmerCoverageTmp", "kmerCoverageTmp - prints the counts for each kmer in the sequence together with the position (0-based)\n\targs = <kmerCount> <fastaFile>\n");
		methodsHelp.put("kmerCoverageFromScratch", "kmerCoverageFromScratch - prints the counts for each kmer in the sequence together with the position (0-based)\n\targs = <fastaFile> <kmerSize>\n");
		methodsHelp.put("kmerize", "kmerize - prints all kmers of size n in each sequence\n\targs = <fastqFile> <n>\n");
		methodsHelp.put("kmerizeEnd", "kmerizeEnd - prints the kmers in the end of a read missing when reducing from n to m\n\targs = <fastqFile> <n> <m>\n");
		methodsHelp.put("reduce", "reduce - reduces the kmers at pos in a file to length n\n\targs = <fastqFile> <pos> <n>\n");
		methodsHelp.put("extractGood", "extractGood - Takes a kmerFile with the kmers at pos (column count starts at 0), and extract all fastq seqs that contains at least min of the good kmers\n\targs = <kmerFile> <pos> <fastqFile> <min> <outPrefix>\n");
		methodsHelp.put("checkEnds", "checkEnds - Takes two seed files. Builds a hash with kmer-length n for both end and beginning for the first file and checks which kmers from the second file that matches this\n\targs = <kmerFile1> <kmerFile2> <kmer length> \n");
		methodsHelp.put("multiKmerMergePipe", "multiKmerMergePipe - takes a pipe of multiple sorted kmer-files and writes a table\n\targs = <number of samples>\n");
		methodsHelp.put("multiKmerMergeJellyfishPipe", "multiKmerMergeJellyfishPipe - takes a pipe of multiple tagged jellyfish kmer-files and writes a table\n\targs = <number of samples>\n");
		methodsHelp.put("betaCall", "betaCall - reads a fasta files and counts the kmers overlapping the targets while collapsing with regard to mehtylation\n\targs = <fasta (CT)> <targets> <methylated> <unmethylated> <kmerSize>\n");
		methodsHelp.put("testPerfectHash", "testPerfectHash - reads the targets while collapsing with regard to mehtylation, and tries to create a perfect hash\n\targs = <targets> <methylated> <unmethylated> <kmerSize>\n");
		methodsHelp.put("blastSNPFilterPipe", "blastSNPFilterPipe - takes the output from blastn -outfmt '6 std qseq sseq' and filters for paired seeds\n\targs = <kmerSize> <tolerance>\n");
		methodsHelp.put("blastRevCompFilterPipe", "blastRevCompFilterPipe - takes the output from blastn -outfmt '6 std qlen slen' and filters for identical seeds. Prints the cleaned fasta file\n\targs = <fasta file>\n");
		methodsHelp.put("blastTransferAnnotationPipe", "blastTransferAnnotationPipe - takes the output from blastn -outfmt '6 std qlen slen' and the output from the SNP/InDel Filters (or another mutation annotation script) and Prints the transfered annotation\n\targs = <annotation file> <mutation column=12>\n");
		methodsHelp.put("blastTransferAnnotationToChrPipe", "blastTransferAnnotationToChrPipe - takes the output from blastn -outfmt '6 std' and the output from the SNP/InDel Filters (or another mutation annotation script) and Prints the transfered annotation\n\targs = <annotation file> <mutation column>\n");
		methodsHelp.put("blastMutationEMSAnnotationPipe", "blastMutationEMSAnnotationPipe - Adds a string of mutationDirections for EMS mutations\n\targs = <annotated blast file>\n");
		methodsHelp.put("blastRemoveRepetitive", "blastRemoveRepetitive(Z) - Takes two blast files from blastSNPFilterPipe and prints only those that are unique (only one hit in each) \n\targs = <annotated blast file1> <annotated blast file2> <out1> <out2>\n");
		methodsHelp.put("blastSeedExtensionPairingPipe", "blastSeedExtensionPairingPipe - takes the output from blastn -outfmt '6 std qlen slen' and a perfect match. If none exists the old sequence is printed. \n\targs = <query fasta file> <target fasta file>\n");
		methodsHelp.put("blastMutatedHits", "blastMutatedHits - takes the output from blastn -outfmt '6 std' and a file with the mutation annotation. Prints all hits that contains a mutation given the flexibility f. \n\targs = <blast file> <mutation annotation file> <flexibility>\n");
		methodsHelp.put("BSpartitionTargetReads", "BSpartitionTargetReads - Reads targets and cluster them with single linkage on common collapsed kmers. Splits into bins with around binSize targets. Then takes reads and sorts them to the correct bin\n\targs = <targetFa> <methylated> <unmethylated> <kmerSize> <binSize> <outPrefix>\n");
		methodsHelp.put("BSstrandFilter", "BSstrandFilter - reads paired fastq files and  the genome. Stores each strand of the genome in a separate bloom filter in order to assign the reads to the separate strands\n\targs = <fastq 1 (CT)> <fastq 2 (GA)> <bloomFilter> <methylated> <unmethylated> <kmerSize> <outPrefix>\n");
		methodsHelp.put("BSstrandFilter_generateBloomFilter", "BSstrandFilter_generateBloomFilter - Generates the bloom filter for the BSstrandFilter \n\targs = <genome> <methylated> <unmethylated> <kmerSize> <outPrefix>\n");
		methodsHelp.put("BSstrandFilter_prepJellyfish", "BSstrandFilter_prepJellyfish - Generates the files to shovle into jellyfish \n\targs = <fastq1> <fastq2> <methylated> <unmethylated> <kmerSize> <outPrefix>\n");
		methodsHelp.put("getGFF3lines", "getGFF3lines - prints lines overlapping the positions in column 2 (chr) and 3 (pos) of the input file. \n\targs = <position file> <GFF3 file>\n");
		methodsHelp.put("extractSeedData", "extractSeedData - extract all reads with pairs covered by a distribution. Also creates separate files for partially covered read sections, completely covered reads and when both reads are covered \n\targs = <kmerFile> <kmer pos> <fastq_1 1> ... <fastq_1 n> <fastq_2 1> ... <fastq_2 n> <outPrefix>\n");
		methodsHelp.put("digiNorm", "digiNorm - titus brown normalization \n\targs = <kmerFile> <kmer pos> <fastq_1 1> <fastq_2 1> <kmer size> <normalization level> <outPrefix>\n");
		methodsHelp.put("kmerTrim", "kmerTrim - takes a kmer distribution (count\\tkmer) from the standard in and trims low coverage kmer ends  \n\targs = <fastq_1 1> <fastq_2 1> <good kmer coverage> <min length> <outPrefix>\n");
		
		//check which method to run
		if(args.length>0){
			if(args[0].equals("sortKmers")&&args.length==5){
				sortKmers(args[1],Integer.parseInt(args[2]),Integer.parseInt(args[3]),args[4]);
			}else if(args[0].equals("sortKmersLink")&&args.length==4){
				sortKmersLink(args[1],Integer.parseInt(args[2]),Integer.parseInt(args[3]));
			}else if(args[0].equals("splitGraph")&&args.length>4){
				splitGraph(new BufferedReader(new FileReader(args[1])),new BufferedReader(new FileReader(args[2])),Integer.parseInt(args[3]),Integer.parseInt(args[4]));
			}else if(args[0].equals("betaCall")&&args.length==6){
				betaCall(args[1],args[2],args[3].charAt(0),args[4].charAt(0),Integer.parseInt(args[5]));
			}else if(args[0].equals("BSpartitionTargetReads")&&args.length==7){
				BSpartitionTargetReads(args[1], args[2].charAt(0), args[3].charAt(0), 
						Integer.parseInt(args[4]), Integer.parseInt(args[5]), args[6]);
			}else if(args[0].equals("BSstrandFilter")&&args.length==8){
				BSstrandFilter(args[1],args[2],args[3],args[4].charAt(0),args[5].charAt(0),Integer.parseInt(args[6]), args[7]);
			}else if(args[0].equals("BSstrandFilter_generateBloomFilter")&&args.length==6){
				BSstrandFilter_generateBloomFilter(args[1],args[2].charAt(0),args[3].charAt(0),Integer.parseInt(args[4]), args[5]);
			}else if(args[0].equals("BSstrandFilter_prepJellyfish")&&args.length==7){
				BSstrandFilter_prepJellyfish(args[1],args[2],args[3].charAt(0),args[4].charAt(0),Integer.parseInt(args[5]), args[6]);
			}else if(args[0].equals("testPerfectHash")&&args.length==5){
				testPerfectHash(args[1],args[2].charAt(0),args[3].charAt(0),Integer.parseInt(args[4]));
			}else if(args[0].equals("digiNorm")&&args.length==6){
				digiNorm(args[1],args[2],Integer.parseInt(args[3]),Integer.parseInt(args[4]),args[5]);
			}else if(args[0].equals("digiNorm2")&&args.length==6){
				digiNorm2(args[1],args[2],args[3],Integer.parseInt(args[4]),args[5]);
			}else if(args[0].equals("generateKmers")&&args.length==2){
				generateKmers(Integer.parseInt(args[1]),"");
			}else if(args[0].equals("generateSeeds")&&args.length==5){
				generateSeeds(args[1],Integer.parseInt(args[2]),Integer.parseInt(args[3]),args[4],true);
			}else if(args[0].equals("generateSeedsNoLong")&&args.length==5){
				generateSeeds(args[1],Integer.parseInt(args[2]),Integer.parseInt(args[3]),args[4],false);
			}else if(args[0].equals("generateSeeds2")&&args.length==6){
				generateSeeds2(args[1],args[2],Integer.parseInt(args[3]),Integer.parseInt(args[4]),args[5]);
			}else if(args[0].equals("generateSeeds3")&&args.length==7){
				generateSeeds3(args[1],args[2],Integer.parseInt(args[3]),Integer.parseInt(args[4]),Integer.parseInt(args[5]),args[6]);
//			}else if(args[0].equals("mapListMutationKmers")&&args.length==5){
//				mapListMutationKmers(args[1],args[2],Integer.parseInt(args[3]),args[4]);
			}else if(args[0].equals("reduce")&&args.length==4){
				reduce(args[1],Integer.parseInt(args[2]),Integer.parseInt(args[3]));
			}else if(args[0].equals("kmerCoverageFromScratch")&&args.length==3){
				kmerCoverageFromScratch(args[1],Integer.parseInt(args[2]));
			}else if(args[0].equals("kmerCoverageFromScratch2")&&args.length==3){
				kmerCoverageFromScratch2(args[1],Integer.parseInt(args[2]));
			}else if(args[0].equals("kmerCoverageFromScratch3")&&args.length==3){
				kmerCoverageFromScratch3(args[1],Integer.parseInt(args[2]));
			}else if(args[0].equals("kmerCoverageFromScratch4")&&args.length==3){
				kmerCoverageFromScratch4(args[1],args[2]);
			}else if(args[0].equals("kmerCoverageTmp")&&args.length==3){
				kmerCoverageTmp(args[1],args[2]);
			}else if(args[0].equals("kmerize")&&args.length==3){
				kmerize(args[1],Integer.parseInt(args[2]));
			}else if(args[0].equals("kmerizeEnd")&&args.length==4){
				kmerizeEnd(args[1],Integer.parseInt(args[2]),Integer.parseInt(args[3]));
			}else if(args[0].equals("kmerTrim")&&args.length==5){
				kmerTrim(args[1],args[2],Integer.parseInt(args[3]),Integer.parseInt(args[4]),args[5]);
			}else if(args[0].equals("extractGood")&&args.length==6){
				extractGood(readKmers_revComp(args[1],Integer.parseInt(args[2])),args[3],Integer.parseInt(args[4]),args[5]);
			}else if(args[0].equals("checkEnds")&&args.length==4){
				checkEnds(args[1],args[2],Integer.parseInt(args[3]));
			}else if(args[0].equals("multiKmerMergePipe")&&args.length>1){
				multiKmerMergePipe(Integer.parseInt(args[1]));
			}else if(args[0].equals("multiKmerMergeJellyfishPipe")&&args.length>1){
				multiKmerMergeJellyfishPipe(Integer.parseInt(args[1]));
			}else if(args[0].equals("multiKmerMergeJellyfishPipe2")&&args.length>1){
				multiKmerMergeJellyfishPipe2(Integer.parseInt(args[1]));
			}else if(args[0].equals("blastMutatedHits")&&args.length>3){
				blastMutationHits(new BufferedReader(new FileReader(args[1])),args[2],Integer.parseInt(args[3]));
			}else if(args[0].equals("blastSNPFilterPipe")&&args.length>2){
				blastSNPFilterPipe(new BufferedReader(new InputStreamReader(System.in)),Integer.parseInt(args[1]),Integer.parseInt(args[2]));
			}else if(args[0].equals("blastRemoveRepetitive")&&args.length>4){
				blastRemoveRepetitive(new BufferedReader(new FileReader(args[1])),new BufferedReader(new FileReader(args[2])),args[3],args[4]);
			}else if(args[0].equals("blastRemoveRepetitiveZ")&&args.length>4){
				BufferedReader in1=new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(args[1]))));
				BufferedReader in2=new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(args[2]))));
				blastRemoveRepetitive(in1,in2,args[3],args[4]);
			}else if(args[0].equals("blastRevCompFilterPipe")&&args.length>1){
				blastRevCompFilterPipe(new BufferedReader(new InputStreamReader(System.in)),args[1]);
			}else if(args[0].equals("blastSeedExtensionPairingPipe")&&args.length>1){
				blastSeedExtensionPairingPipe(new BufferedReader(new InputStreamReader(System.in)),args[1],args[2]);
			}else if(args[0].equals("blastTransferAnnotationPipe")&&args.length>1){
				blastTransferAnnotationPipe(new BufferedReader(new InputStreamReader(System.in)),args[1],12);
			}else if(args[0].equals("blastTransferAnnotationPipe")&&args.length>2){
				blastTransferAnnotationPipe(new BufferedReader(new InputStreamReader(System.in)),args[1],Integer.parseInt(args[2]));
			}else if(args[0].equals("blastTransferAnnotationToChrPipe")&&args.length>2){
				blastTransferAnnotationToChrPipe(new BufferedReader(new InputStreamReader(System.in)),args[1],Integer.parseInt(args[2]));
			}else if(args[0].equals("blastMutationEMSAnnotationPipe")&&args.length>0){
				blastMutationEMSAnnotationPipe(new BufferedReader(new InputStreamReader(System.in)));
			}else if(args[0].equals("getGFF3lines")&&args.length>2){
				getGFF3lines(args[1],args[2]);
			}else if(args[0].equals("extractSeedData")&&args.length>5){
				ArrayList<String> files1= new ArrayList<String>(), files2= new ArrayList<String>();
				int n=(args.length-3)/2;
				for(int i=3;i<n+3;++i){
					files1.add(args[i]);
				}
				for(int i=n+3;i<2*n+3;++i){
					files2.add(args[i]);
				}
				extractSeedData(readKmers(args[1],Integer.parseInt(args[2])), files1, files2, args[args.length-1]);
			}else if(args[0].equals("method2")&&args.length>1){
				
			}else{
				System.err.println(printHelp(args[0]));
			}
		}else{
			System.err.println(printHelp());
		}
	}
	
	public static void kmerTrim(String fastq1, String fastq2, int goodCov,int minLength,String outPrefix) throws Exception{
		
		
	}
	
	
	public static void BSpartitionTargetReads(String targetFa, final char methylated, final char unmethylated, 
			final int kmerSize, final int binSize, String outPrefix)throws Exception{
		
		fastaParser fp= new fastaParser(targetFa);
		FastaSeq fs;
		HashMap<BitSet, Integer> map= new HashMap<BitSet, Integer>(100000000);
		IntArrayList clusters= new IntArrayList();
		KmerSet_binary_utils kt_utils= new KmerSet_binary_utils(kmerSize);
		HashSet<BitSet> kmers=new HashSet<BitSet>();
		HashSet<Integer> putativeClusters= new HashSet<Integer>();
		
		int count=0;
		//cluster targets
		System.err.println("Loading targets ...");
		for(;fp.hasNext();++count){
			if(count%10000==0){
				System.err.print("processed targets: "+count+"\r");
			}
			fs=fp.next();
			if(fs.length()>=kmerSize){
				kmers.clear();
				kt_utils.kmerizeNoNForward(fs.getSeq().replace(methylated, unmethylated), kmers);
				Integer repCluster=clusters.size();
				clusters.add(repCluster);
				putativeClusters.clear();
				putativeClusters.add(repCluster);
				//find cluster
				for(BitSet kmer : kmers){
					Integer curClust=map.get(kmer);
					if(curClust!=null){
						//identify parent cluster
						int nextCluster=clusters.get(curClust);
						while(nextCluster>-1){
							curClust=nextCluster;
							nextCluster=clusters.get(curClust);
						}
						
						if(curClust<repCluster){
							repCluster=curClust;
						}
						putativeClusters.add(curClust);
					}
				}
				//update pointers
				putativeClusters.remove(repCluster);
				clusters.set(repCluster, -1);
				for(Integer cluster : putativeClusters){
					clusters.set(cluster, repCluster);
				}
				for(BitSet kmer : kmers){
					map.put(kmer, repCluster);
				}
			}
		}
		System.err.print("processed targets: "+count+"\n");
		fp.close();
		//Print original cluster association
		fp= new fastaParser(targetFa);
		BufferedWriter statOut= new BufferedWriter(new FileWriter(outPrefix+".stat.csv"));
		System.err.println("Writing clusters ...");
		for(count=0;fp.hasNext();++count){
			if(count%10000==0){
				System.err.print(count+"\r");
			}
			Integer clusterId=count;
			int nextCluster=clusters.get(clusterId);
			while(nextCluster>-1){
				clusterId=nextCluster;
				nextCluster=clusters.get(clusterId);
			}
			statOut.write(fp.next().getQname()+sep+clusterId+'\n');
		}
		System.err.print(count+"\n");
		statOut.close();
		// merge clusters into bins
		OpenIntObjectHashMap seqToClusters= new OpenIntObjectHashMap();
		System.err.println("merging clusters into bins");
		for(int i=0;i<clusters.size();++i){
			if(i%10000==0){
				System.err.print(i+"\r");
			}
			int clusterId=i;
			int nextCluster=clusters.get(clusterId);
			while(nextCluster>-1){
				clusterId=nextCluster;
				nextCluster=clusters.get(clusterId);
			}
			
			if(seqToClusters.containsKey(clusterId)){
				((IntArrayList) seqToClusters.get(clusterId)).add(i);
			}else{
				IntArrayList cluster= new IntArrayList();
				cluster.add(i);
				seqToClusters.put(clusterId, cluster);
			}
		}
		System.err.println(clusters.size());
		IntArrayList clusterIds=seqToClusters.keys();
		int size=binSize+1,repCluster=-1;
		HashMap<Integer, BufferedWriter> targetWriters= new HashMap<Integer, BufferedWriter>();
		HashMap<Integer, BufferedWriter> readWriters= new HashMap<Integer, BufferedWriter>();
		System.err.println("binning ...");
		for(int i=0;i<clusterIds.size();++i){
			IntArrayList cluster= (IntArrayList) seqToClusters.get(clusterIds.get(i));
			if(size>binSize){
				//initialize a new cluster
				repCluster=cluster.get(0);
				cluster.remove(0);
				clusters.set(repCluster, -1);
				size=1;
				targetWriters.put(repCluster, new BufferedWriter(new FileWriter(outPrefix+"."+repCluster+".targets.fa")));
				readWriters.put(repCluster, new BufferedWriter(new FileWriter(outPrefix+"."+repCluster+".reads.fa")));
			}
			size+=cluster.size();
			for(int j=0;j<cluster.size();++j){
				clusters.set(cluster.get(j), repCluster);
			}
		}
		
		//split targets
		fp= new fastaParser(targetFa);
		System.err.println("Writing binned target files ...");
		for(count=0;fp.hasNext();++count){
			if(count%10000==0){
				System.err.print(count+"\r");
			}
			repCluster=clusters.get(count);
			if(repCluster<0){
				repCluster=count;
			}
			targetWriters.get(repCluster).write(fp.next().toString()+'\n');
		}
		for(BufferedWriter out : targetWriters.values()){
			out.close();
		}
		System.err.print(count+"\n");
		fp= new fastaParser(new BufferedReader(new InputStreamReader(System.in)));
		System.err.println("processing reads");
		for(count=0;fp.hasNext();++count){
			if(count%10000==0){
				System.err.print(count+"\r");
			}
			fs=fp.next();
			if(fs.length()>=kmerSize){
				kmers.clear();
				kt_utils.kmerizeNoNForward(fs.getSeq().replace(methylated, unmethylated), kmers);
				putativeClusters.clear();
				//Check kmers and print to clusters
				for(BitSet kmer : kmers){
					Integer curClust= map.get(kmer);
					if(curClust!=null){
						repCluster=clusters.get(curClust);
						if(repCluster<0){
							repCluster=curClust;
						}
						putativeClusters.add(repCluster);
					}
				}
				for(Integer clusterId : putativeClusters){
					readWriters.get(clusterId).write(fs.toString()+'\n');
				}
			}
		}
		System.err.print(count+"\n");
		for(BufferedWriter out : readWriters.values()){
			out.close();
		}
		System.err.println("Done!");
	}
	

	
	public static void testPerfectHash(String targetsFa, final char methylated, final char unmethylated, final int kmerSize)throws Exception{
		fastaParser fp= new fastaParser(new BufferedReader(new FileReader(targetsFa)));
		FastaSeq fs;
		BitSet kmerBin;
		StringBuffer kmerOrig;
		HashSet<BitSet> keys= new HashSet<BitSet>(70000000);
		KmerSet_binary_utils kt_utils= new KmerSet_binary_utils(kmerSize);
		
		int count=0;
		for(;fp.hasNext();++count){
			if(count%1000==0){
				System.err.print("Read #targets: "+count+"\r");
			}
			fs=fp.next();
			final String seq=fs.getSeq();
			if(seq.length()>=kmerSize){
				kmerOrig= new StringBuffer(seq.substring(0, kmerSize));
				kmerBin= kt_utils.stringToBitSet(kmerOrig.toString().replace(methylated, unmethylated));
				keys.add(kmerBin);
				for(int i=kmerSize;i<seq.length();++i){
					kmerOrig.deleteCharAt(0);
					final char curNuc=seq.charAt(i);
					kmerOrig.append(curNuc);
					if(curNuc==methylated){
						kmerBin=kt_utils.shiftBitSet(kmerBin, unmethylated);
					}else{
						kmerBin=kt_utils.shiftBitSet(kmerBin, curNuc);
					}
					keys.add(kmerBin);
				}
			}
			

		}
		System.err.println("Read #targets: "+count);
		System.err.println("There were "+keys.size()+" keys");
		

		int mb = 1024*1024;
		
		System.gc();
		Runtime runtime= Runtime.getRuntime();
		System.err.println("##### Heap utilization statistics [MB] #####");
        
        //Print used memory
        System.err.println("before BMZ: "
            + (runtime.totalMemory() - runtime.freeMemory()) / mb);
 
        EasyPerfectMinimalHashing<BitSet> epmh= (EasyPerfectMinimalHashing<BitSet>) EasyPerfectMinimalHashing.createEasyPerfectMinimalHashing(keys); 
//		BMZ<BitSet> bmz= (BMZ<BitSet>) BMZ.create(keys);
		System.gc();
		System.err.println("After BMZ: "
	            + (runtime.totalMemory() - runtime.freeMemory()) / mb);

		keys.clear();
		System.gc();
		System.err.println("keys dropped: "
	            + (runtime.totalMemory() - runtime.freeMemory()) / mb);
		
		
//		System.out.println("seed1: "+bmz.getSeed1()+'\t'+"seed2: "+bmz.getSeed2());
	}
	
	public static void BSstrandFilter_prepJellyfish(String fastq1, String fastq2, final char methylated, final char unmethylated, final int kmerSize, String outPrefix) throws Exception{
		fastqParser fqp1= new fastqParser(fastq1),
				fqp2= new fastqParser(fastq2);
		int count=0;
		String def="N";
		while(def.length()<kmerSize){
			def+="N";
		}
		
		BufferedWriter outName= new BufferedWriter(new FileWriter(outPrefix+".names")),
				outKmer= new BufferedWriter(new FileWriter(outPrefix+".kmer"));
		
		for(;fqp1.hasNext()&&fqp2.hasNext();++count){
			if(count%10000==0){
				System.err.print("read pairs: "+count+"\r");
			}
			FastqSeq fqs1=fqp1.next(),
					fqs2=fqp1.next().reverseComplement();
			String q1= fqs1.getQname(),
					q2= fqs2.getQname(),
					s1= fqs1.getSeq().toUpperCase(),
					s2= fqs2.getSeq().toUpperCase();
			if(s1.length()>kmerSize){
				StringBuffer k1= new StringBuffer(s1.substring(0, kmerSize).replace(methylated, unmethylated));
				outName.write(q1+"\n");
				outKmer.write(k1.toString()+"\n");
				for(int i=kmerSize;i<s1.length();++i){
					char nuc= s1.charAt(i);
					k1.deleteCharAt(0);
					if(nuc==methylated){
						k1.append(unmethylated);
					}else{
						k1.append(nuc);
					}
					outName.write(q1+"\n");
					outKmer.write(k1.toString()+"\n");
				}
			}
			
			if(s2.length()>kmerSize){
				StringBuffer k2= new StringBuffer(s2.substring(0, kmerSize).replace(methylated, unmethylated));
				outName.write(q2+"\n");
				outKmer.write(k2.toString()+"\n");
				for(int i=kmerSize;i<s2.length();++i){
					char nuc= s2.charAt(i);
					k2.deleteCharAt(0);
					if(nuc==methylated){
						k2.append(unmethylated);
					}else{
						k2.append(nuc);
					}
					outName.write(q1+"\n");
					outKmer.write(k2.toString()+"\n");
				}
			}
		}
		System.err.println("read pairs: "+count);
		
		outKmer.close();
		outName.close();
		fqp1.close();
		fqp2.close();
	}
	
	public static void BSstrandFilter_generateBloomFilter(String genomeFa, final char methylated, final char unmethylated, 
			final int kmerSize, String outFile) throws Exception{
		final int mb= 1024 * 1024;
		Runtime runtime= Runtime.getRuntime();
		System.err.println("Start: "
	            + (runtime.totalMemory() - runtime.freeMemory()) / mb + " Mb");
		Funnel<String> funnel= new Funnel<String>() {

			private static final long serialVersionUID = 1L;

			@Override
			public void funnel(String arg0, PrimitiveSink into) {
				// TODO Auto-generated method stub
				into.putString(arg0);
			}
			
		};
		
		
		BloomFilter<String> forward= BloomFilter.create(funnel,2147483647 , 0.01),
				reverse= BloomFilter.create(funnel,2147483647 , 0.01);
		
		System.err.println("Bloom filters created: "
	            + (runtime.totalMemory() - runtime.freeMemory()) / mb + " Mb");
		
		fastaKmerParser fkp= new fastaKmerParser(new BufferedReader(new FileReader(genomeFa)), kmerSize);
		long count=0L;
		final String regex="[^ATGCatgc]";
		for(;fkp.hasNext();){
			if(count%100000==0){
				System.err.print("Genomic k-mers read: "+count+"\r");
			}
			String kmer=fkp.next();
			if(!kmer.matches(regex)){ //throw away everything that contains non-atgc chars
				++count;
				kmer=kmer.toUpperCase();
				forward.put(kmer.replace(methylated, unmethylated));
				reverse.put(sequenceUtils.reverseComplement(kmer).replace(methylated, unmethylated));
			}
		}
		fkp.close();
		System.err.println("Genomic k-mers read: "+count);
		System.err.println("Bloom filters filled: "
	            + (runtime.totalMemory() - runtime.freeMemory()) / mb + " Mb");
		
		ObjectOutputStream oos= new ObjectOutputStream(new FileOutputStream(outFile));
		oos.writeObject(forward);
		oos.writeObject(reverse);
		oos.close();
		System.err.println("Done!");
		
	}
	
	public static void BSstrandFilter(String fastq1, String fastq2, String bloomFilterFile, 
			final char methylated, final char unmethylated, 
			final int kmerSize, final String outPrefix)throws Exception{
		final int mb = 1024*1024;
		Runtime runtime= Runtime.getRuntime();
		System.err.println("Start: "
	            + (runtime.totalMemory() - runtime.freeMemory()) / mb + " Mb");
//		Funnel<String> funnel= new Funnel<String>() {
//
//			private static final long serialVersionUID = 1L;
//
//			@Override
//			public void funnel(String arg0, PrimitiveSink into) {
//				// TODO Auto-generated method stub
//				into.putString(arg0);
//			}
//			
//		};
		
		ObjectInputStream ois= new ObjectInputStream(new FileInputStream(bloomFilterFile));
		
		BloomFilter<String> forward= (BloomFilter<String>) ois.readObject();
		BloomFilter<String> reverse= (BloomFilter<String>) ois.readObject();
		
		ois.close();
		
//		BloomFilter<String> forward= BloomFilter.create(funnel, Integer.MAX_VALUE, 0.01),
//				reverse= BloomFilter.create(funnel, Integer.MAX_VALUE, 0.01);
		
		System.err.println("Bloom filters created: "
	            + (runtime.totalMemory() - runtime.freeMemory()) / mb + " Mb");
		
		System.err.println("Expected false positive probability: \n\tforward: "+ 
				forward.expectedFpp() + "\n\treverse: "+reverse.expectedFpp());
		
//		fastaKmerParser fkp= new fastaKmerParser(new BufferedReader(new FileReader(genomeFa)), kmerSize);
		long count=0L;
//		final String regex="[^ATGCatgc]";
//		for(;fkp.hasNext();){
//			if(count%100000==0){
//				System.err.print("Genomic k-mers read: "+count+"\r");
//			}
//			String kmer=fkp.next();
//			if(!kmer.matches(regex)){ //throw away everything that contains non-atgc chars
//				++count;
//				kmer=kmer.toUpperCase();
//				forward.put(kmer.replace(methylated, unmethylated));
//				reverse.put(sequenceUtils.reverseComplement(kmer).replace(methylated, unmethylated));
//			}
//		}
//		fkp.close();
//		System.err.println("Genomic k-mers read: "+count);
//		System.err.println("Bloom filters filled: "
//	            + (runtime.totalMemory() - runtime.freeMemory()) / mb + " Mb");
		
		fastqParser fqp1= new fastqParser(fastq1),
				fqp2= new fastqParser(fastq2);
		FastqSeq fqs1,fqs2;
		
		BufferedWriter outF= new BufferedWriter(new FileWriter(outPrefix+".F.fq")),
				outU= new BufferedWriter(new FileWriter(outPrefix+".U.fq")),
				outR= new BufferedWriter(new FileWriter(outPrefix+".R.fq")),
				outStat= new BufferedWriter(new FileWriter(outPrefix+".stat.fq"));
		
		count=0L;
		for(;fqp1.hasNext()&&fqp2.hasNext();++count){
			fqs1=fqp1.next();
			fqs2=fqp2.next().reverseComplement();
			if (fqs1.length()>=kmerSize && fqs2.length()>=kmerSize) {
				int ff = 0, fr = 0, rf = 0, rr = 0;
				//count matching kmers for first read
				String seq = fqs1.getSeq().replace(methylated, unmethylated);
				StringBuffer kmerOrig = new StringBuffer(seq.substring(0,
						kmerSize));
				if (forward.mightContain(kmerOrig.toString())) {
					++ff;
				}
				
				if (reverse.mightContain(kmerOrig.toString())) {
					++fr;
				}
				for (int i = kmerSize; i < seq.length(); ++i) {
					kmerOrig.deleteCharAt(0);
					final char curNuc = seq.charAt(i);
					kmerOrig.append(curNuc);
					if (forward.mightContain(kmerOrig.toString())) {
						++ff;
					}
					if (reverse.mightContain(kmerOrig.toString())) {
						++fr;
					}
				}
				//count matching kmers for second read
				seq = fqs2.getSeq().replace(methylated, unmethylated);
				kmerOrig = new StringBuffer(seq.substring(0, kmerSize));
				if (forward.mightContain(kmerOrig.toString())) {
					++rf;
				}
				if (reverse.mightContain(kmerOrig.toString())) {
					++rr;
				}
				for (int i = kmerSize; i < seq.length(); ++i) {
					kmerOrig.deleteCharAt(0);
					final char curNuc = seq.charAt(i);
					kmerOrig.append(curNuc);
					if (forward.mightContain(kmerOrig.toString())) {
						++rf;
					}
					if (reverse.mightContain(kmerOrig.toString())) {
						++rr;
					}
				}
				//Decide
				
				//TODO: The number of matching kmers could be used to refine the decisions
				outStat.write(fqs1.getQname()+sep+ff+sep+fr+sep+rf+sep+rr);
				
				
				if (ff > fr) {
					//First read forward strand
					if (rf > rr) {
						//second read forward strand
						//pair forward
						outF.write(fqs1.toString() + "\n");
						outF.write(fqs2.toString() + "\n");
						outStat.write(sep+"ff_rf_F\n");
					} else if (rf == rr) {
						//second read undecided
						//pair forward
						outF.write(fqs1.toString() + "\n");
						outF.write(fqs2.toString() + "\n");
						outStat.write(sep+"ff_ru_F\n");
					} else {
						//second read reverse strand
						//pair unknown
						outU.write(fqs1.toString() + "\n");
						outU.write(fqs2.toString() + "\n");
						outStat.write(sep+"ff_rr_U\n");
					}
				} else if (ff == fr) {
					//First read undecided
					if (rf > rr) {
						//second read forward strand
						//pair forward
						outF.write(fqs1.toString() + "\n");
						outF.write(fqs2.toString() + "\n");
						outStat.write(sep+"fu_rf_F\n");
					} else if (rf == rr) {
						//second read undecided
						//pair unknown
						outU.write(fqs1.toString() + "\n");
						outU.write(fqs2.toString() + "\n");
						outStat.write(sep+"fu_ru_U\n");
					} else {
						//second read reverse strand
						//pair reverse
						outR.write(fqs1.toString() + "\n");
						outR.write(fqs2.toString() + "\n");
						outStat.write(sep+"fu_rr_R\n");
					}
				} else {
					//First read reverse strand
					if (rf > rr) {
						//second read forward strand
						//pair unknown
						outU.write(fqs1.toString() + "\n");
						outU.write(fqs2.toString() + "\n");
						outStat.write(sep+"fr_rf_U\n");
					} else if (rf == rr) {
						//second read undecided
						//pair reverse
						outR.write(fqs1.toString() + "\n");
						outR.write(fqs2.toString() + "\n");
						outStat.write(sep+"fr_ru_R\n");
					} else {
						//second read reverse strand
						//pair reverse
						outR.write(fqs1.toString() + "\n");
						outR.write(fqs2.toString() + "\n");
						outStat.write(sep+"fr_rr_R\n");
					}
				}
			}
			
			
		}
		
		outF.close();
		outU.close();
		outR.close();
		outStat.close();
		
		System.err.println("Start: "
		            + (runtime.totalMemory() - runtime.freeMemory()) / mb + " Mb");
	}
	
	public static void betaCall(String fastaFile,String targetsFa, final char methylated, final char unmethylated, final int kmerSize)throws Exception{
		System.err.println("Get keys ...");
		HashSet<BitSet> 		keys= new HashSet<BitSet>(700000);
		final int 				mb = 1024*1024;
		
		BisulfiteKmerSet2 		kmers;
		fastaParser 			fp_targets= new fastaParser(targetsFa);
		BisulfiteTarget 		singleTarget;
		ArrayList<Integer> 		positions;
		BitSet 					kmerBin;
		StringBuffer 			kmerOrig;
		
		KmerSet_binary_utils 	kt_utils= new KmerSet_binary_utils(kmerSize);
		
		FastaSeq 				fs;
		BisulfiteTarget[] 		targets;
//		ArrayList<BisulfiteTarget> targets= new ArrayList<BisulfiteTarget>();
		fastaParser 			fp_reads= new fastaParser(fastaFile);
		
		Runtime 				runtime= Runtime.getRuntime();
		//Print used memory
        System.err.println("Start: "
            + (runtime.totalMemory() - runtime.freeMemory()) / mb + " Mb");
		
		int count=0;
		//goes through the targets and collects all collapsed target kmers
		for(;fp_targets.hasNext();++count){
			if(count%1000==0){
				System.err.print("Read #targets: "+count+"\r");
			}
			fs=fp_targets.next();
			final String seq=fs.getSeq();
			if(seq.length()>=kmerSize){
//				++targetNr;
				kmerOrig= new StringBuffer(seq.substring(0, kmerSize));
				kmerBin= kt_utils.stringToBitSet(kmerOrig.toString().replace(unmethylated, methylated));
				keys.add(kmerBin);
				for(int i=kmerSize;i<seq.length();++i){
					kmerOrig.deleteCharAt(0);
					final char curNuc=seq.charAt(i);
					kmerOrig.append(curNuc);
					if(curNuc==unmethylated){
						kmerBin=kt_utils.shiftBitSet(kmerBin, methylated);
					}else{
						kmerBin=kt_utils.shiftBitSet(kmerBin, curNuc);
					}
					keys.add(kmerBin);
				}
			}
		}
		fp_targets.close();
		targets= new BisulfiteTarget[count];
		
		System.err.println("Read #targets: "+count);
		System.err.println("There were "+keys.size()+" keys");
		
		System.err.println("got keys: "
	            + (runtime.totalMemory() - runtime.freeMemory()) / mb +" Mb");
		kmers= new BisulfiteKmerSet2(kmerSize,EasyPerfectMinimalHashing.createEasyPerfectMinimalHashing(keys),unmethylated,methylated,count); 
		keys.clear();
		
		System.err.println("created kmer structure: "
	            + (runtime.totalMemory() - runtime.freeMemory()) / mb +" Mb");
		
		count=0;
		fp_targets= new fastaParser(new BufferedReader(new FileReader(targetsFa)));
//		long lastmem=runtime.totalMemory() - runtime.freeMemory(),curmem;
		
		//reads the targets again and loads the targets
		for(int targetNr=0;fp_targets.hasNext();++count,++targetNr){
			
			fs=fp_targets.next();
			positions= new ArrayList<Integer>();
			for(int i=0;i<fs.length()-1;++i){
				//TODO: modify to allow other patterns
				if(fs.getSeq().charAt(i)=='C' && fs.getSeq().charAt(i+1)=='G'){
					positions.add(i);
				}
			}
			
			singleTarget=new BisulfiteTarget(fs.getQname(), fs.getSeq(), kmerSize, positions);
//			singleTarget=targets[count];
			final String seq=fs.getSeq();
			if(seq.length()>=kmerSize){
				kmerOrig= new StringBuffer(seq.substring(0, kmerSize));
				kmerBin= kt_utils.stringToBitSet(kmerOrig.toString().replace(unmethylated, methylated));
				singleTarget.add(kmers.addKmerTarget(kmerOrig, kmerBin, new BisulfiteTargetKey(targetNr,0)));
				
				for(int i=kmerSize;i<seq.length();++i){
					kmerOrig.deleteCharAt(0);
					final char curNuc=seq.charAt(i);
					kmerOrig.append(curNuc);
					if(curNuc==unmethylated){
						kmerBin=kt_utils.shiftBitSet(kmerBin, methylated);
					}else{
						kmerBin=kt_utils.shiftBitSet(kmerBin, curNuc);
					}
					singleTarget.add(kmers.addKmerTarget(kmerOrig, kmerBin, new BisulfiteTargetKey(targetNr, i-kmerSize+1)));
				}
			}
			singleTarget.trimToSize(kmers);
			targets[targetNr]=singleTarget;
			if(count%1000==0){
				System.err.print("Read #targets: "+count+"\r");
			}
//			curmem=runtime.totalMemory() - runtime.freeMemory();
//			System.out.println(targetNr+"\t"+(curmem-lastmem));
//			lastmem=curmem;
		}
		fp_targets.close();
		
		
		System.err.println("Read #targets: "+count);
		System.gc();
		System.err.println("targets loaded: "
	            + (runtime.totalMemory() - runtime.freeMemory()) / mb +" Mb");
		
		System.err.println("Trimming and filtering targets ...");
		for(BisulfiteTarget target : targets){
			target.trimToSize(kmers);
		}
		System.gc();
		System.err.println("targets trimmed: "
	            + (runtime.totalMemory() - runtime.freeMemory()) / mb +" Mb");
//		System.exit(616);
		
		kmers.trimToSize();
		
//		kmers.transform();
		System.gc();
		System.err.println("kmers transformed: "
	            + (runtime.totalMemory() - runtime.freeMemory()) / mb +" Mb");
		
	
		
		//generate equivalences
		//add equivalences to targets
		//generate cliques
		//how to store distinct kmers?
		// - read them from jellyfish and do it when building the matrices?
		// - add another layer to BisulfiteCollapsedKmerSet
		
		System.err.println("Processing reads ...");
		count=0;
		int added=0;
		for(;fp_reads.hasNext();++count){
			if(count%10000==0){
				System.err.print(count+"\r");
			}
			fs=fp_reads.next();
			final String seq=fs.getSeq();
			
			if(seq.length()>=kmerSize){
				kmerOrig= new StringBuffer(seq.substring(0, kmerSize));
				kmerBin= kt_utils.stringToBitSet(kmerOrig.toString().replace(unmethylated, methylated));
				kmers.addAvailableKmer(kmerOrig, kmerBin);
				for(int i=kmerSize;i<seq.length();++i){
					kmerOrig.deleteCharAt(0);
					final char curNuc=seq.charAt(i);
					kmerOrig.append(curNuc);
					if(curNuc==unmethylated){
						kmerBin=kt_utils.shiftBitSet(kmerBin, methylated);
					}else{
						kmerBin=kt_utils.shiftBitSet(kmerBin, curNuc);
					}
					added+=kmers.addAvailableKmer(kmerOrig, kmerBin)?1:0;
				}
			}
//			for(String orig : kmerFunctions.kmers(fqs1.getSeq(), kmerSize)){
//				kmers.addAvailableKmer(orig);
//			}
		}
		System.err.println("reads loaded: "
	            + (runtime.totalMemory() - runtime.freeMemory()) / mb +" Mb");
		System.err.println(count+"");
		
		System.err.println("#unique kmers added: "+added);
		
		System.err.println("generating equvivalences");
		kmers.generateEquvivalences();
		
		for(int i=0;i<kmers.equivalences.size();++i){
			IntArrayList group=(IntArrayList) kmers.equivalences.get(i);
			if(group.size()==1){
				//simply add upp the counts
				singleTarget= targets[group.get(0)];
				singleTarget.printToStdOut("+", kmers);
			}else{
				//need to solve the distribution problem
				//TODO: implement
			}
		}
		
		
		System.exit(0);
		
//		System.exit(616);
		System.err.println("Writing output ...");
		count=0;
		for(BisulfiteTarget bt : targets){
			++count;
			if(count%1000==0){
				System.err.print(count+"\r");
			}
			if(bt==null){
				System.err.println(count);
			}
			bt.printToStdOut("F",kmers);
		}
		System.err.println(count+"");
		
		System.err.println("Done!");
		
		
	}
	
	public static void digiNorm2(String fastq1,String fastq2,String kmerFile,int targetCov,String outPrefix) throws Exception{
		HashMap<BitSet, Integer> counts= new HashMap<BitSet, Integer>(1000000000);
		HashSet<BitSet> rep= new HashSet<BitSet>(1000000000);
		BufferedReader kmerIn= new BufferedReader(new FileReader(kmerFile));
		String s=kmerIn.readLine();
		KmerSet_binary_utils ku= new KmerSet_binary_utils(s.length());
		counts.put(ku.kmerToUse(s), targetCov);
		fastqParser fqp1= new fastqParser(fastq1);
		fastqParser fqp2= new fastqParser(fastq2);
		FastqSeq fq1,fq2;
		ArrayList<BitSet> curKmers= new ArrayList<BitSet>();
		BufferedWriter out1=new BufferedWriter(new FileWriter(outPrefix+"red"+targetCov+"_1.fq")),
				out2=new BufferedWriter(new FileWriter(outPrefix+"red"+targetCov+"_2.fq"));
		int repCount,kmerCount,count=1,printed=0;
		
		System.err.println("load rep kmers");
		for(s=kmerIn.readLine();s!=null;s=kmerIn.readLine(),++count){
			if(count%10000==0)
				System.err.print(count+"\r");
			counts.put(ku.kmerToUse(s), targetCov);
		}
		kmerIn.close();
		System.err.println(count);
		
		count=0;
		System.err.println("filter reads");
		for(;fqp1.hasNext() && fqp2.hasNext();++count){
			if(count%10000==0){
				System.err.print(count+"   "+printed+"\r");
			}
			fq1=fqp1.next();
			fq2=fqp2.next();
			curKmers.clear();
			ku.kmerizeNoN(fq1.getSeq(), curKmers);
			ku.kmerizeNoN(fq2.getSeq(), curKmers);
			repCount=0;
			for(BitSet kmer : curKmers){
				if(rep.contains(kmer)){
					++repCount;
				}else if(counts.containsKey(kmer)){
					kmerCount=counts.get(kmer);
					if(kmerCount==0){
						counts.remove(kmer);
						rep.add(kmer);
					}else{
						counts.put(kmer, kmerCount-1);
					}
				}
			}
			if(repCount<curKmers.size()/2){
				++printed;
				out1.write(fq1.toString()+"\n");
				out2.write(fq2.toString()+"\n");
			}
		}
		System.err.println(count+"   "+printed);
		out1.close();
		out2.close();
	}
	
	public static void digiNorm(String fastq1,String fastq2,int kmerSize,int targetCov,String outPrefix)throws Exception{
		HashSet<BitSet> rep= new HashSet<BitSet>(2000000000);
		HashMap<BitSet, Integer> counts= new HashMap<BitSet, Integer>(2000000000);
		fastqParser fqp1= new fastqParser(fastq1);
		fastqParser fqp2= new fastqParser(fastq2);
		FastqSeq fq1,fq2;
		KmerSet_binary_utils ku=new KmerSet_binary_utils(kmerSize);
		ArrayList<BitSet> curKmers= new ArrayList<BitSet>();
		int repCount,kmerCount,count=0,printed=0;
		BufferedWriter out1=new BufferedWriter(new FileWriter(outPrefix+"red"+targetCov+"_1.fq")),
				out2=new BufferedWriter(new FileWriter(outPrefix+"red"+targetCov+"_2.fq"));
		System.err.println("start to normalize ...");
		
		for(;fqp1.hasNext() && fqp2.hasNext();++count){
			if(count%10000==0){
				System.err.print(count+"   "+printed+"\r");
			}
			fq1=fqp1.next();
			fq2=fqp2.next();
			curKmers.clear();
			ku.kmerizeNoN(fq1.getSeq(), curKmers);
			ku.kmerizeNoN(fq2.getSeq(), curKmers);
			repCount=0;
			for(BitSet kmer : curKmers){
				if(rep.contains(kmer)){
					++repCount;
				}else if(counts.containsKey(kmer)){
					kmerCount=counts.get(kmer);
					if(kmerCount==targetCov){
						counts.remove(kmer);
						rep.add(kmer);
					}else{
						counts.put(kmer, kmerCount+1);
					}
				}else{
					counts.put(kmer, 1);
				}
			}
			if(repCount<curKmers.size()/2){
				++printed;
				out1.write(fq1.toString()+"\n");
				out2.write(fq2.toString()+"\n");
			}
		}
		System.err.println(count+"   "+printed);
		out1.close();
		out2.close();
	}
	
	public static void kmerCoverageFromScratch4(String faFile,String kmerCountFile) throws Exception{
		HashSet<String> rep= new HashSet<String>(2000000000);
		final String one="1";
		BufferedReader in= new BufferedReader(new FileReader(kmerCountFile));
		String[] l=in.readLine().split("\t");
		if(!one.equals(l[0])){
			rep.add(l[1]);
		}
		fastaKmerParser fkp= new fastaKmerParser(new BufferedReader(new FileReader(faFile)), l[1].length()); 
		int i=0,curPos=0;
		String curname="",kmer;
		System.err.println("extracting repetitive kmers");
		for(String s=in.readLine();s!=null;s=in.readLine(),++i){
			if(i%100000==0){
				System.err.print("     "+i+"  "+rep.size()+"\r");
			}
			l=s.split("\t");
			if(!one.equals(l[0])){
				rep.add(l[1]);
			}
		}
		System.err.println("     "+i+"  "+rep.size());
		i=1;
		System.err.println("annotating file ...");
		for(;fkp.hasNext();++i){
			if(i%100000==0){
				System.err.print("     "+i+"\r");
			}
			kmer=kmerFunctions.kmerToUse(fkp.next());
			if(!curname.equals(fkp.getQname())){
				curPos=0;
				curname=fkp.getQname();
			}
			System.out.println(curname+sep+(curPos++)+sep+(rep.contains(kmer)?2:1));
		}
		System.err.println("     "+i);
		
	}
	
	public static void kmerCoverageFromScratch3(String faFile, int kmerSize)throws Exception{
		HashSet<String> uniqueA= new HashSet<String>(1000000000),uniqueC= new HashSet<String>(1000000000),
				uniqueG= new HashSet<String>(1000000000),
				uniqueT= new HashSet<String>(1000000000),
				uniqueN= new HashSet<String>(1000000000),
				rep= new HashSet<String>(2000000000);
//		BloomFilter<String> unique= new BloomFilter<String>(1000000000, 2000000000);
//		BloomFilter<String> rep= new BloomFilter<String>(1000000000, 1000000000);
//		System.err.println("Rep FP: "+rep.getFalsePositiveProbability()+" K: "+rep.getK());
//		System.err.println("unique FP: "+unique.getFalsePositiveProbability()+" K: "+unique.getK());
		fastaKmerParser fkp= new fastaKmerParser(new BufferedReader(new FileReader(faFile)), kmerSize);
		String kmer;
		boolean novel;
		
		System.err.println("counting kmers ...");
		int i=0;
		for(;fkp.hasNext();++i){
			if(i%100000==0){
				System.err.print("     "+i+"  "+uniqueA.size()+"  "+rep.size()+"\r");
			}
			kmer=kmerFunctions.kmerToUse(fkp.next());
			switch (kmer.charAt(0)) {
			case 'A':
				novel=uniqueA.add(kmer);
				break;
			case 'C':
				novel=uniqueC.add(kmer);
				break;
			case 'G':
				novel=uniqueG.add(kmer);
				break;
			case 'T':
				novel=uniqueT.add(kmer);
				break;

			default:
				novel=uniqueN.add(kmer);
				break;
			}
			if(!novel){
				rep.add(kmer);
			}
		}
		System.err.println("     "+i+"  "+uniqueA.size()+"  "+rep.size());
		
		uniqueA.clear();
		uniqueC.clear();
		uniqueG.clear();
		uniqueT.clear();
		uniqueN.clear();
		
		fkp= new fastaKmerParser(new BufferedReader(new FileReader(faFile)), kmerSize);
		kmer= kmerFunctions.kmerToUse(fkp.next());
		String curname=fkp.getQname();
		int curPos=0;
		System.out.println(curname+sep+(curPos++)+sep+(rep.contains(kmer)?2:1));
		i=1;
		System.err.println("annotating file ...");
		for(;fkp.hasNext();++i){
			if(i%100000==0){
				System.err.print("     "+i+"\r");
			}
			kmer=kmerFunctions.kmerToUse(fkp.next());
			if(!curname.equals(fkp.getQname())){
				curPos=0;
				curname=fkp.getQname();
			}
			System.out.println(curname+sep+(curPos++)+sep+(rep.contains(kmer)?2:1));
		}
		System.err.println("     "+i);
	}
	
	public static void kmerCoverageFromScratch2(String faFile, int kmerSize)throws Exception{
		HashSet<String> unique= new HashSet<String>(2000000000), rep= new HashSet<String>(2000000000);
//		BloomFilter<String> unique= new BloomFilter<String>(1000000000, 2000000000);
//		BloomFilter<String> rep= new BloomFilter<String>(1000000000, 1000000000);
//		System.err.println("Rep FP: "+rep.getFalsePositiveProbability()+" K: "+rep.getK());
//		System.err.println("unique FP: "+unique.getFalsePositiveProbability()+" K: "+unique.getK());
		fastaKmerParser fkp= new fastaKmerParser(new BufferedReader(new FileReader(faFile)), kmerSize);
		String kmer;
		
		System.err.println("counting kmers ...");
		int i=0;
		for(;fkp.hasNext();++i){
			if(i%100000==0){
				System.err.print("     "+i+"  "+unique.size()+"  "+rep.size()+"\r");
			}
			kmer=kmerFunctions.kmerToUse(fkp.next());
			if(unique.add(kmer)){
				//new kmer
			}else{
				rep.add(kmer);
			}
		}
		System.err.println("     "+i+"  "+unique.size()+"  "+rep.size());
		
		unique.clear();
		fkp= new fastaKmerParser(new BufferedReader(new FileReader(faFile)), kmerSize);
		kmer= kmerFunctions.kmerToUse(fkp.next());
		String curname=fkp.getQname();
		int curPos=0;
		System.out.println(curname+sep+(curPos++)+sep+(rep.contains(kmer)?2:1));
		i=1;
		System.err.println("annotating file ...");
		for(;fkp.hasNext();++i){
			if(i%100000==0){
				System.err.print("     "+i+"\r");
			}
			kmer=kmerFunctions.kmerToUse(fkp.next());
			if(!curname.equals(fkp.getQname())){
				curPos=0;
				curname=fkp.getQname();
			}
			System.out.println(curname+sep+(curPos++)+sep+(rep.contains(kmer)?2:1));
		}
		System.err.println("     "+i);
	}
	
	public static void kmerCoverageFromScratch(String faFile, int kmerSize) throws Exception{
		HashMap<String, Integer> map= new HashMap<String, Integer>(2000000000);
		fastaKmerParser fkp= new fastaKmerParser(new BufferedReader(new FileReader(faFile)), kmerSize);
		String kmer;
		
		System.err.println("counting kmers ...");
		int i=0;
		for(;fkp.hasNext();++i){
			if(i%100000==0){
				System.err.print("     "+i+"\r");
			}
			kmer=kmerFunctions.kmerToUse(fkp.next());
			if(map.containsKey(kmer)){
				map.put(kmer, map.get(kmer)+1);
			}else{
				map.put(kmer, 1);
			}
		}
		System.err.println("     "+i);
		
		
		fkp= new fastaKmerParser(new BufferedReader(new FileReader(faFile)), kmerSize);
		kmer= kmerFunctions.kmerToUse(fkp.next());
		String curname=fkp.getQname();
		int curPos=0;
		System.out.println(curname+sep+(curPos++)+sep+map.get(kmer));
		i=1;
		System.err.println("annotating file ...");
		for(;fkp.hasNext();++i){
			if(i%100000==0){
				System.err.print("     "+i+"\r");
			}
			kmer=kmerFunctions.kmerToUse(fkp.next());
			if(!curname.equals(fkp.getQname())){
				curPos=0;
				curname=fkp.getQname();
			}
			System.out.println(curname+sep+(curPos++)+sep+map.get(kmer));
		}
		System.err.println("     "+i);
	}
	
	public static void kmerCoverageTmp(String kmerCountFile,String faFile) throws Exception{
		HashMap<String, String> map= new HashMap<String, String>(1000000000);
		BufferedReader in=new BufferedReader(new FileReader(kmerCountFile));
		String[] l= in.readLine().split("\t");
		fastaKmerParser fkp= new fastaKmerParser(new BufferedReader(new FileReader(faFile)), l[0].length());
		String kmer= kmerFunctions.kmerToUse(fkp.next());
		String curname=fkp.getQname();
		int curPos=0;
		map.put(l[1], l[0]);
		
		
		
		for (String s=in.readLine();s!=null;s=in.readLine()){
			l=s.split("\t");
			map.put(l[1], l[0]);
		}
		
		System.out.println(curname+sep+(curPos++)+sep+map.get(kmer));
		
		for(;fkp.hasNext();){
			kmer=kmerFunctions.kmerToUse(fkp.next());
			if(!curname.equals(fkp.getQname())){
				curPos=0;
				curname=fkp.getQname();
			}
			System.out.println(curname+sep+(curPos++)+sep+map.get(kmer));
		}
	}
	
	private static int[] readPairs;
//	private static ArrayList<Integer> readPairs;
	
	public static void splitGraph(BufferedReader file1, BufferedReader file2, int kmerSize, int nrOfReadPairs)throws Exception{
		HashMap<BitSet, Integer> map= new HashMap<BitSet, Integer>(1000000000);
		HashSet<Integer> clustersToMerge= new HashSet<Integer>(),clusters= new HashSet<Integer>();
		readPairs= new int[nrOfReadPairs];
//		readPairs=new ArrayList<Integer>(2*nrOfReadPairs);
		fastqParser fqp1= new fastqParser(file1, ""), fqp2= new fastqParser(file2,"");
		FastqSeq fs1,fs2;
		int sequenceCount=0,firstCluster;
		HashSet<BitSet> curKmers= new HashSet<BitSet>();
		KmerSet_binary_utils ku= new KmerSet_binary_utils(kmerSize-1);
		System.err.println("finding union ...");
		boolean verbose=false;
		for(;fqp1.hasNext()&&fqp2.hasNext()&&sequenceCount<nrOfReadPairs;++sequenceCount){
			if(sequenceCount%10000==0){
				System.err.print("     "+sequenceCount+"\r");
			}
//			verbose=sequenceCount>4798;
			fs1=fqp1.next();
			fs2=fqp2.next();
			//kmerize reads
			if(verbose) System.err.println("kmerize");
			curKmers.clear();
			splitGraph_kmerizeNoN(fs1.getSeq(), ku, curKmers);
			splitGraph_kmerizeNoN(fs2.getSeq(), ku, curKmers);
			if(verbose) System.err.println(curKmers.size()+" kmers");
			//find clusters to merge
			clustersToMerge.clear();
			firstCluster=sequenceCount;
			for (BitSet bs : curKmers) {
				if(map.containsKey(bs)){
					clustersToMerge.add(splitGraph_findRoot( map.get(bs),verbose));
				}
			}
			if(verbose) System.err.println(clustersToMerge.size()+" clusters to merge");
			readPairs[firstCluster]=-1;
//			readPairs.add(-1);
			if(clustersToMerge.size()!=0){
				for(Integer cluster2 : clustersToMerge){
					firstCluster= splitGraph_unite( firstCluster, cluster2,verbose);
				}
			}
			if(verbose) System.err.println("updating kmers in map");
			//set the links for all (k-1)-mers
			for (BitSet bs : curKmers){
				map.put(bs, firstCluster);
			}
			if(verbose) System.err.println("map contains "+ map.size()+ " kmers");
		}
		
		System.err.println("     "+sequenceCount);
		map.clear();
		//print
		for(int i=0;i<sequenceCount;++i){
			firstCluster=splitGraph_findRoot( i);
			clusters.add(firstCluster);
			System.out.println(firstCluster);
		}
		System.err.println("Generated: "+clusters.size()+" clusters");
	}
	
	private static int splitGraph_unite( int cluster1, int cluster2,boolean verbose){
		cluster1=splitGraph_findRoot( cluster1);
		cluster2=splitGraph_findRoot( cluster2);
		if(cluster1 == cluster2){
			return cluster1;
		}
//		int c1=readPairs.get(cluster1),c2=readPairs.get(cluster2);
//		if(c1<c2){
//			readPairs.set(cluster1, c1+c2);
//			readPairs.set(cluster2, cluster1);
//			return cluster1;
//		}else{
//			readPairs.set(cluster2, c1+c2);
//			readPairs.set(cluster1, cluster2);
//			return cluster2;
//		}
		
		if(readPairs[cluster1]< readPairs[cluster2]){
			readPairs[cluster1]+=readPairs[cluster2];
			readPairs[cluster2]=cluster1;
			return cluster1;
		}else{
			readPairs[cluster2]+=readPairs[cluster1];
			readPairs[cluster1]=cluster2;
			return cluster2;
		}
	}
	
	private static int splitGraph_findRoot(int start){
		return splitGraph_findRoot(start,false);
	}
	
	private static int splitGraph_findRoot( int start, boolean verbose){
		int i=start,levels=0;
		if(verbose) System.err.println("finding root for: "+ start);
//		while(readPairs.get(i)>0){
//			if(verbose) System.err.println(levels+"\t"+i);
//			++levels;
//			i=readPairs.get(i);
//		}
		
		while (readPairs[i]>=0){
//			readPairs[i]= readPairs[readPairs[i]];
			++levels;
			i= readPairs[i];
		}
		if(levels<2){
			return i;
		}
		if (verbose) System.err.println("simplifying path");
		int end=i;
		i=start;
		int old;
//		while(readPairs.get(i)>0){
//			old=i;
//			i=readPairs.get(i);
//			readPairs.set(old, end);
//		}
		
		while(readPairs[i]>=0){
			old=i;
			i=readPairs[i];
			readPairs[old]=end;
		}
		return end;
	}
	
	private static void splitGraph_kmerizeNoN(String seq, KmerSet_binary_utils ku, HashSet<BitSet> out){
		String[] l= seq.split("N+");
		BitSet bs;
		for(String s: l){
			if(ku.getKmerSize()<=s.length()){
				bs=ku.stringToBitSet(s.substring(0, ku.getKmerSize()));
				out.add(ku.kmerToUse(bs));
				for(int i=ku.getKmerSize();i<s.length();++i){
					bs= ku.shiftBitSet(bs, s.charAt(i));
					out.add(ku.kmerToUse(bs));
				}
			}
		}
	}
	
	public static void extractSeedData(KmerSet_binary kmers, ArrayList<String> files1,ArrayList<String> files2, String outPrefix) throws Exception{
		if(files1.size()>0 && files2.size()>0){
			if(files1.size()==files2.size()){
				System.err.println("gz file postfix -> gzipped file");
				BufferedReader fq1br,fq2br;
				BufferedWriter allPaired1= new BufferedWriter(new FileWriter(outPrefix+"_all_1.fq")),
						allPaired2= new BufferedWriter(new FileWriter(outPrefix+"_all_2.fq")),
						singlet1= new BufferedWriter(new FileWriter(outPrefix+"_singlet_1.fq")),
						singlet2= new BufferedWriter(new FileWriter(outPrefix+"_singlet_2.fq")),
						both1= new BufferedWriter(new FileWriter(outPrefix+"_perfect_1.fq")),
						both2= new BufferedWriter(new FileWriter(outPrefix+"_perfect_2.fq")),
						subseqs= new BufferedWriter(new FileWriter(outPrefix+"_subSeqs.fa"));
				for(int i=0;i<files1.size();++i){
					System.err.println(files1.get(i)+"\t"+files2.get(i));
					if(files1.get(i).endsWith("gz")){
						System.err.println("1: gzipped");
						fq1br= new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(files1.get(i)))));
					}else{
						fq1br= new BufferedReader(new FileReader(files1.get(i)));
					}
					if(files2.get(i).endsWith("gz")){
						System.err.println("2: gzipped");
						fq2br= new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(files2.get(i)))));
					}else{
						fq2br= new BufferedReader(new FileReader(files2.get(i)));
					}
					extractSeedData(kmers,fq1br,fq2br,allPaired1,allPaired2,singlet1,singlet2,both1,both2,subseqs);
				}
				allPaired1.close();
				allPaired2.close();
				singlet1.close();
				singlet2.close();
				both1.close();
				both2.close();
				subseqs.close();
			}else{
				System.err.println("differenly large file lists");
			}
		}else{
			System.err.println("No fastq files");
		}
	}
	
	public static void extractSeedData(KmerSet_binary kmers, BufferedReader fq1br, BufferedReader fq2br,
			BufferedWriter allPaired1, BufferedWriter allPaired2, BufferedWriter singlet1,
			BufferedWriter singlet2, BufferedWriter both1, BufferedWriter both2, BufferedWriter subseqs)throws Exception{
		boolean[] coverage;
		boolean perfect1,perfect2,covered1,covered2;
		int seqCounter=0;
		StringBuffer subseq,seq;
		fastqParser fq1= new fastqParser(fq1br, ""), fq2= new fastqParser(fq2br, "");
		FastqSeq fqs1,fqs2;
		
		for(;fq1.hasNext() && fq2.hasNext();){
			fqs1=fq1.next();
			fqs2=fq2.next();
			coverage=kmers.goodKmers(fqs1.getSeq());
			perfect1=true;
			covered1=false;
			subseq=new StringBuffer();
			seq=new StringBuffer(fqs1.getSeq());
			for(int i=0,j=kmers.getKmerSize();i<coverage.length;++i,++j){
				if(coverage[i]){
					if(!covered1){
						covered1=true;
					}
					if(subseq.length()==0){
						//initiate new subseq
						subseq.append(seq.subSequence(i, j));
					}else{
						//extend current sequence
						subseq.append(seq.charAt(j-1));
					}
				}else{
					if(perfect1) perfect1=false;
					if(subseq.length()>0){
						subseqs.write(">"+(++seqCounter)+"\n"+subseq+"\n");
						subseq= new StringBuffer();
					}
				}
			}
			if(subseq.length()>0){
				if(subseq.length()<seq.length()){
					subseqs.write(">"+(++seqCounter)+"\n"+subseq+"\n");
				}
				subseq= new StringBuffer();
			}
			perfect2=true;
			covered2=false;
			coverage=kmers.goodKmers(fqs2.getSeq());
			seq= new StringBuffer(fqs2.getSeq());
			for(int i=0,j=kmers.getKmerSize();i<coverage.length;++i,++j){
				if(coverage[i]){
					if(!covered2){
						covered2=true;
					}
					if(subseq.length()==0){
						//initiate new subseq
						subseq.append(seq.subSequence(i, j));
					}else{
						//extend current sequence
						subseq.append(seq.charAt(j-1));
					}
				}else{
					if(perfect2) perfect2=false;
					if(subseq.length()>0){
						subseqs.write(">"+(++seqCounter)+"\n"+subseq+"\n");
						subseq= new StringBuffer();
					}
				}
			}
			if(subseq.length()>0){
				if(subseq.length()<seq.length()){
					subseqs.write(">"+(++seqCounter)+"\n"+subseq+"\n");
				}
				subseq= new StringBuffer();
			}
			if(covered1 || covered2){
				allPaired1.write(fqs1.toString()+"\n");
				allPaired2.write(fqs2.toString()+"\n");
				if(perfect1){
					if(perfect2){
						both1.write(fqs1.toString()+"\n");
						both2.write(fqs2.toString()+"\n");
					}else{
						singlet1.write(fqs1.toString()+"\n");
					}
				}else{
					if(perfect2){
						singlet2.write(fqs2.toString()+"\n");
					}
				}
			}
		}
	}
	
	public static void getGFF3lines(String positionFile,String gffFile)throws Exception{
		HashMap<String, HashMap<Integer, ArrayList<String>>> positions= new HashMap<String, HashMap<Integer,ArrayList<String>>>();
		HashMap<Integer, ArrayList<String>> posTmp;
		ArrayList<String> names;
		gffLine gl;
		
		BufferedReader in= new BufferedReader(new FileReader(positionFile));
		String[] l;
		int pos;
		
		for(String s=in.readLine();s!=null;s=in.readLine()){
			l=s.split("\t");
			if(positions.containsKey(l[1])){
				posTmp= positions.get(l[1]);
			}else{
				posTmp= new HashMap<Integer, ArrayList<String>>();
				positions.put(l[1], posTmp);
			}
			pos=Integer.parseInt(l[2]);
			if(posTmp.containsKey(pos)){
				names=posTmp.get(pos);
			}else{
				names= new ArrayList<String>();
				posTmp.put(pos, names);
			}
			names.add(l[0]);
		}
		in.close();
		in= new BufferedReader(new FileReader(gffFile));
		String s=in.readLine();
		while(s.startsWith("#")){
			s=in.readLine();
		}
		gl=new gffLine(s);
		if(gl.getChr().startsWith("Chr") && !positions.containsKey(gl.getChr())){
			//shift the keys of the positions hash
			//will not work for chloroplast and such
			names= new ArrayList<String>(positions.keySet());
			for(String chr : names){
				posTmp=positions.remove(chr);
				positions.put("Chr"+chr, posTmp);
			}
		}
		
		
		
		for(;s!=null;s=in.readLine()){
			if(!s.startsWith("#")){
				gl=new gffLine(s);
				if(positions.containsKey(gl.getChr())){
					posTmp=positions.get(gl.getChr());
					//ugly coding ahead
					for (Integer position : posTmp.keySet()) {
						if(gl.containsPos(position)){
							for(String name : posTmp.get(position)){
								System.out.println(name+sep+gl.toString());
							}
						}
					}
				}
			}
		}
	}
	
	public static void blastMutationEMSAnnotationPipe(BufferedReader in)throws Exception{
		String[] l;
		LocalMutationList lml;
		
		
		for(String s=in.readLine();s!=null;s=in.readLine()){
			l=s.split("\t");
			lml= new LocalMutationList(l[1]);
			System.out.println(s+"\t"+lml.mutationListEMS());
		}
	}
	
	public static void blastMutationHits(BufferedReader blast, String annotationFile, int flexibility) throws Exception{
		BufferedReader in= new BufferedReader(new FileReader(annotationFile));
		HashMap<String, LocalMutationList> mutations= new HashMap<String, LocalMutationList>();
		blastM8Parser bp= new blastM8Parser(blast);
		blastM8Alignment ba;
		LocalMutationList lml;
		String[] l;
		for(String s=in.readLine();s!=null;s=in.readLine()){
			l=s.split("\t");
			mutations.put(l[0], new LocalMutationList(l[1]));
		}
		for(;bp.hasMore();){
			ba=bp.nextHit();
			if(mutations.containsKey(ba.qname)){
				lml=mutations.get(ba.qname);
				if(ba.qstart<ba.qend){
					if(lml.encapsulated(ba.qstart-flexibility, ba.qend+flexibility)){
						System.out.println(ba.getOrigString());
					}
				}else{
					if(lml.encapsulated(ba.qend-flexibility, ba.qstart+flexibility)){
						System.out.println(ba.getOrigString());
					}
				}
			}
		}
	}
	
	public static void blastSeedExtensionPairingPipe(BufferedReader in, String queryFastaFile, String targetFastaFile) throws Exception{
		HashMap<String, FastaSeq> querySeqs= new HashMap<String, FastaSeq>();
		fastaParser fp= new fastaParser(new BufferedReader(new FileReader(queryFastaFile)));
		FastaSeq fs;
		blastM8Alignment ba=new blastM8Alignment(),bestHit=new blastM8Alignment();
		HashMap<String,String> bestHits=new HashMap<String,String>();
		HashSet<String> printed= new HashSet<String>();
		int length=0,start1,start2,end1,end2,edgeCount,qlen,slen;
		String qname="";
		String[] l;
		boolean found=false;
		for (;fp.hasNext();){
			fs=fp.next();
			querySeqs.put(fs.getQname(), fs);
		}
		for(String s=in.readLine();s!=null;s=in.readLine()){
			ba=new blastM8Alignment(s);
			if(!qname.equals(ba.qname)){
				if(found){
					//save for later
					bestHits.put(bestHit.tname,bestHit.qname);
				}else if(length>0){
					System.out.println(querySeqs.get(ba.qname));
				}
				found=false;
				length=querySeqs.get(ba.qname).length();
				qname=ba.qname;
			}
			
			l=s.split("\t");
			qlen=Integer.parseInt(l[12]);
			slen=Integer.parseInt(l[13]);
			
			if(ba.identity==100 && slen>qlen && ba.alignment_length>qlen/2){
				
				if(ba.qstart<ba.qend){
					start1=ba.qstart;
					end1=ba.qend;
				}else{
					start1=ba.qend;
					end1=ba.qstart;
				}
				if(ba.tstart<ba.tend){
					start2=ba.tstart;
					end2=ba.tend;
				}else{
					start2=ba.tend;
					end2=ba.tstart;
				}
				edgeCount=0;
				if(end2==slen){
					++edgeCount;
				}
				if(start2==1){
					++edgeCount;
				}
				if(end1<qlen){
					--edgeCount;
				}
				if(start1>1){
					--edgeCount;
				}
				
				
				if (edgeCount>=0) {
					if (!found) {
						bestHit = new blastM8Alignment(ba);
						found = true;
					} else if (ba.e_value < bestHit.e_value) {
						//shouldn't happen???
						bestHit = new blastM8Alignment(ba);
					}
				}
			}
		}
		if(found){
			//besthit found in the last row
			bestHits.put(bestHit.tname,bestHit.qname);
		}
		
		
		fp= new fastaParser(new BufferedReader(new FileReader(targetFastaFile)));
		for(;fp.hasNext();){
			fs=fp.next();
			if(bestHits.containsKey(fs.getQname())){
				System.out.println(">"+bestHits.get(fs.getQname())+" "+fs.getQname()+"\n"+fs.getSeq());
				printed.add(bestHits.get(fs.getQname()));
			}
		}
		//print the original seed for the seeds without extension
		for (FastaSeq out : querySeqs.values()) {
			if(!printed.contains(out.getQname())){
				System.out.println(out);
			}
		}
	}
	
	public static void blastRemoveRepetitive(BufferedReader in1 ,BufferedReader in2,String outFile1,String outFile2)throws Exception{
		HashMap<String, String> set1= new HashMap<String, String>();
		HashSet<String> rep1= new HashSet<String>();
		BufferedWriter out1= new BufferedWriter(new FileWriter(outFile1)), out2= new BufferedWriter(new FileWriter(outFile2));
		blastM8Alignment ba;
		String repline,line,qname="",tname;
		boolean repetitive=false;
		
		for(String s=in1.readLine();s!=null;s=in1.readLine()){
			ba= new blastM8Alignment(s);
			if(!rep1.contains(ba.qname)){
				if(set1.containsKey(ba.qname)){
					set1.remove(ba.qname);
					rep1.add(ba.qname);
				}else{
					set1.put(ba.qname, s);
				}
			}
		}
		in1.close();
		line=in2.readLine();
		repline=line;
		qname=line.split("\t")[0];
		
		for(line=in2.readLine();line!=null;line=in2.readLine()){
			ba=new blastM8Alignment(line);
			if(qname.equals(ba.qname)){
				repetitive=true;
			}else{
				if(!repetitive){
					//check other direction
					tname=repline.split("\t")[1];
					if(!rep1.contains(tname)){
						if(set1.containsKey(tname)){
							out2.write(repline+"\n");
							out1.write(set1.get(tname)+"\n");
						}
					}
				}
				//restart
				repetitive=false;
				qname=ba.qname;
				repline=line;
			}
		}
		
		if(!repetitive){
			//check other direction
			tname=repline.split("\t")[1];
			if(!rep1.contains(tname)){
				out2.write(repline+"\n");
				if(set1.containsKey(tname)){
					out1.write(set1.get(tname)+"\n");
				}
			}
		}
		
		out1.close();
		out2.close();
	}
	
	public static void blastTransferAnnotationPipe(BufferedReader in, String ablastFile, final int mutPos)throws Exception{
		blastM8Alignment ba;
		BufferedReader blast= new BufferedReader(new FileReader(ablastFile));
		HashMap<String, LocalMutationList> mutations= new HashMap<String, LocalMutationList>();
		LocalMutationList lml;
		int start1,start2,end1,end2,edgeCount,qlen,slen;
		String[] l;
		//Gather mutations
		for(String s=blast.readLine();s!=null;s=blast.readLine()){
			l=s.split("\t");
			mutations.put(l[0], new LocalMutationList(l[mutPos]));
		}
		blast.close();
		//Parse new file
		for(String s=in.readLine();s!=null;s=in.readLine()){
			ba= new blastM8Alignment(s);
			if(ba.qname.equals(ba.tname) && mutations.containsKey(ba.tname) && ba.identity==100){
				//count edges
				l=s.split("\t");
				qlen=Integer.parseInt(l[12]);
				slen=Integer.parseInt(l[13]);
				if(ba.qstart<ba.qend){
					start1=ba.qstart;
					end1=ba.qend;
				}else{
					start1=ba.qend;
					end1=ba.qstart;
				}
				if(ba.tstart<ba.tend){
					start2=ba.tstart;
					end2=ba.tend;
				}else{
					start2=ba.tend;
					end2=ba.tstart;
				}
				edgeCount=0;
				if(end2==slen){
					++edgeCount;
				}
				if(start2==1){
					++edgeCount;
				}
				if(end1<qlen){
					--edgeCount;
				}
				if(start1>1){
					--edgeCount;
				}

				
				if (edgeCount>=0) {
					//shift mutations
					lml = mutations.get(ba.tname);
					if (ba.strand[0] == '-') {
						lml.shiftRevComp(ba.tend + ba.qstart);
					} else {
						lml.shift(ba.tstart - ba.qstart);
					}
					System.out.println(ba.tname + "\t" + lml);
				}
				
			}
		}
	}
	
	
	public static void blastTransferAnnotationToChrPipe(BufferedReader in, String ablastFile, final int mutPos)throws Exception{
		blastM8Alignment ba;
		BufferedReader blast= new BufferedReader(new FileReader(ablastFile));
		HashMap<String, LocalMutationList> mutations= new HashMap<String, LocalMutationList>();
		LocalMutationList lml;
		int direction,qpos,tpos,nextMut,curMutI;
//		int start1,start2,end1,end2,edgeCount,qlen,slen;
		String[] l;
		String chrName;
		char[] qseq,tseq;
		//Gather mutations
		for(String s=blast.readLine();s!=null;s=blast.readLine()){
			l=s.split("\t");
			mutations.put(l[0], new LocalMutationList(l[mutPos]));
		}
		blast.close();
		//Parse new file
		for(String s=in.readLine();s!=null;s=in.readLine()){
			ba= new blastM8Alignment(s);
			if(mutations.containsKey(ba.qname)){
				//shift mutations
				lml = new LocalMutationList(mutations.get(ba.qname));
				lml.sort();
				nextMut=lml.get(0).getPosition();
				curMutI=0;
				qpos=ba.qstart;
				if(ba.strand[0]=='+'){
					direction=1;
					tpos=ba.tstart;
				}else{
					direction=-1;
					tpos=ba.tend;
				}
				l=ba.getOrigString().split("\t");
				qseq=l[12].toCharArray();
				tseq=l[13].toCharArray();
				for(int i=0;i<qseq.length;++i){
					if(qseq[i]!='-'){
						++qpos;
					}
					if(tseq[i]!='-'){
						tpos+=direction;
					}
					if(qpos==nextMut){
						//shift mutation
						if(direction==1){
							lml.get(curMutI).updatePos(tpos);
						}else{
							lml.get(curMutI).updatePosRevComp(tpos);
						}
						//go to next mutation
						++curMutI;
						if(curMutI<lml.size()){
							nextMut=lml.get(curMutI).getPosition();
						}else{
							break;
						}
					}
				}
				if(curMutI<lml.size()){
					System.err.println(ba.qname+" is missing "+(lml.size()-curMutI)+" mutations");
				}
				
				if (curMutI>0) {
					if(ba.tname.startsWith("Chr")){
						chrName= ba.tname.substring(3);
					}else{
						chrName= ba.tname;
					}
					for (LocalMutation lm : lml) {
						System.out.println(ba.qname + "\t" + chrName + "\t"
								+ lm.getPosition() + "\t" + lm.getRef() + "\t"
								+ lm.getMutation() +"\t1\t1\t1\t1\t1");
					}
				}
				
//				}
				
			}
		}
	}

	public static void blastRevCompFilterPipe(BufferedReader in, String faFile)throws Exception{
		blastM8Alignment ba;
		String[] l;
		String qname;
		boolean linkPrinted;
		int qlen,slen;
		HashSet<String> printed= new HashSet<String>(),tmpSet;
		HashMap<String, HashSet<String>> links= new HashMap<String, HashSet<String>>();
		for(String s=in.readLine();s!=null;s=in.readLine()){
			ba= new blastM8Alignment(s);
			if(ba.identity==100 && !ba.qname.equals(ba.tname)){
				l=s.split("\t");
				qlen=Integer.parseInt(l[12]);
				slen=Integer.parseInt(l[13]);
				//get all full length perfect hits
				if(qlen==slen && ba.alignment_length==qlen){
					if(!links.containsKey(ba.qname)){
						links.put(ba.qname, new HashSet<String>());
					}
					links.get(ba.qname).add(ba.tname);
					if(!links.containsKey(ba.tname)){
						links.put(ba.tname, new HashSet<String>());
					}
					links.get(ba.tname).add(ba.qname);
				}
			}
		}
		int added=1,rounds=1;
		while(added>0){
			added=0;
			for (String key : links.keySet()) {
				tmpSet=links.get(key);
				for(String o1 : tmpSet){
					for(String o2 : tmpSet){
						if(links.get(o1).add(o2)){
							++added;
						}
					}
				}
			}
			System.err.println("Round "+rounds+": "+added+" links added");
			++rounds;
		}
		
		fastaParser fp= new fastaParser(new BufferedReader(new FileReader(faFile)));
		FastaSeq fs;
		for(;fp.hasNext();){
			fs=fp.next();
			qname=fs.getQname();
			if(links.containsKey(qname)){
				linkPrinted=false;
				for(String name : links.get(qname)){
					if(printed.contains(name)){
						linkPrinted=true;
						break;
					}
				}
				
				if(!linkPrinted){
					if(printed.add(qname)){
						System.out.println(fs.toString());
					}
				}else{
					printed.add(qname);
				}
			}else if(printed.add(qname)){
				System.out.println(fs.toString());
			}
		}
	}
	
	public static void blastSNPFilterPipe(BufferedReader in,final int kmerSize, final int tolerance)throws Exception{
		blastM8Alignment ba;
		String[] l;
		char[] qseq,sseq;
		LocalMutationList ml= new LocalMutationList();
		final int lengthCutoff= 2*(kmerSize-1)-tolerance-1;
//		System.out.println("#lengthCutoff: "+lengthCutoff);
//		int pos;
		for(String s=in.readLine();s!=null;s=in.readLine()){
			ba= new blastM8Alignment(s);
//			System.out.println("#" +s);
//			System.out.println("#filter-values: "+(ba.alignment_length-ba.mismatches)+"\t"+ba.mismatches+"\t"+ba.gap_openings);
			if(ba.alignment_length-ba.mismatches>lengthCutoff && (ba.mismatches>0 || ba.gap_openings>0)){
				//alignment contains mismatches and no gaps. The alignment is long enough to harbor these changes
				l=s.split("\t");
				qseq=l[12].toCharArray();
				sseq=l[13].toCharArray();
				ml.clear();
				for(int i=0;i<ba.alignment_length;++i){
					if(qseq[i]!=sseq[i]){
						ml.add(new LocalMutation(i+1,sseq[i],qseq[i]));
					}
				}
//				System.out.println("#nr of mismatches: "+mismatches.size());
				if(ml.get(0).getPosition()+(ba.alignment_length-ml.get(ml.size()-1).getPosition())>lengthCutoff){
					System.out.println(ba+"\t"+ml);
				}
			}
		}
	}
	
	public static void multiKmerMergePipe(final int size)throws Exception{
		BufferedReader in = new BufferedReader(new InputStreamReader(System.in));
		ArrayList<String> samples= new ArrayList<String>(size);
		HashMap<String, String> counts= new HashMap<String, String>(size);
		final String zero=sep+"0";
		String s=in.readLine();
		String[] l=s.split("\t");
		String kmer=l[1];
		samples.add(l[2]);
		counts.put(l[2], l[0]);
		//pre-load
		for(s=in.readLine();s!=null;s=in.readLine()){
			l=s.split("\t");
			if(kmer.equals(l[1])){
				//store
				if(!samples.contains(l[2])){
					samples.add(l[2]);
				}
				counts.put(l[2], l[0]);
			}else{
				//print
				System.out.print(kmer);
				for (String sample : samples) {
					if(counts.containsKey(sample)){
						System.out.print(sep+counts.get(sample));
					}else{
						System.out.print(zero);
					}
				}
				for(int i=samples.size();i<size;++i){
					System.out.print(zero);
				}
				System.out.println();
				counts.clear();
				kmer=l[1];
				if(!samples.contains(l[2])){
					samples.add(l[2]);
				}
				counts.put(l[2], l[0]);
				if(samples.size()==size){
					break;
				}
			}
		}
		//main loop
		for(s=in.readLine();s!=null;s=in.readLine()){
			l=s.split("\t");
			if(kmer.equals(l[1])){
				counts.put(l[2], l[0]);
			}else{
				//print
				System.out.print(kmer);
				for (String sample : samples) {
					if(counts.containsKey(sample)){
						System.out.print(sep+counts.get(sample));
					}else{
						System.out.print(zero);
					}
				}
				System.out.println();
				counts.clear();
				kmer=l[1];
				counts.put(l[2], l[0]);
			}
		}
		//flush the last kmer
		System.out.print(kmer);
		for (String sample : samples) {
			if(counts.containsKey(sample)){
				System.out.print(sep+counts.get(sample));
			}else{
				System.out.print(zero);
			}
		}
		System.out.println();
	}

	public static void multiKmerMergeJellyfishPipe(final int size)throws Exception{
		BufferedReader in = new BufferedReader(new InputStreamReader(System.in));
		HashMap<String, String[]> kmers= new HashMap<String, String[]>();

		final String zero=sep+"0";
		String s=in.readLine();
		String[] l=s.split("\t"),pair= new String[size];
		String cur=l[1];
		pair[Integer.parseInt(l[3])]=l[0];
		kmers.put(l[2], pair);
		
		for(s=in.readLine();s!=null;s=in.readLine()){
			l=s.split("\t");
			if(cur.equals(l[1])){
				if(kmers.containsKey(l[2])){
					kmers.get(l[2])[Integer.parseInt(l[3])]=l[0];
				}else{
					pair= new String[size];
					pair[Integer.parseInt(l[3])]=l[0];
					kmers.put(l[2], pair);
				}
			}else{
				//print
				for(String kmer : kmers.keySet()){
					pair=kmers.get(kmer);
					System.out.print(kmer);
					for(int i=0;i<size;++i){
						if(pair[i]==null){
							System.out.print(zero);
						}else{
							System.out.print(sep+pair[i]);
						}
					}
					System.out.println();
				}
				//restart
				kmers.clear();
				cur=l[1];
				pair= new String[size];
				pair[Integer.parseInt(l[3])]=l[0];
				kmers.put(l[2], pair);
			}
		}
		for(String kmer : kmers.keySet()){
			pair=kmers.get(kmer);
			System.out.print(kmer);
			for(int i=0;i<size;++i){
				if(pair[i]==null){
					System.out.print(zero);
				}else{
					System.out.print(sep+pair[i]);
				}
			}
			System.out.println();

		}
	}
	
	public static void multiKmerMergeJellyfishPipe2(final int size)throws Exception{
		multiKmerMergeJellyfishPipe2(size,1000);
	}
	
	public static void multiKmerMergeJellyfishPipe2(final int size,final int maximumCollisions)throws Exception{
		if(size>10){
			throw new Exception("size's above 10 is not yet supported");
		}
		BufferedReader in = new BufferedReader(new InputStreamReader(System.in));
		
//		final String zero="0";
		String[][] buffer= new String[maximumCollisions][size];

		
//		HashMap<String,Integer> kmers= new HashMap<String, Integer>(maximumCollisions*2);
		String[] kmers= new String[maximumCollisions];
		
		String s=in.readLine();
		String[] l=s.split("\t");
		String cur=l[1];
		int curSize=1,curRow;
		kmers[0]=l[2];
		for(int i=0;i<size;++i){
			buffer[0][i]="0";
		}
		buffer[0][multiKmerMergeJellyfishPipe2_parseInt(l[3])]=l[0];
		
//		kmers.put(l[2], pair);

		for(s=in.readLine();s!=null;s=in.readLine()){
			l=s.split("\t");
			if(cur.equals(l[1])){
				curRow=-1;
				for(int i=0;i<curSize&&curRow==-1;++i){
					if(l[2].equals(kmers[i])){
						curRow=i;
					}
				}
				if(curRow==-1){
					//add new
					for(int i=0;i<size;++i){
						buffer[curSize][i]="0";
					}
					buffer[curSize][multiKmerMergeJellyfishPipe2_parseInt(l[3])]=l[0];
					kmers[curSize]=l[2];
					++curSize;
				}else{
					buffer[curRow][multiKmerMergeJellyfishPipe2_parseInt(l[3])]=l[0];
				}
			}else{
				//print
				for(int i=0;i<curSize;++i){
					System.out.print(kmers[i]);
					for(int j=0;j<size;++j){
						System.out.print(sep+buffer[i][j]);
					}
					System.out.println();
				}
				System.out.flush();
				
				//restart
				cur=l[1];
				curSize=1;
				for(int i=0;i<size;++i){
					buffer[0][i]="0";
				}
				buffer[0][multiKmerMergeJellyfishPipe2_parseInt(l[3])]=l[0];
				kmers[0]=l[2];
			}
		}
		for(int i=0;i<curSize;++i){
			System.out.print(kmers[i]);
			for(int j=0;j<size;++j){
				System.out.print(sep+buffer[i][j]);
			}
			System.out.println();
		}
	}	
	
	public static int multiKmerMergeJellyfishPipe2_parseInt(String s)throws Exception{
//		if(s.length()==1){
			switch (s.charAt(0)) {
			case '0':
				return 0;
			case '1':
				return 1;
			case '2':
				return 2;
			case '3':
				return 3;
			case '4':
				return 4;
			case '5':
				return 5;
			case '6':
				return 6;
			case '7':
				return 7;
			case '8':
				return 8;
			case '9':
				return 9;
			default:
				throw new Exception(s+" is not a number");
			}
//		}else{
//			throw new Exception(s+" is an unknown column value");
//		}
	}

	
	public static void checkEnds(String kmerFile1, String kmerFile2, final int length)throws Exception{
		HashMap<String, ArrayList<String[]>> starts= new HashMap<String, ArrayList<String[]>>(),ends= new HashMap<String, ArrayList<String[]>>();
		ArrayList<String[]> alstmp;
		fastaParser fp= new fastaParser(new BufferedReader(new FileReader(kmerFile1)));
		FastaSeq fs;
		String start, end;
		final int sizeLimit=100;
		
		
//		String kmer=s.split("\t")[pos];
//		
//		//prep the hashmaps with the data from the first file
//		final int length=kmer.length();
//		final int endStart=length-n;
//		start=kmer.substring(0, n);
//		end=kmer.substring(endStart);
//		alstmp=new ArrayList<String[]>();
//		alstmp.add(new String[] {end,kmer});
//		starts.put(start, alstmp);
//		alstmp=new ArrayList<String[]>();
//		alstmp.add(new String[] {start,kmer});
//		ends.put(end, alstmp);
		System.err.println("Parsing first file: "+kmerFile1);
		for(;fp.hasNext();){
			fs=fp.next();
			start=kmerFunctions.kmerToUse(fs.getSeq(0, length));
			end=kmerFunctions.kmerToUse(fs.getSeq(fs.length()-length, fs.length()));
			if(starts.containsKey(start)){
				alstmp=starts.get(start);
			}else{
				alstmp= new ArrayList<String[]>();
				starts.put(start, alstmp);
			}
			alstmp.add(new String[] {end,fs.getQname()});
			
			if(ends.containsKey(end)){
				alstmp=ends.get(end);
			}else{
				alstmp= new ArrayList<String[]>();
				ends.put(end, alstmp);
			}
			alstmp.add(new String[] {start,fs.getQname()});
		}
		
		//Analyse the data in the second file with respect to the hashmaps
		fp=new fastaParser(new BufferedReader(new FileReader(kmerFile2)));
		System.err.println("Parsing second file: "+kmerFile2);
		for(;fp.hasNext();){
			fs=fp.next();
			start=kmerFunctions.kmerToUse(fs.getSeq(0, length));
			end=kmerFunctions.kmerToUse(fs.getSeq(fs.length()-length, fs.length()));
			
			if(starts.containsKey(start)){
				alstmp=starts.get(start);
				if (alstmp.size()<sizeLimit) {
					for (String[] strings : alstmp) {
						System.out.println(strings[1] + sep + fs.getQname()
								+ sep + "ss");
					}
				}
			}
			if(starts.containsKey(end)){
				alstmp=starts.get(end);
				if (alstmp.size()<sizeLimit) {
					for (String[] strings : alstmp) {
						System.out.println(strings[1] + sep + fs.getQname()
								+ sep + "se");
					}
				}
			}
			if(ends.containsKey(start)){
				alstmp=ends.get(start);
				if (alstmp.size()<sizeLimit) {
					for (String[] strings : alstmp) {
						System.out.println(strings[1] + sep + fs.getQname()
								+ sep + "es");
					}
				}
			}
			if(ends.containsKey(end)){
				alstmp=ends.get(end);
				if (alstmp.size()<sizeLimit) {
					for (String[] strings : alstmp) {
						System.out.println(strings[1] + sep + fs.getQname()
								+ sep + "ee");
					}
				}
			}
		}
		System.err.println("Done!");
	}
	
	public static void extractGood(KmerSet_binary kmers,String fastqFile, int min, String outPrefix)throws Exception{
		BufferedWriter fastqout= new BufferedWriter(new FileWriter(outPrefix+"_cleaned.fastq"));
		BufferedWriter countout= new BufferedWriter(new FileWriter(outPrefix+"_count.csv"));
		fastqParser fqp= new fastqParser(new BufferedReader(new FileReader(fastqFile)), "");
		FastqSeq fqs;
//		int kmerSize=kmers.getKmerSize();
		int count;
		System.err.println("Parsing fastq file...");
		for(int i=0;fqp.hasNext();++i){
			if(i%10000==0){
				System.err.print("    "+i+"           \r");
			}
			fqs=fqp.next();
			count=kmers.count(fqs.getSeq());
//			for (String kmer : kmerizeString(fqs.getSeq(),kmerSize)) {
//				if(kmers.exists(kmer)){
//					++count;
//				}
//			}
			if(count>=min){
				fastqout.write(fqs.toString()+"\n");
				countout.write(fqs.getQname()+"\t"+count+"\n");
			}
		}
		
		fastqout.close();
		countout.close();
	}
	
//	private static ArrayList<String> kmerizeString(String seq, int length){
//		ArrayList<String> seqs = new ArrayList<String>();
//		if(seq.length()>=length){
//			String kmer=seq.substring(0, length);
//			seqs.add(kmerToUse(kmer));
//			for(int i=length;i<seq.length();i++){
//				kmer=kmer.substring(1)+seq.charAt(i);
//				seqs.add(kmerToUse(kmer));
//			}
//		}
//		return seqs;
//	}
	
	private static KmerSet_binary readKmers_revComp(String kmerCountFile, int pos)throws Exception{
		BufferedReader in= new BufferedReader(new FileReader(kmerCountFile));
		String[] l=in.readLine().split("\t");
		KmerSet_binary kmers= new KmerSet_binary(l[pos]);
		int line=1;
		for(String s=in.readLine();s!=null;s=in.readLine(),line++){
			if(line%10000==0)
				System.err.print("\t"+line+"   "+kmers.size()+"\r");
			l=s.split("\t");
			kmers.addSeq(l[pos]);
			kmers.addSeq(reverseComplementSeq(l[pos]));
		}
		
		return kmers;
	}
	
	private static KmerSet_binary readKmers(String kmerCountFile, int pos)throws Exception{
		BufferedReader in= new BufferedReader(new FileReader(kmerCountFile));
		String[] l=in.readLine().split("\t");
		KmerSet_binary kmers= new KmerSet_binary(l[pos]);
		int line=1;
		for(String s=in.readLine();s!=null;s=in.readLine(),line++){
			if(line%10000==0)
				System.err.print("\t"+line+"   "+kmers.size()+"\r");
			l=s.split("\t");
			kmers.addSeq(l[pos]);
		}
		
		return kmers;
	}
	
	public static void reduce(String kmerFile,int pos,int length)throws Exception{
		BufferedReader in= new BufferedReader(new FileReader(kmerFile));
		for(String s=in.readLine();s!=null;s=in.readLine()){
			System.out.println((new kmerData(pos, s)).reduce(length));
		}
	}
	
	private static String printHelp(){
		String help="";
		ArrayList<String> cmds= new ArrayList<String>(methodsHelp.keySet());
		Collections.sort(cmds);
		for (String s : cmds) {
			help+=methodsHelp.get(s);
		}
		
		return help;
	}
	
	private static String printHelp(String method){
		if(methodsHelp.containsKey(method)){
			return methodsHelp.get(method);
		}
		return printHelp();
	}
	
	protected static String reverseComplementSeq(String seq){
		String revCompSeq="";
		for(int i=0;i<seq.length();i++){
			revCompSeq= complement(seq.charAt(i))+revCompSeq;
		}
		return revCompSeq;
	}
	
	protected static String reverse(String seq){
		String rev="";
		for(int i=0;i<seq.length();i++){
			rev=seq.charAt(i)+rev;
		}
		return rev;
	}
	
	protected static String complement(String seq){
		String comp="";
		for(int i=0;i<seq.length();i++){
			comp=comp+complement(seq.charAt(i));
		}
		return comp;
	}
	
	private static char complement(char c){
		switch (c) {
		case 'A':
			return 'T';
		case 'T':
			return 'A';
		case 'G':
			return 'C';
		case 'C':
			return 'G';
		case 'N':
			return 'N';
		default:
			System.err.println("unknown nucleotide: "+c+", will return "+c);
			return c;
		}
	}
	
//	protected static String kmerToUse(String kmer){
//		String revKmer=reverseComplementSeq(kmer);
//		if(kmer.compareTo(revKmer)<0){
//			return kmer;
//		}else{
//			return revKmer;
//		}
//	}
	
	public static String kmerToUse(String kmer){
		String kmerRev= "";
		int j=kmer.length()-1,result;
		char c,r;
		for(int i=0;i<kmer.length();i++,j--){
			c=kmer.charAt(i);
			r=complement(kmer.charAt(j));
			if((result=c-r)!=0){
				if(result<0){
					return kmer;
				}else{
					kmerRev+=""+r;
					++i;
					break;
				}
			}else{
				kmerRev=r+kmerRev;
			}
		}
		for(;j>=0;j--){
			kmerRev+=complement(kmer.charAt(j))+"";
		}
		return kmerRev;
	}
	
	private static void printKmer(String kmer){
		System.out.println(kmerToUse(kmer));
	}
	
	public static void kmerize(String fastqFile, int n)throws Exception{
		fastqParser fqp=new fastqParser(new BufferedReader(new FileReader(fastqFile)),getPrefix(fastqFile));
		String seq;
		for(;fqp.hasNext();){
			seq=fqp.next().getSeq();
			for(int i=0,j=n;j<=seq.length();i++,j++){
				printKmer(seq.substring(i, j));
			}
		}
	}
	
	public static void kmerizeEnd(String fastqFile, int n, int m)throws Exception{
		if(m>=n){
			throw new Exception("m ("+m+") must be smaller than n ("+n+")");
		}
		fastqParser fqp=new fastqParser(new BufferedReader(new FileReader(fastqFile)),getPrefix(fastqFile));
		for(;fqp.hasNext();){
			String seq= fqp.next().getSeq();
			seq=seq.substring(seq.length()-n+1);
			for(int i=0,j=m;j<=seq.length();i++,j++){
				printKmer(seq.substring(i, j));
			}
		}
	}
	
//	private static void mapListMutationKmers(String mapListFile,String mutationFile,int n, String outPrefix) throws Exception{
//		BufferedReader mutationReader = new BufferedReader(new FileReader(mutationFile));
//		HashMap<String, ArrayList<Mutation>> mutations= new HashMap<String, ArrayList<Mutation>>();
//		Mutation mutation;
//		for(String s=mutationReader.readLine();s!=null;s=mutationReader.readLine()){
//			String l[]=s.split("\t");
//			mutation=new Mutation(l[0], Integer.parseInt(l[1]), l[2].charAt(0));
//			if(!mutations.containsKey(mutation.getChr())){
//				mutations.put(mutation.getChr(), new ArrayList<Mutation>());
//			}
//			mutations.get(mutation.getChr()).add(mutation);
//		}
//		mutationReader.close();
//		BufferedReader mapListReader = new BufferedReader(new FileReader(mapListFile));
//		ShoreMapLine sml;
//		BufferedWriter origWriter= new BufferedWriter(new FileWriter(outPrefix+"_orig.kmers"));
//		BufferedWriter mutatedWriter= new BufferedWriter(new FileWriter(outPrefix+"_mutated.kmers"));
//		for(String s= mapListReader.readLine();s!=null;s=mapListReader.readLine()){
//			sml= new ShoreMapLine(s);
//			if(mutations.containsKey(sml.getChr())){
//				//check for mutation... this has to be rewritten if there are a lot of mutations or an unfiltered map.list
//				String orig=sml.getSeq();
//				String mutated=sml.getMutatedSeq(mutations.get(sml.getChr()));
//				if(!orig.equals(mutated)){
//					//kmerize
//					for(int i=0,j=n;j<orig.length();i++,j++){
//						origWriter.write(kmerToUse(orig.substring(i,j))+"\n");
//						mutatedWriter.write(kmerToUse(mutated.substring(i, j))+"\n");
//					}
//				}
//			}
//		}
//		mapListReader.close();
//		origWriter.close();
//		mutatedWriter.close();
//	}
	
	private static void generateKmers(int n,String s)throws Exception{
		if(n==0){
			System.out.println(s);
		}else{
			for(int i=0;i<4;i++){
				generateKmers(n-1, nucleotides[i]+s);
			}
		}
	}
	
	private static void generateSeeds3(String kmerFile1, String kmerFile2, int pos, int minSeedSize, int lmer,String outPrefix)throws Exception{
		boolean verbose= false,added;
		BufferedReader in;
		HashMap<String,HashSet<kmerData>> kmers= new HashMap<String, HashSet<kmerData>>(50000000);
		HashSet<String> suffixes= new HashSet<String>(70000000);
		ArrayList<seedData> seeds= new ArrayList<seedData>();
		int counter=0;
		ArrayList<Character> alphabet= new ArrayList<Character>();
		String alp="ACGT";
		HashSet<String> kmersToCheck= new HashSet<String>();
		kmersToCheck.add("TCTCTTTCTGATTCTCGAAGAATGAGGGGCG");
		kmersToCheck.add("CGCCCCTCATTCTTCGAGAATCAGAAAGAGA");

		for(int i=0; i < alp.length();++i)
		{
			alphabet.add(alp.charAt(i));
		}
//		Entropy ent = new Entropy(alphabet);
		
		HashSet<kmerData> putativeStartPoints= new HashSet<kmerData>(),trueStartPoints= new HashSet<kmerData>();
		
		System.err.println("reading background...");
		kmerData kmer;
		
		KmerSet_binary lmers=new KmerSet_binary(lmer, 10000000);
		
		in= new BufferedReader(new FileReader(kmerFile2));
//		int count=0;
		for(String s= in.readLine();s!=null;s=in.readLine(),++counter){
			if(counter%10000==0){
				System.err.print("         "+counter+"\r");
			}
			kmer=new kmerData(pos, s);
			lmers.addSeq(kmer.getKmer());
//			if(count>10)
//				System.exit(0);
		}
		System.err.println("         "+counter);
		
		System.err.println("reading dataset ");
		counter=0;
		in= new BufferedReader(new FileReader(kmerFile1));
		for(String s=in.readLine();s!=null;s=in.readLine(),++counter){
			if(counter%10000==0){
				System.err.print("         "+counter+"\r");
			}
			kmer=new kmerData(pos, s);
			verbose= kmersToCheck.contains(kmer.getKmer());
			if(verbose) System.err.println("#"+kmer.getKmer()+" exists...");
			if(!kmers.containsKey(kmer.prefix())){
				kmers.put(kmer.prefix(), new HashSet<kmerData>());
			}
			if(!kmers.containsKey(kmer.prefixRevComp())){
				kmers.put(kmer.prefixRevComp(), new HashSet<kmerData>());
			}
			kmers.get(kmer.prefix()).add(kmer);
			kmers.get(kmer.prefixRevComp()).add(kmer);
			suffixes.add(kmer.suffix());
			suffixes.add(kmer.suffixRevComp());
			
			//find kmers that contains an lmer from the other sample
			if(lmers.covered(kmer.getKmer())){
				if(verbose) System.err.println("#"+kmer.getKmer()+" added as possible start point");
				putativeStartPoints.add(kmer);
			}
//			System.exit(0);
		}
		lmers=new KmerSet_binary(lmer, false);
		System.err.println("         "+counter);
		verbose=false;
//		System.err.println(putativeStartPoints.size());
		
		System.err.println("finding startpoints");
		counter=0;
		//find the startPoints, which have unique start points
		for (kmerData kmerData : putativeStartPoints) {
			++counter;
			if(counter%10000==0){
				System.err.print("         "+counter+"\r");
			}
			verbose= kmersToCheck.contains(kmerData.getKmer());
			if(verbose) System.err.println("#"+kmerData.getKmer()+" is start");
			added=false;
			if(!suffixes.contains(kmerData.prefix())){
				if(verbose) System.err.println("#"+kmerData.getKmer()+" added");
				seeds.add(new seedData(kmerData));
				added=true;
			}
			if(!suffixes.contains(kmerData.prefixRevComp())){
				if(verbose) System.err.println("#"+kmerData.getKmer()+" added rev comp");
				seeds.add(new seedData(kmerData.revComp()));
				added=true;
			}
			if(added){
				trueStartPoints.add(kmerData);
//				System.out.println(kmerData.getKmer()+sep+ent.calcEntropy(kmerData.getKmer()));
			}
		}
		System.err.println("         "+counter);
		verbose=false;
		
//		System.err.println(seeds.size()+"\t"+trueStartPoints.size()+"\t"+putativeStartPoints.size());
		putativeStartPoints.clear();
		suffixes.clear();
		
		System.err.println("generating seeds ("+seeds.size()+")");
		counter=0;
		int seedCount=0;
		BufferedWriter out= new BufferedWriter(new FileWriter(outPrefix+".fa"));
		BufferedWriter outLong= new BufferedWriter(new FileWriter(outPrefix+"_long.txt"));
//		HashSet<kmerData> usedKmers=new HashSet<kmerData>();
		
		BufferedWriter outSeeds= new BufferedWriter(new FileWriter(outPrefix+"_initialSeeds.txt"));
		
		for (seedData initialSeed : seeds) {
			outSeeds.write(initialSeed.consensus()+"\n");
			++counter;
			if(counter%1000==0){
				System.err.print("         "+counter+"\r");
			}
			verbose= kmersToCheck.contains(initialSeed.consensus());
			if(verbose) System.err.println("#"+initialSeed.consensus()+" extending");
			ArrayList<seedData> finalSeeds= generateSeeds_extend(initialSeed, kmers, trueStartPoints,kmersToCheck);
			final int batchSize=finalSeeds.size();
			for (seedData finalSeed : finalSeeds) {
				if(verbose) System.err.println("#"+initialSeed.consensus()+" "+finalSeed.consensus());
//				usedKmers.addAll(finalSeed.getKmers());
				if(finalSeed.size()>=minSeedSize){
					seedCount++;
					out.write(">"+seedCount+" "+finalSeed.size()+" "+batchSize+"\n"+finalSeed.consensus()+"\n");
					outLong.write(">"+seedCount+" "+finalSeed.size()+" "+batchSize+"\n"+finalSeed.toString()+"\n");
				}
			}
		}
		System.err.println("         "+counter);
		verbose=false;
		out.close();
		outLong.close();
		outSeeds.close();
		
//		BufferedWriter outKmer= new BufferedWriter(new FileWriter(outPrefix+"_used.kmerDiff"));
//		for (kmerData kmerData : usedKmers) {
//			outKmer.write(kmerData.toString()+"\n");
//		}
		
//		outKmer.close();
		System.err.println("Done!");

		
	}
	
	private static void generateSeeds2(String kmerFile1,String kmerFile2,int pos,int minSeedSize,String outPrefix)throws Exception{
		BufferedReader in;
		HashMap<String,HashSet<kmerData>> kmers= new HashMap<String, HashSet<kmerData>>();
		HashSet<String> suffixes= new HashSet<String>(1000000), prefixesOther= new HashSet<String>();
		ArrayList<seedData> seeds= new ArrayList<seedData>();
		
		HashSet<kmerData> startPoints= new HashSet<kmerData>();
		
		kmerData kmer;
		
		in=new BufferedReader(new FileReader(kmerFile2));
		for(String s= in.readLine();s!=null;s=in.readLine()){
			kmer=new kmerData(pos, s);
			prefixesOther.add(kmer.prefix());
			prefixesOther.add(kmer.prefixRevComp());
			prefixesOther.add(kmer.suffix());
			prefixesOther.add(kmer.suffixRevComp());
		}
		
		in= new BufferedReader(new FileReader(kmerFile1));
		for(String s=in.readLine();s!=null;s=in.readLine()){
			kmer=new kmerData(pos, s);
			if(!kmers.containsKey(kmer.prefix())){
				kmers.put(kmer.prefix(), new HashSet<kmerData>());
			}
			if(!kmers.containsKey(kmer.prefixRevComp())){
				kmers.put(kmer.prefixRevComp(), new HashSet<kmerData>());
			}
			kmers.get(kmer.prefix()).add(kmer);
			kmers.get(kmer.prefixRevComp()).add(kmer);
			suffixes.add(kmer.suffix());
			suffixes.add(kmer.suffixRevComp());
			if(prefixesOther.contains(kmer.prefix())){
				startPoints.add(kmer);
			}
			if(prefixesOther.contains(kmer.prefixRevComp())){
				startPoints.add(kmer.revComp());
			}
//			startPoints.add(kmer);
		}
		
		prefixesOther.clear();
		//Find starting points
//		in=new BufferedReader(new FileReader(kmerFile2));
//		for(String s= in.readLine();s!=null;s=in.readLine()){
//			kmer=new kmerData(pos, s);
//			if(kmers.containsKey(kmer.prefix())){
//				startPoints.addAll(kmers.get(kmer.prefix()));
//				
//			}
//			if(kmers.containsKey(kmer.prefixRevComp())){
//				for (kmerData kmerData : kmers.get(kmer.prefixRevComp())) {
//					startPoints.add(kmerData.revComp());
//				}
//				startPoints.addAll(kmers.get(kmer.prefixRevComp()));
//			}
//		}
		
		for (kmerData kmerData : startPoints) {
			seeds.add(new seedData(kmerData));
		}
		
//		Collections.sort(kmers);
		//Find seeds
		int seedCount=0;
		BufferedWriter out= new BufferedWriter(new FileWriter(outPrefix+".fa"));
		BufferedWriter outLong= new BufferedWriter(new FileWriter(outPrefix+"_long.txt"));
		HashSet<kmerData> usedKmers=new HashSet<kmerData>();
		
		for (seedData initialSeed : seeds) {
//			boolean print=(initialSeed.consensus().equals("CACCCTTACTCTCAGCTTCAACAAGTTTTTA")||initialSeed.consensus().equals("GATTAGCGCGGGTCTTCACGTAATAATCAGC"));
//			if(print) System.err.println(initialSeed.consensus());
			for (seedData finalSeed : generateSeeds_extend(initialSeed, kmers, startPoints)) {
//				if(print) System.err.println(finalSeed.toString());
				usedKmers.addAll(finalSeed.getKmers());
				if(finalSeed.size()>=minSeedSize){
					seedCount++;
					out.write(">"+seedCount+" "+finalSeed.size()+"\n"+finalSeed.consensus()+"\n");
					outLong.write(">"+seedCount+" "+finalSeed.size()+"\n"+finalSeed.toString()+"\n");
				}
			}
		}
		out.close();
		outLong.close();
		
		BufferedWriter outKmer= new BufferedWriter(new FileWriter(outPrefix+"_used.kmerDiff"));
		for (kmerData kmerData : usedKmers) {
			outKmer.write(kmerData.toString()+"\n");
		}
		
		outKmer.close();
//		return seeds;
		System.err.println("Done!");
	}
	
	private static ArrayList<seedData> generateSeeds_extend(seedData seed, HashMap<String,HashSet<kmerData>> kmers, HashSet<kmerData> endPoints){
		return generateSeeds_extend(seed,kmers,endPoints, new HashSet<String>());
	}
	
	private static ArrayList<seedData> generateSeeds_extend(seedData seed, HashMap<String,HashSet<kmerData>> kmers, HashSet<kmerData> endPoints,HashSet<String> kmersToCheck){
		ArrayList<seedData> seeds= new ArrayList<seedData>(),nextSeeds, doneSeeds= new ArrayList<seedData>();
		boolean fin;
		seedData nextSeed;
//		Entropy ent= new Entropy();
//		String broken="complete";
		
		seeds.add(seed);
		boolean print=kmersToCheck.contains(seed.consensus());
		String kp="";
		if(print) kp="#"+seed.consensus()+" ";
		int i;
		
		for(i=0;seeds.size()>0;++i){
			nextSeeds= new ArrayList<seedData>();
			for(seedData curSeed: seeds){
				//try to extend
//				if(print)System.err.println(kp+curSeed.suffix());
				if(kmers.containsKey(curSeed.suffix())){
					if(print)System.err.println(kp+"extend "+kmers.get(curSeed.suffix()).size());
					fin=true;
					for (kmerData kmer : kmers.get(curSeed.suffix())) {
						nextSeed= new seedData(curSeed);
						if(nextSeed.addRight(kmer)){
							if(endPoints.contains(kmer)){
								if(print)System.err.println(kp+"found endpoint");
								doneSeeds.add(nextSeed);
							}else{
								if(print)System.err.println(kp+"continue "+nextSeed.consensus());
								nextSeeds.add(nextSeed);
							}
							fin=false;
						}
					}
					if(fin){
						if(print)System.err.println(kp+"no further path");
						doneSeeds.add(curSeed);
					}
				}else{
					if(print)System.err.println(kp+"no more extension");
					doneSeeds.add(curSeed);
				}
			}
			if(i>150 || nextSeeds.size()>250){ //length cutoff, size cutoff
				if(print)System.err.println(kp+"length cutoff");
				doneSeeds.addAll(nextSeeds);
				nextSeeds.clear();
			}
			seeds.clear();
			seeds.addAll(nextSeeds);
		}
//		System.out.println(broken+"\t"+seed.getKmers().get(0).getKmer()+"\t"+doneSeeds.size()+"\t"+i+"\t"+ent.calcEntropy(seed.getKmers().get(0).getKmer()));
		
		return doneSeeds;
	}
	
	private static void generateSeeds(String kmerFile,int pos,int minSeedSize,String outPrefix,boolean longOutput)throws Exception{
		BufferedReader in= new BufferedReader(new FileReader(kmerFile));
		HashMap<String,HashSet<kmerData>> kmers= new HashMap<String, HashSet<kmerData>>(60000000);
		HashSet<String> suffixes= new HashSet<String>(1000000);
		
		HashSet<kmerData> uniqueKmers= new HashSet<kmerData>();
		int counter=0;
		kmerData kmer;
		System.err.println("loading data...");
		for(String s=in.readLine();s!=null;s=in.readLine(),++counter){
			if(counter%10000==0) System.err.print("      "+counter+"\r");
			kmer=new kmerData(pos, s);
			if(!kmers.containsKey(kmer.prefix())){
				kmers.put(kmer.prefix(), new HashSet<kmerData>());
			}
			if(!kmers.containsKey(kmer.prefixRevComp())){
				kmers.put(kmer.prefixRevComp(), new HashSet<kmerData>());
			}
			kmers.get(kmer.prefix()).add(kmer);
			kmers.get(kmer.prefixRevComp()).add(kmer);
			suffixes.add(kmer.suffix());
			suffixes.add(kmer.suffixRevComp());
			uniqueKmers.add(kmer);
		}
		System.err.print("      "+counter+"\n");
		
//		Collections.sort(kmers);
		//Find seeds
		int seedCount=0;
		ArrayList<seedData> seeds= new ArrayList<seedData>();
		BufferedWriter out= new BufferedWriter(new FileWriter(outPrefix+".fa"));
		BufferedWriter outLong= new BufferedWriter(new FileWriter(outPrefix+"_long.txt"));
		HashSet<kmerData> usedKmers;
		boolean found=true;
		System.err.println("Unique seeds...");
		while(found){
			found=false;			
//			System.err.print("    "+kmers.size()+"\r");
			//find unique seeds
			System.err.println("finding startpoints...");
			counter=0;
			seeds=new ArrayList<seedData>();
			for (kmerData kmerData : uniqueKmers) {
				if((++counter)%10000==0) System.err.print("     "+counter+"\r");
				if(!suffixes.contains(kmerData.prefix())){
					seeds.add(new seedData(kmerData));
				}
				if(!suffixes.contains(kmerData.prefixRevComp())){
					seeds.add(new seedData(kmerData.revComp()));
				}
			}
			System.err.print("     "+counter+"\n");
//			System.err.print("    "+kmers.size()+" "+seeds.size()+"      "+"\r");
			found=seeds.size()>0;
			//extend seeds
//			usedKmers=new HashSet<kmerData>();
			System.err.println("extending seeds ("+seeds.size()+")...");
			counter=0;
			for (seedData initialSeed : seeds) {
				if((++counter)%1000==0) System.err.print("     "+counter+"\r");
				ArrayList<seedData> finalSeeds= generateSeeds_extend(initialSeed, kmers);
				final int batchSize= finalSeeds.size();
				for (seedData finalSeed : finalSeeds) {
//					usedKmers.addAll(finalSeed.getKmers());
					if(finalSeed.size()>=minSeedSize){
						seedCount++;
						out.write(">"+seedCount+" "+finalSeed.size()+" "+batchSize+"\n"+finalSeed.consensus()+"\n");
						if(longOutput) outLong.write(">"+seedCount+" "+finalSeed.size()+" "+batchSize+"\n"+finalSeed.toString()+"\n");
					}
				}
			}
			System.err.print("     "+counter+"\n");
			//clean
//			uniqueKmers.removeAll(usedKmers);
			
//			for (kmerData kmerData : usedKmers) {
//				kmers.remove(kmerData.prefix());
//				kmers.remove(kmerData.prefixRevComp());
//			}
			//Rebuild the suffix list
//			suffixes= new HashSet<String>(1000000);
//			for (kmerData kmerData : uniqueKmers) {
//				suffixes.add(kmerData.suffix());
//				suffixes.add(kmerData.suffixRevComp());
//			}
			found=false;
			kmers.clear();
		}
		System.err.println("Repetitive structures...");
		while(kmers.size()>0){ // not used
			HashSet<kmerData> onePrefix= new HashSet<kmerData>();
			for (HashSet<kmerData> tmp : kmers.values()) {
				onePrefix=new HashSet<kmerData>(tmp);
				break;
			}
			usedKmers= new HashSet<kmerData>();
			for (kmerData kmerData : onePrefix) {
				for(seedData finalSeed : generateSeeds_extend(new seedData(kmerData), kmers)){
					usedKmers.addAll(finalSeed.getKmers());
					if(finalSeed.size()>=minSeedSize){
						seedCount++;
						out.write(">"+seedCount+" "+finalSeed.size()+" cyclic\n"+finalSeed.consensus()+"\n");
						if(longOutput) outLong.write(">"+seedCount+" "+finalSeed.size()+" cyclic\n"+finalSeed.toString()+"\n");
					}
				}
			}
			//clean
			uniqueKmers.removeAll(usedKmers);
			
			for (kmerData kmerData : usedKmers) {
				kmers.remove(kmerData.prefix());
				kmers.remove(kmerData.prefixRevComp());
			}
		}
		out.close();
		outLong.close();
//		return seeds;
		System.err.println("Done!");
	}
	
	private static ArrayList<seedData> generateSeeds_extend(seedData seed, HashMap<String,HashSet<kmerData>> kmers){
		return generateSeeds_extend(seed, kmers,new HashSet<kmerData>());
//		ArrayList<seedData> seeds= new ArrayList<seedData>(),nextSeeds, doneSeeds= new ArrayList<seedData>();
//		boolean fin;
//		seedData nextSeed;
//		
//		seeds.add(seed);
//		
//		for(int i=0;seeds.size()>0;++i){
//			if(i%10==0){
//				StringBuffer out=new StringBuffer(i+"  "+seeds.size()+"\r");
//				while(out.length()<20){
//					out.insert(0, ' ');
//				}
//				System.err.print(out);
//			}
//			nextSeeds= new ArrayList<seedData>();
//			for(seedData curSeed: seeds){
//				//try to extend
//				if(kmers.containsKey(curSeed.suffix())){
//					fin=true;
//					for (kmerData kmer : kmers.get(curSeed.suffix())) {
//						nextSeed= new seedData(curSeed);
//						if(nextSeed.addRight(kmer)){
//							nextSeeds.add(nextSeed);
//							fin=false;
//						}
//					}
//					if(fin){
//						doneSeeds.add(curSeed);
//					}
//				}else{
//					doneSeeds.add(curSeed);
//				}
//			}
//			if(i>150){
//				doneSeeds.addAll(nextSeeds);
//				nextSeeds.clear();
//			}
//			seeds.clear();
//			seeds.addAll(nextSeeds);
//		}
//		return doneSeeds;
	}
	
	
	private static void sortKmersLink(String kmerFile, int pos, int maxJump)throws Exception{
//		System.err.println("kmerFile:"+kmerFile+"\npos: "+pos+"\nmaxJump: "+maxJump);
		BufferedReader in= new BufferedReader(new FileReader(kmerFile));
		ArrayList<kmerData> kmers= new ArrayList<kmerData>();
//		System.err.println("Reading "+kmerFile);
//		int line=0;
		for(String s=in.readLine();s!=null;s=in.readLine()){
//			if(line%10000==0){
//				System.out.print("\t"+line+"\r");
//			}
			kmers.add(new kmerData(pos, s));
//			++line;
		}
//		System.err.println(line+"\nsorting...");
		Collections.sort(kmers);
//		System.err.println("printing and idenifying links");
		String lastKmer="XXXXXXXXXXX";
		while(lastKmer.length()<maxJump+1);{
			lastKmer+="X";
		}
		boolean linked;
		for (kmerData kmerData : kmers) {
			System.out.print(kmerData);
			linked=false;
			for(int i=1;i<=maxJump&&!linked;i++){
				linked=kmerData.getKmer().startsWith(lastKmer.substring(i));
			}
			if(linked){
				System.out.println(sep+"1");
			}else{
				System.out.println(sep+"0");
			}
			lastKmer=kmerData.getKmer();
		}
	}
	
	private static void sortKmers(String kmerFile,int pos, int binSize,String tmpPrefix)throws Exception{
		BufferedReader in= new BufferedReader(new FileReader(kmerFile));
		ArrayList<kmerData> kmers= new ArrayList<kmerData>();
		for(String s=in.readLine();s!=null;s=in.readLine()){
			kmers.add(new kmerData(pos, s));
			if(kmers.size()>=binSize){
				Collections.sort(kmers);
				//push to disk
				
			}
		}
		Collections.sort(kmers);
		//do mergeSort...
		
		for (kmerData kmerData : kmers) {
			System.out.println(kmerData);
		}
	}
	public static String getPrefix(String fastqFile)throws Exception{
		fastqParser fqp;
		FastqSeq fqs;
		String prefix="";
		boolean done=true;
		for (int i=6;i>-1&&done;i--){
			fqp= new fastqParser(new BufferedReader(new FileReader(fastqFile)),"");
			fqs=fqp.next();
			prefix=fqs.getQname().substring(0,i);
			done=true;
			for(int j=0;j<10&&fqp.hasNext();j++){
				fqs=fqp.next();
				if(!prefix.equals(fqs.getQname().subSequence(0, i))){
					done=false;
				}
			}
		}
		return prefix;
	}
}

class seedData{
	
	private ArrayList<kmerData> kmers;
	private ArrayList<Integer> strand;
	private String consensus;
	
	private seedData(){
		kmers= new ArrayList<kmerData>();
		strand= new ArrayList<Integer>();
		consensus="";
	}
	
	public seedData(kmerData first){
		this();
		kmers.add(first);
		strand.add(0);
		consensus=first.getKmer();
	}
	
	public seedData(seedData old){
		this();
		this.kmers.addAll(old.kmers);
		this.strand.addAll(old.strand);
		this.consensus=old.consensus;
	}
	
	public String suffix(int n){
//		System.err.println(n+";"+consensus.length());
		return consensus.substring(consensus.length()-n,consensus.length());
	}
	
	public String suffix(){
		return suffix(kmers.get(0).getKmer().length()-1);
	}
	
	public ArrayList<kmerData> getKmers(){
		return kmers;
	}
	
	public String consensus(){
		//presumes that the kmers are of the same length
//		String seq=kmers.get(0).getKmer();
//		int lastPos=seq.length()-1;
//		for(int i=1;i<kmers.size();i++){
//			seq+=""+kmers.get(i).getKmer().charAt(lastPos);
//		}
//		return seq;
		return consensus;
	}
	
	public boolean addRight(kmerData kd){
		if(this.contains(kd)){
			return false;
		}
		//same strand
		if(consensus.endsWith(kd.prefix())){
			kmers.add(kd);
			strand.add(0);
			consensus+=kd.tail();
			return true;
		}
		//opposite strand
		if(consensus.endsWith(kd.prefixRevComp())){
			kmers.add(kd);
			strand.add(1);
			consensus+=kmerUtils.complement(kd.head());
			return true;
		}
		return false;
	}
	
	public boolean addLeft(kmerData kd){
		if(this.contains(kd)){
			return false;
		}
		//same strand
		if(consensus.startsWith(kd.suffix())){
			kmers.add(0, kd);
			strand.add(0,0);
			consensus=kd.head()+consensus;
			return true;
		}
		//opposite strand
		if(consensus.startsWith(kd.suffixRevComp())){
			kmers.add(0, kd);
			strand.add(0,1);
			consensus=kmerUtils.complement(kd.tail())+consensus;
		}
		return false;
	}
	
//	public boolean addRight(kmerData ext){
//		if(extendRight(ext)&&!contains(ext)){
//			kmers.add(ext);
//			return true;
//		}
//		return false;
//	}
//	
//	public boolean addLeft(kmerData ext){
//		if(extendLeft(ext)&&!contains(ext)){
//			kmers.add(0, ext);
//			return true;
//		}
//		return false;
//	}
//	
//	public boolean extendRight(kmerData ext){
//		return kmers.get(kmers.size()-1).extendRight(ext);
//	}
//	
//	public boolean extendLeft(kmerData ext){
//		return kmers.get(0).extendLeft(ext);
//	}
	
	public boolean contains(kmerData kd){
		return kmers.contains(kd);
	}
	
	public boolean contains(seedData seed){
		return this.kmers.containsAll(seed.kmers);
	}
	
	public int size(){
		return kmers.size();
	}
	
	private StringBuffer toStringBuffer(){
		StringBuffer out=new StringBuffer(consensus());
		for (int i=0;i<kmers.size();i++){
			if(strand.get(i)==0){
				out.append("\n");
				out.append(kmers.get(i).toStringBuffer());
			}else{
				out.append("\n");
				out.append(kmers.get(i).toStringBufferRevComp());
			}
		}
		
//		String lastSuffix=kmers.get(0).suffix();
//		kmerData kmer;
//		for(int i=1;i<kmers.size();i++){
//			kmer=kmers.get(i);
//			if(kmer.getKmer().startsWith(lastSuffix)){
//				out+="\n"+kmer.toString();
//			}else{
//				out+="\n"+kmer.toStringRevComp();
//			}
//			lastSuffix=kmer.suffix();
//		}
		return out;
	}
	
	public String toString(){
		return this.toStringBuffer().toString();
	}
}

class kmerData implements Comparable<kmerData>,Serializable{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private int pos;
	private int permutations;
	private String rep;
	private String kmerSeqRevComp;
	private String[] data;
	
	public kmerData(int pos, String s) throws Exception{
		this(pos, s.split("\t"));
	}
	
	public kmerData(int pos, String[] data) throws Exception {
		super();
		this.pos = pos;
		this.data = data;
		kmerSeqRevComp=kmerUtils.reverseComplementSeq(data[pos]);
		rep=data[pos];
		permutations=0;
		String tmp=rep;
		for(int i=0;i<rep.length();i++){
			tmp=tmp.substring(1)+tmp.charAt(0);
			if(tmp.compareTo(rep)<=0){
				rep=tmp;
				permutations=0;
			}else{
				permutations++;
			}
			if(data[pos].equals(tmp)){
				break;
			}
		}
	}
	
	public kmerData revComp()throws Exception{
		return new kmerData(pos, this.toStringRevComp());
	}
	
	public kmerData reduce(int n) throws Exception{
		String[] newData=new String[data.length];
		for(int i=0;i<data.length;i++){
			newData[i]=data[i];
		}
		newData[pos]=kmerUtils.kmerToUse(data[pos]);
		return new kmerData(pos,newData);
	}
	
	public String suffixRevComp(){
		return suffixRevComp(1);
	}
	
	public String suffixRevComp(int n){
		return getKmerRevComp().substring(n);
	}
	
	public String prefixRevComp(int n){
		return getKmerRevComp().substring(0, getKmer().length()-n);
	}
	
	public String prefixRevComp(){
		return prefixRevComp(1);
	}
	
	public String tail(){
		return ""+getKmer().charAt(getKmer().length()-1);
	}
	
	public String head(){
		return ""+getKmer().charAt(0);
	}
	
	public String suffix(int n){
		return getKmer().substring(n);
	}
	
	public String suffix(){
		return suffix(1);
	}
	
	public String prefix(int n){
		return getKmer().substring(0, getKmer().length()-n);
	}
	
	public String prefix(){
		return prefix(1);
	}
	
	private boolean extendLeft(String prefix){
		return getKmer().startsWith(prefix)||getKmerRevComp().startsWith(prefix);
	}
	
	public boolean extendLeft(kmerData prefix,int n){
		return extendLeft(prefix.suffix(n));
//		return startsWith(prefix.getKmer().substring(n));
	}
	
	public boolean extendLeft(kmerData prefix){
		return extendLeft(prefix,1);
	}
	
	private boolean extendRight(String suffix){
		return getKmer().endsWith(suffix)||getKmerRevComp().endsWith(suffix);
	}
	
	public boolean extendRight(kmerData suffix, int n){
		return extendRight(suffix.prefix(n));
		//return endsWith(suffix.getKmer().substring(0, suffix.getKmer().length()-n));
	}
	
	public boolean extendRight(kmerData suffix){
		return extendRight(suffix,1);
	}
	
	public String getKmer(){
		return data[pos];
	}
	
	private String getKmerRevComp(){
		return kmerSeqRevComp;
	}

	public int compareTo(kmerData o) {
		if(this.rep.equals(o.rep))
			return this.permutations-o.permutations;
		return this.rep.compareTo(o.rep);
	}
	
	public StringBuffer toStringBuffer(){
		return this.toStringBuffer("\t");
	}
	
	public StringBuffer toStringBufferRevComp(){
		return this.toStringBufferRevComp("\t");
	}
	
	public String toString(){
		return this.toStringBuffer().toString();
	}
	
	public String toStringRevComp(){
		return toStringBufferRevComp().toString();
	}
	
	public String toStringRevComp(String sep){
		return this.toStringBufferRevComp(sep).toString();
	}
	
	
	public StringBuffer toStringBufferRevComp(String sep){
		StringBuffer s=new StringBuffer("");
		if(pos==0){
			s=new StringBuffer(getKmerRevComp());
		}else{
			s=new StringBuffer(data[0]);
		}
		for(int i=1;i<data.length;i++){
			if(i==pos){
				s.append(sep);
				s.append(getKmerRevComp());
			}else{
				s.append(sep);
				s.append(data[i]);
			}
		}
		return s;
	}
	
	public String toString(String sep){
		return this.toStringBuffer(sep).toString();
	}
	
	public StringBuffer toStringBuffer(String sep){
		StringBuffer s=new StringBuffer(data[0]);
		for(int i=1;i<data.length;i++){
			s.append(sep);
			s.append(data[i]);
		}
		return s;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + Arrays.hashCode(data);
		result = prime * result + pos;
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		kmerData other = (kmerData) obj;
		if (!Arrays.equals(data, other.data))
			return false;
		if (pos != other.pos)
			return false;
		return true;
	}
	
	
}

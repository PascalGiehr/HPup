package tools.fastq;

//import information.StringContent;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.EOFException;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;
import java.util.regex.MatchResult;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import tools.fasta.FastaSeq;
import tools.fasta.fastaParser;
//import tools.kmer.KmerMap_binary;
//import tools.kmer.KmerSet_binary;
import tools.kmer.KmerSet_binary_utils;
//import tools.kmer.KmerTable;
//import tools.kmer.Kmer_count_tuple;
//import tools.kmer.kmerFunctions;
//import tools.reptile.ReptileLine;
//import tools.rocheQual.RocheQualParser;
//import tools.rocheQual.RocheQualSeq;
import tools.sequences.sequenceUtils;
import tools.utils.general;

public class fastqUtils {

	private static final String sep="\t";
	static final int match_award      = 10;
	static final int mismatch_penalty = -5;
	static final int gap_penalty      = -5;
	static int score[][];
	static int pointer[][];
	
	public static void main(String[] args) throws Exception{
		if(args.length>0){
			if(args[0].equals("sub")&&args.length>2){
//				ArrayList<String> fastqFiles= new ArrayList<String>();
				ArrayList<BufferedReader> fastqFiles= new ArrayList<BufferedReader>();
				for(int i=1;i<args.length-1;i++){
					fastqFiles.add(new BufferedReader(new FileReader(args[i])));
				}
				sub(fastqFiles, new BufferedReader(new FileReader(args[args.length-1])), true);
			}else if(args[0].equals("subPipe")&&args.length>1){
//				ArrayList<String> fastqFiles= new ArrayList<String>();
				ArrayList<BufferedReader> fastqFiles= new ArrayList<BufferedReader>();
				for(int i=1;i<args.length;i++){
//					fastqFiles.add(new BufferedReader(new FileReader(args[i])));
					fastqFiles.add(general.getBufferedReader(args[i]));
				}
				sub(fastqFiles, new BufferedReader(new InputStreamReader(System.in)), true);
			}else if(args[0].equals("subZ")&&args.length>2){
//				ArrayList<String> fastqFiles= new ArrayList<String>();
				ArrayList<BufferedReader> fastqFiles= new ArrayList<BufferedReader>();
				for(int i=1;i<args.length-1;i++){
					fastqFiles.add(new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(args[i])))));
				}
				sub(fastqFiles, new BufferedReader(new FileReader(args[args.length-1])), true);
			}else if(args[0].equals("subToFile")&&args.length>2){
				ArrayList<String> fastqFiles= new ArrayList<String>();
				for(int i=1;i<args.length-1;i++){
					fastqFiles.add(args[i]);
				}
				subToFile(fastqFiles, new BufferedReader(new FileReader(args[args.length-1])), true);
			}else if(args[0].equals("antisub")&&args.length>2){
//				ArrayList<String> fastqFiles= new ArrayList<String>();
				ArrayList<BufferedReader> fastqFiles= new ArrayList<BufferedReader>();
				for(int i=1;i<args.length-1;i++){
					fastqFiles.add(new BufferedReader(new FileReader(args[i])));
				}
				sub(fastqFiles, new BufferedReader(new FileReader(args[args.length-1])), false);
			}else if(args[0].equals("antisubZ")&&args.length>2){
//				ArrayList<String> fastqFiles= new ArrayList<String>();
				ArrayList<BufferedReader> fastqFiles= new ArrayList<BufferedReader>();
				for(int i=1;i<args.length-1;i++){
					fastqFiles.add(new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(args[i])))));
				}
				sub(fastqFiles, new BufferedReader(new FileReader(args[args.length-1])), true);
			}else if(args[0].equals("subReg")&&args.length>2){
				subReg(args[1],args[2],true);
			}else if(args[0].equals("antisubReg")&&args.length>2){
				subReg(args[1],args[2],false);
//			}else if(args[0].equals("clean")&&args.length>6){
//				clean(args[1],args[2],Integer.parseInt(args[3]),Integer.parseInt(args[4]),args[5],Integer.parseInt(args[6]));
			}else if(args[0].equals("complement")&&args.length>1){
				complement(general.getBufferedReader(args[1]));
//			}else if(args[0].equals("conversionRateNOMe")&&args.length>1){
//				conversionRateNOMe(args[1], false);
//			}else if(args[0].equals("conversionRateNOMeZ")&&args.length>1){
//				conversionRateNOMe(args[1], true);
//			}else if(args[0].equals("randomsub")&&args.length>2){
//				randomsub(args[1],Integer.parseInt(args[2]));
			}else if(args[0].equals("reverseComplement")&&args.length>1){
				reverseComplement(general.getBufferedReader(args[1]));
			}else if(args[0].equals("removeWindowPaired")&&args.length>6){
				removeWindowPaired(args[1],args[2],Integer.parseInt(args[3]),Integer.parseInt(args[4]),args[5],args[6],false);
			}else if(args[0].equals("removeWindowPairedZ")&&args.length>6){
				removeWindowPaired(args[1],args[2],Integer.parseInt(args[3]),Integer.parseInt(args[4]),args[5],args[6],true);
			}else if(args[0].equals("removeShort")&&args.length>2){
				removeShort(args[1],Integer.parseInt(args[2]));
//			}else if(args[0].equals("removeShortVelvetPair")&&args.length>4){
//				removeShortVelvetPair(args[1],Integer.parseInt(args[2]),args[3],args[4]);
//			}else if(args[0].equals("convertNumberToPhredQual")&&args.length>1){
//				convertNumberToPhredQual(args[1]);
//			}else if(args[0].equals("randomVelvetPairedSubset")&&args.length>2){
//				randomVelvetPairedSubset(args[1],Double.parseDouble(args[2]));
			}else if(args[0].equals("getPaired")&&args.length>4){
				getPaired(general.getBufferedReader(args[1]), 
						general.getBufferedReader(args[2]), 
						general.getBufferedReader(args[3]),
						args[4]);
			}else if(args[0].equals("hpKmerTrim")&&args.length>4){
				ArrayList<String> read1= new ArrayList<String>(),
						read2= new ArrayList<String>();
				for(int i=4;i<args.length-1;i+=2){
					read1.add(args[i-1]);
					read2.add(args[i]);
				}
				hpKmerTrim(args[1],Integer.parseInt(args[2]),read1,read2,args[args.length-1],true,2);
			}else if(args[0].equals("hpSeqSync")&&args.length>6){
				
				ArrayList<String> read1= new ArrayList<String>(),
						read2= new ArrayList<String>();
				String outPrefix="";
				if(args.length%2==0){
					for(int i=6;i<args.length-1;i+=2){
						read1.add(args[i-1]);
						read2.add(args[i]);
					}
					outPrefix=args[args.length-1];
					//System.err.println("outprefix");
				}else{
					for(int i=6;i<args.length;i+=2){
						read1.add(args[i-1]);
						read2.add(args[i]);
					}
					outPrefix=read1.get(0);
					outPrefix=outPrefix.replace("_R1", "");
					outPrefix=outPrefix.replace("_val_1", "");
					outPrefix=outPrefix.replace(".fq.gz", "");
				}
				//hpSeqSync(Boolean.parseBoolean(args[1]), Integer.parseInt(args[2]), Integer.parseInt(args[3]), Integer.parseInt(args[4]), read1, read2, outPrefix);
				hpSeqSync(Boolean.parseBoolean(args[1]), Integer.parseInt(args[2]), Integer.parseInt(args[3]), Integer.parseInt(args[4]), read1, read2, outPrefix);
			}else if(args[0].equals("hpSeqRemerge")&&args.length>2){
				hpSeqRemerge(args[1],args[2]);
//			}else if(args[0].equals("kmerCount")&&args.length>2){
//				kmerCount(args[1],Integer.parseInt(args[2]));
			}else if(args[0].equals("indexHammingDist")&&args.length==3){
				indexHammingDist(general.getBufferedReader(args[1]), args[2]);
			}else if(args[0].equals("indexHammingFilter")&&args.length==5){
				indexHammingFilterSE(general.getBufferedReader(args[1]), args[2], Integer.parseInt(args[3]), args[4]);
			}else if(args[0].equals("indexHammingFilter")&&args.length==7){
				indexHammingFilterPE(general.getBufferedReader(args[1]), general.getBufferedReader(args[2]), args[3], Integer.parseInt(args[4]), args[5], args[6]);
//			}else if(args[0].equals("kmerCount_convert")&&args.length>1){
//				kmerCount_convert(args[1]);
//			}else if(args[0].equals("kmerCount_histogram")&&args.length>2){
//				kmerCount_histogram(args[1], Integer.parseInt(args[2]));
//			}else if(args[0].equals("kmerCount_merge")&&args.length>3){
//				ArrayList<String> files= new ArrayList<String>();
//				for(int i=1;i<args.length-1;i++){
//					files.add(args[i]);
//				}
//				kmerCount_merge(files,args[args.length-1]);
//			}else if(args[0].equals("kmerize")&&args.length>2){
//				kmerize(new BufferedReader(new FileReader(args[1])),Integer.parseInt(args[2]));
//			}else if(args[0].equals("kmerizePipe")&&args.length>1){
//				kmerize(new BufferedReader(new InputStreamReader(System.in)),Integer.parseInt(args[1]));
//			}else if(args[0].equals("kmerCoverageToKmerHits")&&args.length==3){
//				kmerCoverageToKmerHits(args[1], Integer.parseInt(args[2]));
//			}else if(args[0].equals("kmerCoverage")&&args.length>3){
//				ArrayList<String> fastqFiles=new ArrayList<String>();
//				for(int i=3;i<args.length;i++){
//					fastqFiles.add(args[i]);
//				}
//				kmerCoverage(fastqFiles, kmerFunctions.readKmers(args[1], Integer.parseInt(args[2])));
//			}else if(args[0].equals("kmerCoverageZ")&&args.length>3){
//				ArrayList<String> fastqFiles=new ArrayList<String>();
//				for(int i=3;i<args.length;i++){
//					fastqFiles.add(args[i]);
//				}
//				kmerCoverage(fastqFiles, kmerFunctions.readKmersZ(args[1], Integer.parseInt(args[2])));
//			}else if(args[0].equals("kmerCoverage_count")&&args.length>3){
//				ArrayList<String> fastqFiles=new ArrayList<String>();
//				for(int i=3;i<args.length;i++){
//					fastqFiles.add(args[i]);
//				}
//				kmerCoverage_count(fastqFiles, kmerFunctions.readKmers(args[1], Integer.parseInt(args[2])));
//			}else if(args[0].equals("kmerCoverage_stored")&&args.length>2){
//				ArrayList<String> fastqFiles=new ArrayList<String>();
//				for(int i=2;i<args.length;i++){
//					fastqFiles.add(args[i]);
//				}
//				kmerCoverage(fastqFiles, kmerFunctions.readKmers(args[1]));
//			}else if(args[0].equals("storeKmerSet")&&args.length>3){
//				storeKmerSet(args[1],Integer.parseInt(args[2]),args[3]);
//			}else if(args[0].equals("kmerToUse")&&args.length>1){
//				kmerToUse(args[1]);
//			}else if(args[0].equals("kmerCoverage_generateCorrection")&&args.length>2){
//				kmerCoverage_generateCorrection(args[1],Integer.parseInt(args[2]));
//			}else if(args[0].equals("kmerCoverage_extractCovered")&&args.length>5){
//				ArrayList<String> file1=new ArrayList<String>(), file2= new ArrayList<String>();
//				String kmerDistFile=args[args.length-3];
//				int kmerPos=Integer.parseInt(args[args.length-2]);
//				String outPrefix= args[args.length-1];
//				
//				for(int i=1;i<1+(args.length-4)/2;++i){
//					file1.add(args[i]);
//				}
//				for(int i=1+(args.length-4)/2;i<args.length-3;++i){
//					file2.add(args[i]);
//				}
//				
//				System.err.println(file1.size()+"\t"+file2.size());
//				
//				kmerCoverage_extractCovered(file1,file2,kmerDistFile,kmerPos,outPrefix,false);
//			}else if(args[0].equals("kmerCoverage_extractCoveredZ")&&args.length>5){
//				ArrayList<String> file1=new ArrayList<String>(), file2= new ArrayList<String>();
//				String kmerDistFile=args[args.length-3];
//				int kmerPos=Integer.parseInt(args[args.length-2]);
//				String outPrefix= args[args.length-1];
//				
//				for(int i=1;i<1+(args.length-4)/2;++i){
//					file1.add(args[i]);
//				}
//				for(int i=1+(args.length-4)/2;i<args.length-3;++i){
//					file2.add(args[i]);
//				}
//				
//				System.err.println(file1.size()+"\t"+file2.size());
//				kmerCoverage_extractCovered(file1,file2,kmerDistFile,kmerPos,outPrefix,true);
//			}else if(args[0].equals("kmerCoverage_extractCoveredSingle")&&args.length>4){
//				ArrayList<String> file1=new ArrayList<String>();
//				String kmerDistFile=args[args.length-3];
//				int kmerPos=Integer.parseInt(args[args.length-2]);
//				String outPrefix= args[args.length-1];
//				
//				for(int i=1;i<args.length-3;++i){
//					file1.add(args[i]);
//				}
//				
//				System.err.println(file1.size());
//				kmerCoverage_extractCoveredSingle(file1,kmerDistFile,kmerPos,outPrefix,false);
//			}else if(args[0].equals("kmerCoverage_extractCoveredSingleZ")&&args.length>4){
//				ArrayList<String> file1=new ArrayList<String>();
//				String kmerDistFile=args[args.length-3];
//				int kmerPos=Integer.parseInt(args[args.length-2]);
//				String outPrefix= args[args.length-1];
//				
//				for(int i=1;i<args.length-3;++i){
//					file1.add(args[i]);
//				}
//				
//				System.err.println(file1.size());
//				kmerCoverage_extractCoveredSingle(file1,kmerDistFile,kmerPos,outPrefix,false);
//			}else if(args[0].equals("kmerCoverage_extractPerfect")&&args.length>3){
//				kmerCoverage_extractPerfect(args[1],args[2],Integer.parseInt(args[3]));
//			}else if(args[0].equals("kmerCoverage_extractPerfectDirect")&&args.length>3){
//				ArrayList<String> fastqFiles= new ArrayList<String>();
//				for(int i=3;i<args.length;++i){
//					fastqFiles.add(args[i]);
//				}
//				kmerCoverage_extractPerfectDirect(kmerFunctions.readKmers(args[1], Integer.parseInt(args[2])),fastqFiles);
			}else if(args[0].equals("lengths")&&args.length>1){
				lengths(general.getBufferedReader(args[1]));
//			}else if(args[0].equals("reptileCorrection")&&args.length>2){
//				reptileCorrection(args[1],args[2]);
//			}else if(args[0].equals("normalizeKmer")&&args.length==6){
//				normalizeKmer(args[1],args[2], Integer.parseInt(args[3]),Integer.parseInt(args[4]),args[5]);
//			}else if(args[0].equals("normalizeKmer2")&&args.length==6){
//				normalizeKmer2(args[1],args[2], args[3],Integer.parseInt(args[4]),args[5]);
//			}else if(args[0].equals("cleanVelvetFile")&&args.length==2){
//				cleanVelvetFile(args[1]);
//			}else if(args[0].equals("qualitySangerToIllumina")&&args.length==2){
//				qualitySangerToIllumina(args[1]);
//			}else if(args[0].equals("qualityIlluminaToSanger")&&args.length==2){
//				qualityIlluminaToSanger(args[1]);
			}else if(args[0].equals("reverse")&&args.length==2){
				reverse(general.getBufferedReader(args[1]));
			}else if(args[0].equals("samplePaired")&&args.length==5){
				samplePaired(args[1],args[2],Integer.parseInt(args[3]),args[4]);
			}else if(args[0].equals("shuffle")&&args.length>2){
				ArrayList<String> file1=new ArrayList<String>(), file2= new ArrayList<String>();
				
				for(int i=1;i<1+(args.length-1)/2;++i){
					file1.add(args[i]);
				}
				for(int i=1+(args.length-1)/2;i<args.length;++i){
					file2.add(args[i]);
				}
				
				System.err.println(file1.size()+"\t"+file2.size());
				
				shuffle(file1,file2,false);
			}else if(args[0].equals("shuffleUgly")&&args.length>2){
				ArrayList<String> file1=new ArrayList<String>(), file2= new ArrayList<String>();
				
				for(int i=1;i<1+(args.length-1)/2;++i){
					file1.add(args[i]);
				}
				for(int i=1+(args.length-1)/2;i<args.length;++i){
					file2.add(args[i]);
				}
				
				System.err.println(file1.size()+"\t"+file2.size());
				
				shuffle(file1,file2,true);
			}else if(args[0].equals("splitZ")&&args.length==4){
				split(args[1],Integer.parseInt(args[2]),args[3],true);
			}else if(args[0].equals("splitPaired")&&args.length==5){
				splitPaired(args[1],args[2],args[3],args[4]);
			}else if(args[0].equals("subFiles")&&args.length==3){
				subFiles(args[1],args[2]);
			}else if(args[0].equals("toFasta")&&args.length==2){
				toFasta(new BufferedReader(new FileReader(args[1])));
			}else if(args[0].equals("toFastaPipe")&&args.length==1){
				toFasta(new BufferedReader(new InputStreamReader(System.in)));
//			}else if(args[0].equals("toNewblerHeader")&&args.length==4){
//				toNewblerHeader(args[1],args[2],args[3]);
			}else if(args[0].equals("trimPipe")&&args.length==2){
				trim(new BufferedReader(new InputStreamReader(System.in)),Integer.parseInt(args[1]));
//			}else if(args[0].equals("generateReadPairsFrom454")&&args.length==5){
//				generateReadPairsFrom454(args[1],Integer.parseInt(args[2]),Integer.parseInt(args[3]),Integer.parseInt(args[4]));
			}else{
				System.err.println(printHelp());
				System.exit(616);
			}
		}else{
			System.err.println(printHelp());
			System.exit(616);
		}
	}
	
	private static String printHelp(){
		String help="Usage: illuminaUtils <cmd> <input>\n";
		help+="where <cmd> is:\n";
		help+="information - \n";
		help+="\t<input> = <fastqFile> <burnin> <max nr to sample> <sample  distance> <min word length (inclusive)> <max word length (exclusive)> <outPrefix>\n";
		help+="randomsub - extracts a random subset from a fastq file. Each sequence will be given one chance in N.\n";
		help+="\t<input> = <fastqFile> <N>\n";
		help+="(anti)sub(Z) - extracts a subset from a fastq file\n";
		help+="\t<input> = <fastqFile1> ...<fastqFileN> <subset file>\n";
		help+="subPipe - extracts a subset from a fastq file with a subset file piped through stdin\n";
		help+="\t<input> = <fastqFile1> ...<fastqFileN>\n";
		help+="subToFile - extracts a subset from a fastq file\n";
		help+="\t<input> = <fastqFile1> ...<fastqFileN> <subset file>\n";
		help+="(anti)subReg - extracts a subset from a fastq file\n";
		help+="\t<input> = <fastqFile> <regular expression>\n";
		help+="clean - cleans a fastq file from stretches of the letters in the word. That is, XAT will search for stretches X's, A's or T's with a minimum length of word size in the fasta file and return the reads from the fastq file. Also ditches reads with too many N's\n";
		help+="\t<input> = <fastqFile> <fastaFile> <min length> <word size> <word> <max number of N's>\n";
		help+="cleanVelvetFile - Only prints the complete pairs from a velvet interleaved file\n";
		help+="\t<input> = <fastqFile>\n";
		help+="reverseComplement - reverse complements a fastqfile\n";
		help+="\t<input> = <fastqFile> \n";
		help+="removeShort - drop sequences shorter than n\n";
		help+="\t<input> = <fastqFile> <n>\n";
		help+="removeWindowPaired(Z) - drop sequences between n and m (n<=length<=m\n";
		help+="\t<input> = <fastqFile1> <fastqFile2> <n> <m> <outFile1> <outFile2\n";
		help+="removeShortVelvetPair - drop sequences shorter than n\n";
		help+="\t<input> = <fastqFile> <n> <singletFile> <pairFile>\n";
		help+="convertNumberToPhredQual - converts number qualities to Phred values\n";
		help+="\t<input> = <fastqFile>\n";
		help+="randomVelvetPairedSubset - Selects a random number of reads corresponding to the given fraction (0-1.0) from a paired (interleaved velvet style) fastq file \n";
		help+="\t<input> = <fastqFile> <fraction>\n";
		help+="getPaired - extracts paired sequences from file R1 and R2. These must be ordered in the same way as the raw file and the names must be the same for R1 and R2 and raw\n";
		help+="\t<input> = <fastqFile_raw> <fastqFile R1> <fastqFile R2> <outprefix>\n";
		help+="hpKmerTrim - Trims with respect to the reference. Only leaves sections flanked by kmers present in the reference. The section is padded by two bp in both directions\n";
		help+="\t<input> = <fastaFile ref> <kmer length> <fastqFile R1> <fastqFile R2> <outprefix>\n";
		help+="hpSeqSync - Given n sets of fastq files, they will be syncronized and trimmed to the same length\n";
		help+="\t<input> = <maxShift> <minLength> <maxError> <fastqFile 1 R1> <fastqFile 1 R2> ... <fastqFile n R1> <fastqFile n R2> <outprefix>\n";
		help+="hpSeqRemerge - Replaces the sequences in a sam file (- for stdin) with those from the fastq file\n";
		help+="\t<input> = <sam file> <fastq file>\n";
		help+="indexHammingDist - calculates the Hamming distance between the index of each read (bcl2fastq formatted reads) and the given sequence\n";
		help+="\t<input> = <bcl2fastq fastq file> <sequence>\n";
		help+="indexHammingFilter - Filter reads if the Hamming distance to the given index is larger than the cutoff\n";
		help+="\t<input> = <bcl2fastq file1> (<bcl2fastq file2>) <index> <max distance> <outfile1.gz> (<outfile2.gz>)\n";
		help+="kmerCount - counts the kmers in a fastqfile and prints to stdout. OBS! DOES NOT HANDLE Ns PROPERLY (positions must be merged afterwards) OBS!\n";
		help+="\t<input> = <fastqFile> <kmerSize>\n";
		help+="kmerCount_converts - converts the kmers to their representative kmers and prints \n";
		help+="\t<input> = <kmerCount (kmer\\tcount)>\n";
//		help+="kmerCount_merge - merges several kmerize|sort|uniq -c|awk '{print $1\"\\t\"$2}' files\n";
//		help+="\t<input> = <kmer uniq file 1> <...> <kmer uniq file n> <outFile>\n";
		
		
		
//		help+="kmerCount_histogram - makes a histogram over the kmerCounts\n";
//		help+="\t<input> = <kmer obj file> <bin size>\n";
		
		help+="kmerize - prints all kmers of size n in each sequence\n";
		help+="\t<input> = <fastqFile> <n>\n";
		help+="kmerizePipe - prints all kmers of size n in each sequence\n";
		help+="\t<input> = <n>\n";
		help+="kmerCoverage(Z) - prints kmer-coverage of all sequences in a file with respect to a KmerSet file\n";
		help+="\t<input> = <KmerCountFile (gz)> <cutoff> <fastqFile 1> ... <fastqFile n>\n";
		help+="kmerCoverageToKmerHits - prints kmerHits of all sequences in a file with respect to a Kmer-coverage file\n";
		help+="\t<input> = <KmerCoverageFile> <kmerSize> \n";
		help+="kmerCoverage_stored - prints kmer-coverage of all sequences in a file with respect to a KmerSet file\n";
		help+="\t<input> = <KmerCountFile> <fastqFile 1> ... <fastqFile n>\n";
		help+="kmerCoverage_count - prints total kmer-coverage (count) of all sequences in a file with respect to a KmerSet file\n";
		help+="\t<input> = <KmerCountFile> <cutoff> <fastqFile 1> ... <fastqFile n>\n";
		help+="kmerCoverage_generateCorrection - reads the output from kmerCoverage and generates a file with the longest fully covered subsequence\n";
		help+="\t<input> = <KmerCoverage file> <kmerSize>\n";
		help+="kmerCoverage_extractCovered(Single)(Z) - prints the all pairs with at least one hit against the given kmer distribution. Add Z to accept gzipped files, and single for non-paired files\n";
		help+="\t<input> = <fastqFile(s)1> (<fastqFile(s)2>) <KmerCoverage file> <kmerPos> <outPrefix>\n";
		help+="kmerCoverage_extractPerfect - prints the qname of all sequences with perfect kmer-coverage\n";
		help+="\t<input> = <fastqFile> <KmerCoverage file> <kmerSize>\n";
		help+="kmerCoverage_extractPerfectDirect - prints all sequences with perfect kmer-coverage\n";
		help+="\t<input> = <KmerCount file> <kmerPos> <fastqFile 1> ... <fastqFile n>\n";
		help+="kmerToUse - reverse complements if neccessary\n";
		help+="\t<input> = <kmer-txt file>\n";
		help+="lengths - prints length of sequence\n";
		help+="\t<input> = <fastqFile>\n";
		help+="normalizeKmer - only prints read pairs with kmers that occured less than limit in the previous reads\n";
		help+="\t<input> = <fastq 1> <fastq 2> <kmer size> <limit> <outprefix>\n";
		help+="normalizeKmer2 - only prints read pairs with kmers that occured less than limit in the previous reads\n";
		help+="\t<input> = <fastq 1> <fastq 2> <repetitve kmers> <limit> <outprefix>\n";
		help+="storeKmerSet - stores the kmers (count\\tkmer) that occcurs at least <min cutoff> times in a kmerSet on the disk for future use\n";
		help+="\t<input> = <Kmer count file> <min cutoff> <outfile>\n";
//		help+="reptileCorrection - Applies the corrections calculated by reptile (Assumes that the read index starts at 0)\n";
//		help+="\t<input> = <fastqFile> <reptileFile>\n";
		help+="complement - complements the sequences in the fastqFile\n";
		help+="\t<input> = <fastqFile>\n";
		help+="conversionRateNOMe(Z) - counts the number of C's and G's not affected by NOMe conversion\n";
		help+="\t<input> = <fastqFile>\n";
		help+="generateReadPairsFrom454 - Generates read pairs of length l bp, with the insert size i (can be negative) A read pair tiling starts every 40 bp and once from the end. (The total length of the frament is 2l+i)\n";
		help+="\t<input> = <fastqFile> <l> <i> <cov>\n";
		help+="reverse - reverses a fastq file\n";
		help+="\t<input> = <fastqFile>\n";
		help+="samplePaired - samples a paired fastqFile into n different parts\n";
		help+="\t<input> = <fastqFile1> <fastqFile2> <n> <outPrefix>\n";
		help+="shuffle(Ugly) - generates a interleaved fastq file. If a file name ends with gz it is interpreated as a gzipped archive. Ugly adds /1 /2 to the ends of the names\n";
		help+="\t<input> = <fastqFile1 1> .. <fastqFile1 n> <fastqFile2 1> .. <fastqFile2 n>\n";
		help+="splitZ - splits a gzipped (only thing implemented) into gzipped subfiles containing N reads\n";
		help+="\t<input> = <fastqFile1> <N> <outPrefix>\n";
		help+="splitPaired - splits paired fastq files given a file with the cluster number on each line\n";
		help+="\t<input> = <fastqFile1> <fastqFile2> <clusterFile> <outPrefix>\n";
		help+="subFiles - takes a fastq file and a giTable with filenames in the second column. Writes the ids to the respective files\n";
		help+="\t<input> = <fastqFile> <giTable>\n";
		help+="toFasta(Pipe) - drops the quality values and converts to fasta format\n";
		help+="\t<input> = <fastqFile> \n";
		help+="toNewblerHeader - Shifts the header to include the newbler pairing information\n";
		help+="\t<input> = <fastqFile> <library> <direction (F or R)> \n";
		help+="trimPipe - Trims the reads to the given length\n";
		help+="\t<input> = <length> \n";
		help+="qualityIlluminaToSanger - Shifts the quality score from base 64 to 33... no math\n";
		help+="\t<input> = <fastqFile>\n";
		help+="qualitySangerToIllumina - Shifts the quality score from base 33 to 64... no math\n";
		help+="\t<input> = <fastqFile>\n";
		
		
		
		return help;
	}
	
	private static void indexHammingDist(BufferedReader fqr1, final String index) throws Exception{
		fastqParser fqp1= new fastqParser(fqr1,"");
		FastqSeq fqs1;
		final char[] indexArray= index.toCharArray();
		
		for(;fqp1.hasNext();){
			fqs1=fqp1.next();
			final String header=fqs1.getHeader();
			System.out.println(fqs1.indexHammingDist(indexArray)+sep+index+sep+header.substring(header.lastIndexOf(':')+1));
		}
		
	}
	
	private static void indexHammingFilterSE(BufferedReader fqr1, final String index, final int maxDist, String outFile1) throws Exception{
		fastqParser fqp1= new fastqParser(fqr1,"");
		FastqSeq fqs1;
		final char[] indexArray= index.toCharArray();
		int count=0,
				count2=0;
		BufferedWriter r1Out= new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outFile1))));
		
		for(;fqp1.hasNext();++count){
			if(count%10000==0){
				System.err.print(count+" ### "+count2+"\r");
			}
			fqs1=fqp1.next();
			if(fqs1.indexHammingDist(indexArray)<=maxDist){
				++count2;
				r1Out.write(fqs1.toString()+'\n');
			}
		}
		System.err.println(count+" ### "+count2);
		
		r1Out.close();
	}
	
	private static void indexHammingFilterPE(BufferedReader fqr1, BufferedReader fqr2, final String index, final int maxDist, String outFile1, String outFile2) throws Exception{
		fastqParser fqp1= new fastqParser(fqr1, ""),
				fqp2= new fastqParser(fqr2, "");
		FastqSeq fqs1, fqs2;
		final char[] indexArray= index.toCharArray();
		int count=0,
				count2=0;
		BufferedWriter r1Out= new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outFile1)))),
				r2Out= new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outFile2))));
		
		for(;fqp1.hasNext() && fqp2.hasNext();++count){
			if(count%10000==0){
				System.err.print(count+" ### "+count2+"\r");
			}
			fqs1=fqp1.next();
			fqs2=fqp2.next();
			if(fqs1.indexHammingDist(indexArray)<=maxDist){
				++count2;
				r1Out.write(fqs1.toString()+'\n');
				r2Out.write(fqs2.toString()+'\n');
			}
		}
		System.err.println(count+" ### "+count2);
		
		if(fqp1.hasNext() || fqp2.hasNext()){
			System.err.println("WARNING: One of the files contains more reads than the other");
		}
		
		r1Out.close();
		r2Out.close();
	}
	
	private static void hpSeqRemerge(String samFileName, String fqFileName) throws Exception{
		BufferedReader samFile= general.getBufferedReader(samFileName);
		fastqParser fqp= new fastqParser(general.getBufferedReader(fqFileName),"");
		FastqSeq fq;
		HashMap<String, String> seqStore= new HashMap<String, String>(5500000);
		int count=0,
				count2=0;
		
		System.err.println("Storing sequences...");
		for(;fqp.hasNext();++count){
			if(count%10000==0){
				System.err.print(count+"\r");
			}
			fq=fqp.next();
			seqStore.put(fq.getQname(), fq.getSeq());
		}
		System.err.println(count+"");
		
		count=0;
		System.err.println("Parsing alignment and replacing sequences");
		for(String s=samFile.readLine();s!=null;s=samFile.readLine(),++count){
			if(count%10000==0){
				System.err.print(count+"\r");
			}
			if(s.startsWith("@") || s.length()<1){
//				let header lines and empty lines pass through 
				System.out.println(s);
			}else{
				String[] l=s.split("\t");
				if(l.length<10){
//					lines with too few columns also pass through
					System.out.println(s);
				}else{
					if(seqStore.containsKey(l[0])){
						l[9]=seqStore.get(l[0]);
						System.out.print(l[0]);
						for(int i=1;i<l.length;++i){
							System.out.print("\t"+l[i]);
						}
						System.out.print("\n");
						++count2;
					}else{
//						go back to original line when there isn't a new sequence
						System.out.println(s);
					}
				}
			}
		}
		System.err.println(count+"");
		System.err.println("LOGG reset "+count2+" sequences");
		System.err.println("DONE!");
	}
	//	private static void hpSeqSync(final boolean decision,final int maxShift, final int minLength, final int maxError, ArrayList<String> read1List, ArrayList<String> read2List, String outPrefix) throws Exception{
	private static void hpSeqSync(final boolean decision, final int maxShift, final int minLength, final int maxError, ArrayList<String> read1List, ArrayList<String> read2List, String outPrefix) throws Exception{
		if(minLength<1){
			throw new Error("ERROR: (hpSeqSync) minLength ("+minLength+") must be at least 1");
		}
		fastqParser fqp1, fqp2;
		FastqSeq fq1, fq2;
		BufferedWriter r1Out= new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outPrefix+"_R1_sync.fq.gz")))),
				r2Out= new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outPrefix+"_R2_sync.fq.gz")))),
				stat= new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outPrefix+"_sync.stat.csv.gz")))),
				ru1Out= new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outPrefix+"_R1_unsync.fq.gz")))),
				ru2Out= new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outPrefix+"_R2_unsync.fq.gz"))));	
		int rp=0,
			count=0,
			count2=0,
			errorCount=0,
			shift;
		StringBuffer newSeq1,newSeq2,newQual1,newQual2,errorBuffer;
		String read1, read2, qual1, qual2;
		char c1,c2,p1,p2,q1,q2;
		shift=0;
		HashSet<String> allowedCombinations= new HashSet<String>();
		allowedCombinations.add("AA");
		allowedCombinations.add("CC");
		allowedCombinations.add("GG");
		allowedCombinations.add("TT");
		allowedCombinations.add("GA");
		allowedCombinations.add("TC");
		
		for(;rp<read1List.size() && rp<read2List.size();++rp){
			fqp1= new fastqParser(general.getBufferedReader(read1List.get(rp)),"");
			fqp2= new fastqParser(general.getBufferedReader(read2List.get(rp)),"");
			for(;fqp1.hasNext()&&fqp2.hasNext();++count){
				if(count%10000==0){
					System.err.print(rp+"/"+count+"\r");
				}
				fq1= fqp1.next();
				fq2= fqp2.next();
				read1=fq1.getSeq();
				read2=fq2.getSeq();
				newSeq1= new StringBuffer();
				newSeq2= new StringBuffer();
				newQual1= new StringBuffer();
				newQual2= new StringBuffer();
				errorBuffer= new StringBuffer();
				ArrayList<Object> alignreads = water(read1,read2);
				int s1= (int) alignreads.get(0);
				int s2= (int) alignreads.get(1);
				errorCount= geterrorCount(alignreads.get(2).toString(),alignreads.get(3).toString());
				qual1=fq1.getQuality();
				qual2=fq2.getQuality();
				
				//System.out.println("Read1:");
				//System.out.println(alignreads.get(0));
				if(errorCount<=maxError && alignreads.get(2).toString().length()>=minLength){
					for(int i=0;i<alignreads.get(2).toString().length();++i){
						c1= alignreads.get(2).toString().charAt(i);
						c2= alignreads.get(3).toString().charAt(i);
						if(c1==c2){
							newSeq1.append(c1);
							newSeq2.append(c1);
							newQual1.append(qual1.charAt(s1+i));
							newQual2.append(qual2.charAt(s2+i));
						}else if(c1=='G' && c2=='A'){
							newSeq1.append('G');
							newSeq2.append('R');
							newQual1.append(qual1.charAt(s1+i));
							newQual2.append(qual2.charAt(s2+i));
						}else if(c1=='T' && c2=='C'){
							newSeq1.append('C');
							newSeq2.append('Y');
							newQual1.append(qual1.charAt(s1+i));
							newQual2.append(qual2.charAt(s2+i));
						}else if(c1=='-'){
//							TODO: handle gap in first sequence
							newSeq1.append("");
							s1--;
						}else if(c2=='-'){
//							TODO: handle gap in second sequence
							newSeq2.append("");
							s2--;
						}else{
							newSeq1.append('N');
							newSeq2.append('N');
							errorBuffer.append(c1+""+c2+"\n");
							newQual1.append(qual1.charAt(s1+i));
							newQual2.append(qual2.charAt(s2+i));
						}
					}	
//						check dimers	
					boolean previousN=true;
					p1=newSeq1.charAt(0);
					p2=newSeq2.charAt(0);
					for(int i=1;i<newSeq1.length()-1;++i){
						c1=newSeq1.charAt(i);
						c2=newSeq2.charAt(i);
						if(p1=='C' && c1=='G'){
//							strong symmetry
							if(p2=='C' && c2=='G'){
								newSeq2.setCharAt(i-1, 'F');
							}else if(p2=='C'){
								newSeq2.setCharAt(i-1, 'P');
							}else if(c2=='G'){
								newSeq2.setCharAt(i-1, 'M');
							}else{
								newSeq2.setCharAt(i-1, 'U');
							}
						}else if(p1=='G' && c1=='C'){
//							weak symmetry
							if(p2=='G' && c2=='C'){
								newSeq2.setCharAt(i-1, 'f');
							}else if(p2=='G'){
								newSeq2.setCharAt(i-1, 'p');
							}else if(c2=='C'){
								newSeq2.setCharAt(i-1, 'm');
							}else{
								newSeq2.setCharAt(i-1, 'u');
							}
						}else if(p1=='C' && previousN){
							if(p2=='Y'){
								newSeq2.setCharAt(i-1, 'y');
							}else if(p2=='C'){
								newSeq2.setCharAt(i-1, 'c');
							}
						}else if(p1=='G' && previousN){
							if(p2=='R'){
								newSeq2.setCharAt(i-1, 'r');
							}else if(p2=='G'){
								newSeq2.setCharAt(i-1, 'g');
							}
						}	
						previousN=(p1=='N');
						p1=c1;
						p2=c2;
					}
					stat.write(fq1.getQname()+"\t"+errorCount+"\t"+shift+"\t"+newSeq1.length()+"\n");
					r1Out.write(fq1.getHeader()+"\n"+newSeq1.toString()+"\n+\n"+newQual1.toString()+"\n");
					r2Out.write(fq2.getHeader()+"\n"+newSeq2.toString()+"\n+\n"+newQual2.toString()+"\n");
					++count2;
				}else{
					ru1Out.write(fq1.toString()+"\n");
					ru2Out.write(fq2.toString()+"\n");
				}	
			}
		}
		System.err.print(rp+"/"+count+"\n");
		System.err.println("LOGG: synced "+count2+" sequences");
		stat.close();
		r1Out.close();
		r2Out.close();
		ru1Out.close();
		ru2Out.close();
	}
		
	private static int hpSeqSync_scoreShift(String read1, String read2, int shift,HashSet<String> allowedCombinations){
		//do "alignment"
		//TODO: replace with smith-waterman...
		int errorCount=0;
		for(int i=0,j=shift;i<read1.length() && j<read2.length();++i,++j){
			if(j>-1){
				//String s=read1.charAt(i)+""+read2.charAt(j);
				if(!allowedCombinations.contains(read1.charAt(i)+""+read2.charAt(j))){
					++errorCount;
					//System.err.println("Mismatches: "+ read1.charAt(i)+""+read2.charAt(j));
				}
			}
		}
		return errorCount;
	}
	/**
	 * Score function of Waterman & Smith algorithm
	 * return score of a match, mismatch and gap penalty
	 */
	private static int match_score(char x,char y){
		if(x == y){
			return match_award;
		}
		else if(x == 'G' && y == 'A'){
			return match_award;
		}
		else if(x == 'T' && y == 'C'){
			return match_award;
		}
		else if(x == '-' || y == '-'){
			return gap_penalty;
		}
		else {
			return mismatch_penalty;
		}
	}
	
	private static int geterrorCount(String read1, String read2){
		HashSet<String> allowedCombinations= new HashSet<String>();
		allowedCombinations.add("AA");
		allowedCombinations.add("CC");
		allowedCombinations.add("GG");
		allowedCombinations.add("TT");
		allowedCombinations.add("GA");
		allowedCombinations.add("TC");
		int errorCount=0;
		for(int i=0,j=0;i<read1.length() && j<read2.length();++i,++j){
			if(!allowedCombinations.contains(read1.charAt(i)+""+read2.charAt(j))){
				++errorCount;
			}
		}
		return errorCount;
	}
	private static ArrayList<Object> water(String seq1, String seq2){
		/**
		* Waterman & Smith algorithm:
		* Create Dynamic Programming table and trace back with pointers
		*/
		String result1 ="";
		String result2 ="";
		int length1 =seq1.length();
		int length2 =seq2.length();
		score = new int[length1 + 1][length2 + 1];
		pointer = new int[length1 + 1][length2 + 1];
		int max_score = 0; 
		int max_i = 0;
		int max_j = 0;
		for(int i = 1; i < (length1+1); ++i){
			for(int j =1; j < (length2+1); ++j){
				int score_diagonal = score[i-1][j-1] + match_score(seq1.charAt(i-1),seq2.charAt(j-1));
				int score_up = score[i][j-1] + gap_penalty;
				int score_left = score[i-1][j] + gap_penalty;
				score[i][j] =  Math.max(Math.max(0, score_left),Math.max(score_up, score_diagonal));
				if(score[i][j]==0){
					pointer[i][j] = 0; //0 means end of the path
				}
				if(score[i][j] == score_left){
					pointer[i][j] = 1; //1 means trace left
				}
				if(score[i][j] == score_up){
					pointer[i][j] = 2; //2 means trace up
				}
				if(score[i][j] == score_diagonal){
					pointer[i][j] = 3; //3 means trace diagonal
				}
				if(score[i][j] >= max_score){
					max_i = i;
					max_j = j;
					max_score = score[i][j];
				}
			}
		}
		String align1 = "";
		String align2 = "";
		int i = max_i;
		int j = max_j;
		while(pointer[i][j]!=0){
			if (pointer[i][j] == 3){
				align1 += seq1.charAt(i-1);
				align2 += seq2.charAt(j-1);
				i-=1;
				j-=1;
			}else if(pointer[i][j] == 2){
				align1 += "-";
				align2 += seq2.charAt(j-1);
				j-=1;
			}else if(pointer[i][j] == 1){
				align1 += seq1.charAt(i-1);
				align2 += "-";
				i-=1;
			}
		}
		ArrayList<Object> result = new ArrayList<Object>();
		result.add(i);
		result.add(j);
		result1 = new StringBuffer(align1).reverse().toString();
		result.add(result1.toString());
		result2 = new StringBuffer(align2).reverse().toString();
		result.add(result2.toString());
		return result;
	}
	
	
	
	private static void hpKmerTrim(String faName, int kmerSize,ArrayList<String> read1List, ArrayList<String> read2List, String outPrefix, boolean isHP, int padding) throws Exception{
		final int minLength=40;
		fastaParser fp= new fastaParser(general.getBufferedReader(faName));
		FastaSeq fs;
		fastqParser fqp1, fqp2;
		FastqSeq fq1, fq2;
		KmerSet_binary_utils ku= new KmerSet_binary_utils(kmerSize);
//		might need to be made larger
		HashSet<BitSet> refKmers= new HashSet<BitSet>();
		ArrayList<BitSet> read1Kmers= new ArrayList<BitSet>(),
				read2Kmers= new ArrayList<BitSet>();
		int count=0,rp=0,start,end;
		BufferedWriter r1Out, r2Out;
		if(isHP){
			r1Out= new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outPrefix+"_R1.fq.gz"))));
			r2Out= new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outPrefix+"_R2.fq.gz"))));
		}else{
			r1Out= new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outPrefix+"_merged.fq.gz"))));
			r2Out= new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outPrefix+"_emptyPlaceholderFile.gz"))));
		}
		System.err.println("Processed reference sequences:");
		for(;fp.hasNext();++count){
			if(count%1000==0){
				System.err.print(count+"\r");
			}
			fs= fp.next();
//			add the bisulfite converted forward strand
			ku.kmerizeNoN(fs.getSeq().replace('C', 'T'), refKmers);
//			add the bisulfite converted reverse strand
			ku.kmerizeNoN(sequenceUtils.reverseComplement(fs.getSeq()).replace('C', 'T'), refKmers);
		}
		System.err.println(count);
		System.err.println("These contained "+ refKmers.size() + " unique bisulfite kmers");
		count=0;
		System.err.println("Processed files/read pairs:");
		for(;rp<read1List.size() && rp<read2List.size();++rp){
			fqp1= new fastqParser(general.getBufferedReader(read1List.get(rp)),"");
			fqp2= new fastqParser(general.getBufferedReader(read2List.get(rp)),"");
			for(;fqp1.hasNext()&&fqp2.hasNext();++count){
				if(count%10000==0){
					System.err.print(rp+"/"+count+"\r");
				}
				fq1= fqp1.next();
				fq2= fqp2.next();
				read1Kmers.clear();
				read2Kmers.clear();
				ku.kmerize(fq1.getSeq().replace('C', 'T'), read1Kmers);
				ku.kmerize(sequenceUtils.reverseComplement(fq2.getSeq()).replace('C', 'T'), read2Kmers);
				if(isHP){
//					assume the reads are in sync...
//					and trim them together
					
//					find first position
					start=0;
					for(;start<read1Kmers.size() && start<read2Kmers.size();++start){
						if(refKmers.contains(read1Kmers.get(start)) || refKmers.contains(read2Kmers.get(read2Kmers.size()-1-start))){
//							found start position
//							find end point
							end=start;
//							this would better be done by searching from the end of the read, but then it would be better to align them first... different lengths etc.
							for(int j=start+1;j<read1Kmers.size() && j<read2Kmers.size();++j){
								if(refKmers.contains(read1Kmers.get(j)) || refKmers.contains(read2Kmers.get(read2Kmers.size()-1-j))){
									end=j;
								}
							}
							start=start<padding?0:start-padding;
							end= end+kmerSize+padding;
							if(end-start>minLength){
								r1Out.write(fq1.getSubFastqSeq(start, end>fq1.length()?fq1.length():end).toString()+"\n");
								r2Out.write(fq2.getSubFastqSeq(start, end>fq2.length()?fq2.length():end).toString()+"\n");
							}
							break;
						}
					}
				}else{
//					trim reads individually
//					the files will be completely out of sync... should be aligned separately
					
//					assume the reads are in sync...
//					and trim them together
					
//					process first read
//					find first position
					start=0;
					for(;start<read1Kmers.size();++start){
						if(refKmers.contains(read1Kmers.get(start))){
//							found start position
//							find end point
							end=start;
//							this would better be done by searching from the end of the read, but then it would be better to align them first... different lengths etc.
							for(int j=start+1;j<read1Kmers.size();++j){
								if(refKmers.contains(read1Kmers.get(j))){
									end=j;
								}
							}
							start=start<padding?0:start-padding;
							end= end+kmerSize+padding;
							if(end-start>minLength)
								r1Out.write(fq1.getSubFastqSeq(start, end>fq1.length()?fq1.length():end).toString("@"+fq1.getQname()+"_1")+"\n");
							break;
						}
					}
//					process second read
//					find first position
					start=0;
					for(;start<read2Kmers.size();++start){
						if(refKmers.contains(read2Kmers.get(start))){
//							found start position
//							find end point
							end=start;
//							this would better be done by searching from the end of the read, but then it would be better to align them first... different lengths etc.
							for(int j=start+1;j<read2Kmers.size();++j){
								if(refKmers.contains(read2Kmers.get(j))){
									end=j;
								}
							}
							start=start<padding?0:start-padding;
							end= end+kmerSize+padding;
							if(end-start>minLength)
								r1Out.write(fq2.getSubFastqSeq(start, end>fq2.length()?fq2.length():end).reverseComplement().toString("@"+fq2.getQname()+"_2")+"\n");
							break;
						}
					}

					
				}
			}
		}
		System.err.println(rp+"/"+count);
		r1Out.close();
		r2Out.close();
	}
	
	
	private static void getPaired(BufferedReader rawReader, BufferedReader r1Reader, BufferedReader r2Reader, String outprefix) throws Exception{
		
		//TODO: add support for keeping singlet reads
		fastqParser rawParser= new fastqParser(rawReader, ""),
				r1Parser= new fastqParser(r1Reader, ""),
				r2Parser= new fastqParser(r2Reader, "");
		BufferedWriter r1Out= new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outprefix+".R1.fq.gz")))),
				r2Out= new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outprefix+".R2.fq.gz"))));
		FastqSeq raw= rawParser.next(),
				r1= r1Parser.next(),
				r2= r2Parser.next();
		
		while(raw!=null && r1!=null && r2!=null){
			if(r1.getQname().equals(r2.getQname())){
				r1Out.write(r1.toString()+'\n');
				r1= r1Parser.hasNext()?r1Parser.next():null;
				r2Out.write(r2.toString()+'\n');
				r2= r2Parser.hasNext()?r2Parser.next():null;
			}else{
				//find the side that is lagging behind
				if(r1.getQname().equals(raw.getQname())){
					r1= r1Parser.hasNext()?r1Parser.next():null;
				}else if(r2.getQname().equals(raw.getQname())){
					r2= r2Parser.hasNext()?r2Parser.next():null;
				}else{
					raw= rawParser.hasNext()?rawParser.next():null;
				}
			}
		}
		
		r1Out.close();
		r2Out.close();
	}
	
	private static void complement(BufferedReader fastq)throws Exception {
		fastqParser fp= new fastqParser(fastq,"");
		FastqSeq fs;
		for(;fp.hasNext();){
			fs=fp.next();
			fs.setSeq(sequenceUtils.Complement(fs.getSeq()));
			System.out.println(fs.toString());
		}
	}
	
	private static void lengths(BufferedReader fastq)throws Exception{
		fastqParser fp= new fastqParser(fastq, "");
		FastqSeq fs;
		for(;fp.hasNext();){
			fs=fp.next();
			System.out.println(fs.getQname()+sep+fs.getSeq().length());
		}
	}
	
	private static void reverse(BufferedReader fastq)throws Exception{
		fastqParser fp= new fastqParser(fastq, "");
		FastqSeq fs;
		int count=0;
		for(;fp.hasNext();++count){
			if(count%10000==0){
				System.err.print(count+"\r");
			}
			fs=fp.next();
			fs.reverseThis();
			System.out.println(fs.toString());
		}
		System.err.print(count+"\n");
	}
	
	public static void removeWindowPaired(String file1,String file2,int n,int m,String outFile1, String outFile2,boolean gzipped)throws Exception{
		fastqParser fp1= new fastqParser(file1, gzipped), fp2= new fastqParser(file2, gzipped);
		FastqSeq fqs1,fqs2;
		BufferedWriter out1=new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outFile1)))),
				out2=new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outFile2))));
		int count=0,kept=0;
		for(;fp1.hasNext()&&fp2.hasNext();++count){
			fqs1=fp1.next();
			fqs2=fp2.next();
			if((fqs1.length()<n || fqs1.length()>m) && (fqs2.length()<n || fqs2.length()>m)){
				++kept;
				out1.write(fqs1.toString()+'\n');
				out2.write(fqs2.toString()+'\n');
			}
		}
		System.err.println(count+sep+kept);
		out1.close();
		out2.close();
	}
	
	public static void conversionRateNOMe(String inFile, boolean gzipped)throws Exception{
		fastqParser fp= new fastqParser(inFile,gzipped);
		FastqSeq fqs;
//		String Cexp= "[ATC]C[ATC]", Gexp="[ATG]G[ATG]";
		long Ccount=0L, Gcount=0L;
		for(;fp.hasNext();){
			fqs=fp.next();
			if(fqs.length()>3){
				final char[] seq= fqs.getSeq().toCharArray();
				for( int i=1;i<seq.length-1;++i){
					switch (seq[i]) {
					case 'C':
						if(seq[i-1]!='G' && seq[i+1]!='G'){
							++Ccount;
						}
						break;
					case 'G':
						if(seq[i-1]!='C' && seq[i+1]!='C'){
							++Gcount;
						}
						break;
					default:
						break;
					}
				}
			}
		}
		System.out.println(inFile+sep+Ccount+sep+Gcount);
	}
	

	
	public static void split(String inFile,int N, String outPrefix,boolean gzipped)throws Exception{
		fastqParser fp= new fastqParser(inFile,gzipped);
//		 bufferedWriter = new BufferedWriter(
//                 new OutputStreamWriter(
//                     new GZIPOutputStream(new FileOutputStream(outFilename))
//                 ));
		int fileCounter=1;
		String fileSuffix="0001";
		BufferedWriter out=new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outPrefix+"_"+fileSuffix+".fq.gz"))));
		for(int readCounter=0;fp.hasNext();++readCounter){
			if(readCounter==N){
				//Start new file
				out.close();
				++fileCounter;
				fileSuffix=""+fileCounter;
				while(fileSuffix.length()<4){
					fileSuffix="0"+fileSuffix;
				}
				out=new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outPrefix+"_"+fileSuffix+".fq.gz"))));
				readCounter=0;
			}
			out.write(fp.next().toString()+"\n");
		}
		out.close();
	}
	
	public static void subFiles(String fastqFile,String giTableFile)throws Exception{
		HashMap<String, BufferedWriter> outputMap= new HashMap<String, BufferedWriter>();
		HashMap<String, BufferedWriter> idMap= new HashMap<String, BufferedWriter>(6000000);
		BufferedReader giTable= new BufferedReader(new FileReader(giTableFile));
		fastqParser fqp= new fastqParser(fastqFile);
		FastqSeq fqs;
		String[] l;
		
		for(String s=giTable.readLine();s!=null;s=giTable.readLine()){
			l=s.split("\t");
			if(!outputMap.containsKey(l[1])){
				outputMap.put(l[1], new BufferedWriter(new FileWriter(l[1])));
			}
			idMap.put(l[0], outputMap.get(l[1]));
		}
		giTable.close();
		
		for(;fqp.hasNext();){
			fqs=fqp.next();
			if(idMap.containsKey(fqs.getQname())){
				idMap.get(fqs.getQname()).write(fqs.toString()+"\n");
			}
		}
		for(BufferedWriter toClose : outputMap.values()){
			toClose.close();
		}
	}
	
//	public static void normalizeKmer2(String file1,String file2,String repKmerfile,int limit,String outPrefix)throws Exception{
//		fastqParser fqp1=new fastqParser(new BufferedReader(new FileReader(file1)), "");
//		fastqParser fqp2=new fastqParser(new BufferedReader(new FileReader(file2)), "");
//		HashMap<BitSet, Integer> stillOK= new HashMap<BitSet, Integer>(1000000000);
//		HashSet<BitSet> rep= new HashSet<BitSet>(1000000000);
//		BufferedReader kmers=new BufferedReader(new FileReader(repKmerfile));
//		String kmer=kmers.readLine();
//		ArrayList<BitSet> curKmers=new ArrayList<BitSet>();
//		KmerSet_binary_utils ku= new KmerSet_binary_utils(kmer.length());
//		stillOK.put(ku.stringToBitSet(kmer), limit);
//		boolean print;
//		int count,i=0,printed=0;
//		FastqSeq fs1,fs2;
//		BufferedWriter out1=new BufferedWriter(new FileWriter(outPrefix+"_1.fq")),
//				out2=new BufferedWriter(new FileWriter(outPrefix+"_2.fq"));
//		System.err.println("reading repetitive kmers ...");
//		
//		for(kmer=kmers.readLine();kmer!=null;kmer=kmers.readLine(),++i){
//			if(i%1000==0){
//				System.err.print("     "+i+"\r");
//			}
//			stillOK.put(ku.stringToBitSet(kmer), limit);
//		}
//		kmers.close();
//		System.err.println("     "+i);
//		System.err.println("parsing fastq file");
//		i=0;
//		for(;fqp1.hasNext()&&fqp2.hasNext();++i){
//			if(i%1000==0){
//				System.err.print("     "+i+"   "+printed+"\r");
//			}
//			fs1=fqp1.next();
//			fs2=fqp2.next();
//			print=false;
//			curKmers.clear();
//			ku.kmerizeNoN(fs1.getSeq(), curKmers);
//			ku.kmerizeNoN(fs2.getSeq(), curKmers);
//			for(BitSet k : curKmers){
//				
//				if (!rep.contains(k)) {
//					if (stillOK.containsKey(k)) {
//						count = stillOK.get(k);
//						if (count == 0) {
//							rep.add(k);
//						} else {
//							print = true;
//							stillOK.put(k, count - 1);
//						}
//					} else {
//						print = true;
//					}
//				}
//			}
//			if(print){
//				++printed;
//				out1.write(fs1.toString());
//				out2.write(fs2.toString());
//			}
//			
//		}
//		System.err.println("     "+i+"   "+printed);
//		out1.close();
//		out2.close();
//		
//	}
//	
//	public static void normalizeKmer(String file1, String file2, int kmerSize, int limit,String outPrefix)throws Exception{
//		fastqParser fqp1=new fastqParser(new BufferedReader(new FileReader(file1)), "");
//		fastqParser fqp2=new fastqParser(new BufferedReader(new FileReader(file2)), "");
//		HashSet<BitSet> rep= new HashSet<BitSet>(1000000000);
//		ArrayList<BitSet> curKmers=new ArrayList<BitSet>();
//		HashMap<BitSet, Integer> ok= new HashMap<BitSet, Integer>(1000000000);
//		KmerSet_binary_utils ku= new KmerSet_binary_utils(kmerSize);
//		boolean print;
//		int count,i=0,printed=0;
//		FastqSeq fs1,fs2;
//		BufferedWriter out1=new BufferedWriter(new FileWriter(outPrefix+"_1.fq")),
//				out2=new BufferedWriter(new FileWriter(outPrefix+"_2.fq"));
//		for(;fqp1.hasNext()&&fqp2.hasNext();++i){
//			if(i%1000==0){
//				System.err.print("     "+i+"   "+printed+"   "+rep.size()+"   "+ok.size()+"\r");
//			}
//			fs1=fqp1.next();
//			fs2=fqp2.next();
//			print=false;
//			curKmers.clear();
//			ku.kmerizeNoN(fs1.getSeq(), curKmers);
//			ku.kmerizeNoN(fs2.getSeq(), curKmers);
//			for(BitSet kmer : curKmers){
//				if(rep.contains(kmer)){
//					
//				}else{
//					print=true;
//					if(ok.containsKey(kmer)){
//						count=ok.get(kmer);
//						if(count==limit){
//							ok.remove(kmer);
//							rep.add(kmer);
//						}else{
//							ok.put(kmer, count+1);
//						}
//					}else{
//						ok.put(kmer, 2);
//					}
//				}
//			}
//			if(print){
//				++printed;
//				out1.write(fs1.toString());
//				out2.write(fs2.toString());
//			}
//		}
//		System.err.println("     "+i+"   "+printed+"   "+rep.size()+"   "+ok.size());
//		out1.close();
//		out2.close();
//	}
	
	public static void splitPaired(String file1, String file2, String clusterFile, String outPrefix)throws Exception{
		fastqParser fqp1=new fastqParser(new BufferedReader(new FileReader(file1)), "");
		fastqParser fqp2=new fastqParser(new BufferedReader(new FileReader(file2)), "");
		BufferedReader clusters= new BufferedReader(new FileReader(clusterFile));
		BufferedWriter out1,out2;
//		FastqSeq fs1,fs2;
//		String last="";
		for(String s=clusters.readLine();s!=null;s=clusters.readLine()){
			out1=new BufferedWriter(new FileWriter(outPrefix+"_"+s+"_R1.fq"));
			out2=new BufferedWriter(new FileWriter(outPrefix+"_"+s+"_R2.fq"));

			out1.append(fqp1.next().toString());
			out2.append(fqp2.next().toString());
			out1.close();
			out2.close();
			
		}
		clusters.close();
	}
	
	public static void shuffle(ArrayList<String> file1, ArrayList<String> file2,final boolean ugly)throws Exception {
		fastqParser fqp1,fqp2;
		FastqSeq fqs1,fqs2;
		
		for(int i=0;i<file1.size();++i){
			System.err.println(file1.get(i)+"\t"+file2.get(i));
			if(file1.get(i).endsWith("gz")){
				fqp1= new fastqParser(new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file1.get(i))))), "");
			}else{
				fqp1= new fastqParser(new BufferedReader(new FileReader(file1.get(i))), "");
			}
			if(file2.get(i).endsWith("gz")){
				fqp2= new fastqParser(new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file1.get(i))))), "");
				}else{
				fqp2= new fastqParser(new BufferedReader(new FileReader(file1.get(i))), "");
			}
			for(;fqp1.hasNext()&&fqp2.hasNext();){
				if(ugly){
					fqs1=fqp1.next();
					fqs2=fqp2.next();
					fqs1.addToQname("/1");
					fqs2.addToQname("/2");
					System.out.println(fqs1);
					System.out.println(fqs2);
				}else{
					System.out.println(fqp1.next());
					System.out.println(fqp2.next());
				}
			}
		}
	}
	
	public static void trim(BufferedReader in, int length)throws Exception{
		fastqParser fqp= new fastqParser(in, "");
		FastqSeq fs;
		for(;fqp.hasNext();){
			fs=fqp.next();
			if(fs.length()>length)
				fs.trimThis(length);
			System.out.println(fs.toString());
		}
	}
	
//	public static void kmerToUse(String kmerFile) throws Exception{
//		BufferedReader in= new BufferedReader(new FileReader(kmerFile));
//		for(String s=in.readLine();s!=null;s=in.readLine()){
//			System.out.println(kmerFunctions.kmerToUse(s));
//		}
//		in.close();
//	}
	
	public static void samplePaired(String fastqFile1, String fastqFile2,int n,String outprefix)throws Exception{
		fastqParser fqp1= new fastqParser(new BufferedReader(new FileReader(fastqFile1)),"");
		fastqParser fqp2= new fastqParser(new BufferedReader(new FileReader(fastqFile2)),"");
		BufferedWriter[] out1= new BufferedWriter[n],out2=new BufferedWriter[n];
		for(int i=0;i<n;++i){
			String nr="0"+i;
			while(nr.length()<3){
				nr="0"+nr;
			}
			out1[i]=new BufferedWriter(new FileWriter(outprefix+"_"+nr+"_1.fq"));
			out2[i]=new BufferedWriter(new FileWriter(outprefix+"_"+nr+"_2.fq"));
		}
		Random rand= new Random();
		int cur;
		for(;fqp1.hasNext()&&fqp2.hasNext();){
			cur=rand.nextInt(n);
			out1[cur].write(fqp1.next()+"\n");
			out2[cur].write(fqp2.next()+"\n");
		}
		for(int i=0;i<n;++i){
			out1[i].close();
			out2[i].close();
		}
	}
	
	public static void toFasta(BufferedReader in)throws Exception{
		fastqParser fqp= new fastqParser(in,"");
		for(;fqp.hasNext();){
			System.out.println(fqp.next().toFastaSeq());
		}
		
	}
	
//	public static void toNewblerHeader(String fastqFile,String library, String direction)throws Exception{
//		fastqParser fqp= new fastqParser(new BufferedReader(new FileReader(fastqFile)),"");
//		FastqSeq fqs;
//		for(;fqp.hasNext();){
//			fqs=fqp.next();
//			fqs.toNewblerHeader(library, direction);
//			System.out.println(fqs);
//		}
//	}
	
//	public static void kmerCount_convert(String kmerFile)throws Exception{
//		BufferedReader in= new BufferedReader(new FileReader(kmerFile));
//		String[] l;
//		for(String s=in.readLine();s!=null;s=in.readLine()){
//			l=s.split("\t");
//			System.out.println(l[1]+"\t"+kmerFunctions.kmerToUse(l[0]));
//		}
//		in.close();
//	}
	
//	public static void qualitySangerToIllumina(String fastqFile)throws Exception{
//		fastqParser fqp= new fastqParser(new BufferedReader(new FileReader(fastqFile)),"");
//		for(;fqp.hasNext();){
//			System.out.println(fqp.next().changeQualBase(33, 64,104));
//		}
//	}
//	
//	public static void qualityIlluminaToSanger(String fastqFile)throws Exception{
//		fastqParser fqp= new fastqParser(new BufferedReader(new FileReader(fastqFile)),"");
//		for(;fqp.hasNext();){
//			System.out.println(fqp.next().changeQualBase(64, 33,73));
//		}
//	}
	
//	public static void cleanVelvetFile(String fastqFile)throws Exception{
//		fastqParser fqp= new fastqParser(new BufferedReader(new FileReader(fastqFile)),"");
//		FastqSeq previous=new FastqSeq("@QWEQRTZZZZF", "", ""),current;
//		for(int i=0,printed=0;fqp.hasNext();i++){
//			if(i%10000==0||printed%10000==0){
//				System.err.print("    "+i+"  "+printed+"\r");
//			}
//			current=fqp.next();
////			System.err.print(current.getQname()+"\r");
//			if(current.getQname().startsWith(previous.getQname().substring(0,previous.getQname().length()-2))){
//				System.out.println(previous+"\n"+current);
//				printed++;
//			}else{
//				previous=new FastqSeq(current);
//			}
//		}
//	}
	
//	public static void generateReadPairsFrom454(String fastqFile,int readLength,int insertSize,int coverage)throws Exception{
////		System.out.println(getPrefix(fastqFile));
//		fastqParser fqp= new fastqParser(new BufferedReader(new FileReader(fastqFile)),"");
//		FastqSeq fqs;
//		int fragmentSize=2*readLength+insertSize;
//		for(;fqp.hasNext();){
//			fqs=fqp.next().changeQualBase(33, 64,104);
//			if(fqs.length()>=fragmentSize){
//				int count=0;
//				for(int i=0;i*40<fragmentSize;i++){
//					for(int start=i*40;start+fragmentSize<=fqs.length();start+=fragmentSize,count++){
//						System.out.println(fqs.toString(start, start+readLength, fqs.getHeader().substring(1)+"_"+count+"\\1"));
//						System.out.println(fqs.getSubFastqSeq(start+fragmentSize-readLength, start+fragmentSize).reverseComplement().toString(fqs.getHeader().substring(1)+"_"+count+"\\2"));
//					}
////					if(i%2==0){
////						//even: start from the beginning
////						for(int start=i*fragmentSize/coverage/2;start+fragmentSize<=fqs.length();start+=fragmentSize,count++){
////							System.out.println(fqs.toString(start, start+readLength, fqs.getHeader().substring(1)+"_"+count+"\\1"));
////							System.out.println(fqs.getSubFastqSeq(start+fragmentSize-readLength, start+fragmentSize).reverseComplement().toString(fqs.getHeader().substring(1)+"_"+count+"\\2"));
////						}
////					}else{
////						//odd: start from the end
////						for(int start=fqs.length()-(i-1)*fragmentSize/coverage/2-fragmentSize;start>=0;start-=fragmentSize,count++){
////							System.out.println(fqs.toString(start, start+readLength, fqs.getHeader().substring(1)+"_"+count+"\\1"));
////							System.out.println(fqs.getSubFastqSeq(start+fragmentSize-readLength, start+fragmentSize).reverseComplement().toString(fqs.getHeader().substring(1)+"_"+count+"\\2"));
////						}
////					}
//				}
//				//Do one run from the end
//				for(int start=fqs.length()-fragmentSize;start>=0;start-=fragmentSize,count++){
//					System.out.println(fqs.toString(start, start+readLength, fqs.getHeader().substring(1)+"_"+count+"\\1"));
//					System.out.println(fqs.getSubFastqSeq(start+fragmentSize-readLength, start+fragmentSize).reverseComplement().toString(fqs.getHeader().substring(1)+"_"+count+"\\2"));
//				}
//			}
//		}
//	}
	
//	private static void generateReadPairsFrom454_print(FastqSeq fqs,int start, int readLength,int fragmentSize,int count)throws Exception{
//		System.out.println(fqs.changeQualBase(33, 64).toString(start, start+readLength, fqs.getHeader().substring(1)+"_"+count+"\\1"));
//		System.out.println(fqs.changeQualBase(33, 64).getSubFastqSeq(start+fragmentSize-readLength, start+fragmentSize).reverseComplement().toString(fqs.getHeader().substring(1)+"_"+count+"\\2"));
//	}


//	public static void kmerCoverageToKmerHits(String kmerCoverageFile,int kmerSize)throws Exception{
//		RocheQualParser rqp= new RocheQualParser(new BufferedReader(new FileReader(kmerCoverageFile)));
//		RocheQualSeq rqs;
//		int[] cov;
//		for(int i=0;rqp.hasNext();++i){
//			if(i%10000==0){
//				System.err.print("    "+i+"        \r");
//			}
//			rqs=rqp.next();
//			cov=rqs.getIntSeq();
//			System.out.print(rqs.getHeader()+"\n"+cov[0]);
//			for(int j=1;j<=cov.length-kmerSize;j++){
//				if(cov[j]>cov[j-1]||cov[j]==kmerSize){
//					System.out.print(sep+"1");
//				}else{
//					System.out.print(sep+"0");
//				}
//			}
//			System.out.println();
//		}
//	}
//	
//		
//	public static void kmerCoverage_extractPerfectDirect(KmerSet_binary kmers,ArrayList<String> fastqFiles)throws Exception{
//		fastqParser fqp;
//		FastqSeq fqs;
//		BufferedWriter out;
//		int count;
//		System.err.println("parsing fastq files");
//		for (String fileName : fastqFiles) {
//			count=0;
//			out= new BufferedWriter(new FileWriter(fileName+".perfect.fq"));
//			fqp=new fastqParser(new BufferedReader(new FileReader(fileName)), "");
//			for(;fqp.hasNext();++count){
//				if(count%10000==0){
//					System.err.print(fileName+":  "+count+"\r");
//				}
//				fqs=fqp.next();
//				if(kmers.isPerfectND(fqs.getSeq())){
//					fqs.write(out);
//				}
//			}
//			System.err.println(fileName+":  "+count);
//			out.close();
//		}
//	}
//
//	public static void kmerCoverage_extractPerfect(String fastqFile, String coverageFile,int kmerSize)throws Exception{
//		RocheQualParser rqp= new RocheQualParser(new BufferedReader(new FileReader(coverageFile)));
//		RocheQualSeq rqs;
//		fastqParser fqp= new fastqParser(new BufferedReader(new FileReader(fastqFile)),"");
//		FastqSeq fqs;
//		String[] quals;
////		boolean perfect;
//		int leftCov,maxCov;
//		for(int i=0;rqp.hasNext()&&fqp.hasNext();i++){
//			if(i%10000==0){
//				System.err.print("\t"+i+"\r");
//			}
//			rqs=rqp.next();
//			
//			quals=rqs.getSeq().trim().split("\\s+");
////			perfect=true;
//			fqs=fqp.next();
//			maxCov=fqs.length()+1-kmerSize<kmerSize?fqs.length()+1-kmerSize:kmerSize;
//			if(Integer.parseInt(quals[0])==1&&Integer.parseInt(quals[quals.length-1])==1&&Integer.parseInt(quals[quals.length-kmerSize])==maxCov){
//				leftCov=maxCov*(fqs.length()+1-maxCov);
//				for(int j=0;j<quals.length;j++){
//					leftCov-=Integer.parseInt(quals[j]);
//				}
//				if(leftCov==0){
//					System.out.println(fqs);
//				}
////				for(int j=kmerSize;j<quals.length-kmerSize;j+=kmerSize-1){
////					if(Integer.parseInt(quals[j])!=kmerSize){
////						perfect=false;
////						break;
////					}
////				}
////				if(perfect){
////					System.out.println(fqs);
////				}
//			}else{
//				if(!rqs.getQname().equals(fqs.getQname())){
//					System.err.println("Out of sync, q: "+rqs.getQname()+", fq: "+fqs.getQname());
//				}
//			}
//
//		}
//	}
//
//	
//	
//	public static void kmerCoverage_generateCorrection(String coverageFile,int kmerSize)throws Exception{
//		RocheQualParser rqp= new RocheQualParser(new BufferedReader(new FileReader(coverageFile)));
//		RocheQualSeq rqs;
//		int[] coords,quals;
//		for(;rqp.hasNext();){
//			rqs=rqp.next();
//			quals=rqs.getIntSeq();
//			coords=kmerCoverage_generateCorrection_longestMaxSub(quals);
//			System.out.println(rqs.getQname()+sep+quals[coords[0]]+sep+coords[0]+sep+coords[1]+sep+kmerCoverage_count(quals,kmerSize));
//		}
//	}
//	
//	public static int kmerCoverage_count(int[] quals,int kmerSize){
//		int count=quals[0];
//		for(int i=1;i<quals.length;i++){
//			if(quals[i]>quals[i-1]||quals[i]==kmerSize){
//				count++;
//			}
//		}
//		return count;
//	}
//	
//	
//	public static int[] kmerCoverage_generateCorrection_longestMaxSub(int[] quals){
//		int[] coords= new int[]{0,0};
//		int[] tmpcoords=new int[]{-1,-1};
//		for(int i=0;i<quals.length;i++){
//			if(quals[i]>quals[coords[0]]){
//				//found a new max
//				coords[0]=i;
//				coords[1]=i+1;
//				tmpcoords=new int[]{-1,-1};
//			}else if(quals[i]==quals[coords[0]]){
//				if(i==coords[1]){
//					coords[1]++;
//				}else if(tmpcoords[0]!=-1){
//					//start a second stretch
//					tmpcoords[0]=i;
//					tmpcoords[1]=i+1;
//				}else if(tmpcoords[1]==i){
//					tmpcoords[1]++;
//				}else{
//					//check if a longer stretch is found and start a new
//					if(tmpcoords[1]-tmpcoords[0]>coords[1]-coords[0]){
//						coords[0]=tmpcoords[0];
//						coords[1]=tmpcoords[1];
//					}
//					tmpcoords[0]=i;
//					tmpcoords[1]=i+1;
//				}
//			}
//		}
//		if(tmpcoords[1]-tmpcoords[0]>coords[1]-coords[0]){
//			coords[0]=tmpcoords[0];
//			coords[1]=tmpcoords[1];
//		}
//		return coords;
//	}
//	
////	public static void reptileCorrection(String fastqFile,String reptileFile)throws Exception{
////		fastqParser fqp= new fastqParser(new BufferedReader(new FileReader(fastqFile)),"");
////		FastqSeq fqs;
////		int nextToBeRead=0,changed=0,base=0;
////		ReptileLine rl;
////		BufferedReader reptileIn= new BufferedReader(new FileReader(reptileFile));
////		for(String s=reptileIn.readLine();s!=null;s=reptileIn.readLine()){
////			rl=new ReptileLine(s);
////			for(;fqp.hasNext()&&nextToBeRead<rl.getReadId()-base;){
////				System.out.println(fqp.next());
////				nextToBeRead++;
////			}
////			if(fqp.hasNext()){
////				fqs=fqp.next();
////				nextToBeRead++;
////				fqs.correct(rl.getCorrections());
////				changed++;
////				System.out.println(fqs);
////			}else{
////				reptileIn.close();
////				throw new Exception("No reads to cover for error with readID: "+rl.getReadId());
////			}
////		}
////		reptileIn.close();
////		for(;fqp.hasNext();){
////			System.out.println(fqp.next());
////		}
////		System.err.println("Corrected "+changed+" sequences");
////	}
//	
//	
//	public static void kmerCoverage_count(ArrayList<String> fastqFiles,KmerSet_binary kmers) throws Exception{
//		fastqParser fqp;
//		FastqSeq fqs;
//		BufferedWriter out;
//		for (String fastqFile : fastqFiles) {
//			System.err.println(fastqFile);
//			fqp= new fastqParser(new BufferedReader(new FileReader(fastqFile)),"");
//			out= new BufferedWriter(new FileWriter(fastqFile+".kmerCoverageCount"));
//			for(int i=0;fqp.hasNext();i++){
//				if(i%10000==0){
//					System.err.print("\t"+i+"\r");
//				}
//				fqs=fqp.next();
//				out.write(fqs.getQname()+"\t"+kmers.count(fqs.getSeq())+"\n");
//			}
//		}
//	}
//	
//	public static void kmerCoverage_extractCoveredSingle(ArrayList<String> file1, String kmerFile, int kmerPos, String outPrefix,boolean gzipped)throws Exception{
//		fastqParser fqp1;
//		BufferedWriter out1= new BufferedWriter(new FileWriter(outPrefix+"_1.fastq"));
//		
//		FastqSeq fqs1;
//		BufferedReader in= new BufferedReader(new FileReader(kmerFile));
//		String s= in.readLine();
//		String[]l= s.split("\t");
//		int kmerLength=l[kmerPos].length();
//		KmerSet_binary kmers= new KmerSet_binary(kmerLength,false);
//		kmers.addKmer(l[kmerPos]);
//		for(s=in.readLine();s!=null;s=in.readLine()){
//			if(s.length()>0){
//				l=s.split("\t");
//				if(l.length>kmerPos){
//					kmers.addKmer(l[kmerPos]);
//					kmers.addKmer(sequenceUtils.reverseComplement(l[kmerPos]));
//				}
//			}
//		}
//		in.close();
//		for(int i=0;i<file1.size();++i){
//			System.err.println(file1.get(i));
//			if(gzipped){
//				fqp1= new fastqParser(new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file1.get(i))))), "");
//			}else{
//				fqp1= new fastqParser(new BufferedReader(new FileReader(file1.get(i))), "");
//			}
//			for(;fqp1.hasNext();){
//				fqs1=fqp1.next();
//				if(kmers.covered(fqs1.getSeq())){
//					//print
//					out1.write(fqs1+"\n");
//				}
//			}
//		}
//		out1.close();
//	}
//
//	
//	public static void kmerCoverage_extractCovered(ArrayList<String> file1, ArrayList<String> file2, String kmerFile, int kmerPos, String outPrefix,boolean gzipped)throws Exception{
//		fastqParser fqp1;
//		fastqParser fqp2;
//		BufferedWriter out1= new BufferedWriter(new FileWriter(outPrefix+"_1.fastq"));
//		BufferedWriter out2= new BufferedWriter(new FileWriter(outPrefix+"_2.fastq"));
//		
//		FastqSeq fqs1,fqs2;
//		BufferedReader in= new BufferedReader(new FileReader(kmerFile));
//		String s= in.readLine();
//		String[]l= s.split("\t");
//		int kmerLength=l[kmerPos].length();
//		KmerSet_binary kmers= new KmerSet_binary(kmerLength,false);
//		kmers.addKmer(l[kmerPos]);
//		for(s=in.readLine();s!=null;s=in.readLine()){
//			if(s.length()>0){
//				l=s.split("\t");
//				if(l.length>kmerPos){
//					kmers.addKmer(l[kmerPos]);
//					kmers.addKmer(sequenceUtils.reverseComplement(l[kmerPos]));
//				}
//			}
//		}
//		in.close();
//		for(int i=0;i<file1.size();++i){
//			System.err.println(file1.get(i)+"\t"+file2.get(i));
//			if(gzipped){
//				fqp1= new fastqParser(new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file1.get(i))))), "");
//				fqp2= new fastqParser(new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file2.get(i))))), "");
//			}else{
//				fqp1= new fastqParser(new BufferedReader(new FileReader(file1.get(i))), "");
//				fqp2= new fastqParser(new BufferedReader(new FileReader(file2.get(i))), "");
//			}
//			for(;fqp1.hasNext()&&fqp2.hasNext();){
//				fqs1=fqp1.next();
//				fqs2=fqp2.next();
//				if(kmers.covered(fqs1.getSeq()) || kmers.covered(fqs2.getSeq())){
//					//print
//					out1.write(fqs1+"\n");
//					out2.write(fqs2+"\n");
//				}
//			}
//		}
//		out1.close();
//		out2.close();
//	}
//	
//	public static void kmerCoverage(ArrayList<String> fastqFiles,KmerSet_binary kmers)throws Exception{
//		System.err.println("Reading background file...");
////		KmerSet kmers= (KmerSet)(new ObjectInputStream(new FileInputStream(kmerBackgroundFile))).readObject();
//		System.err.println("Parsing sequence file...");
//		fastqParser fqp;
//		FastqSeq fqs;
//		int[] cov;
//		BufferedWriter out;
//		for (String fastqFile : fastqFiles) {
//			System.err.println(fastqFile);
//			fqp= new fastqParser(new BufferedReader(new FileReader(fastqFile)),"");
//			out= new BufferedWriter(new FileWriter(fastqFile+".kmerCoverage"));
//			for(int i=0;fqp.hasNext();i++){
//				if(i%10000==0){
//					System.err.print("\t"+i+"\r");
//				}
//				fqs=fqp.next();
//				if(fqs.length()>=kmers.getKmerSize()){
//					out.write(">"+fqs.getQname()+"\n");
//					cov=kmers.coverage(fqs.getSeq());
//					out.write(""+cov[0]);
//					for(int j=1;j<cov.length;j++){
//						out.write(" "+cov[j]);
//					}
//					out.write("\n");
//				}
//			}
//			out.close();
//		}
//		System.err.println("\nDone!");
//	}
//	
//	public static int[] kmerCoverage_coverage(String seq,KmerSet_binary kmers){
//		int kmerSize=kmers.getKmerSize();
//		int[] cov=new int[seq.length()];
//		boolean[] goodKmer= new boolean[seq.length()-kmerSize+1];
//		String kmer=seq.substring(0, kmerSize);
//		goodKmer[0]=kmers.exists(kmer);
//		cov[0]=goodKmer[0]?1:0;
//		for(int i=1,j=kmerSize;j<seq.length();i++,j++){
//			kmer=kmer.substring(1)+seq.charAt(j);
//			goodKmer[i]=kmers.exists(kmer);
//			cov[i]=cov[i-1]+(goodKmer[i]?1:0)-(i>=kmerSize?(goodKmer[i-kmerSize]?1:0):0);
//		}
//		for(int i=seq.length()-kmerSize+1;i<seq.length();i++){
//			cov[i]=cov[i-1]-(goodKmer[i-kmerSize]?1:0);
//		}
//		return cov;
//	}
//	
//	
//	public static void storeKmerSet(String kmerCountFile,int cutoff,String outfile)throws Exception{
//		KmerSet_binary kmers= kmerFunctions.readKmers(kmerCountFile, cutoff);
//		System.err.println("\nwriting to file...");
//		ObjectOutputStream out= new ObjectOutputStream(new FileOutputStream(outfile));
//		out.writeObject(kmers);
//		out.close();
//	}
//	
//	public static void kmerize(BufferedReader fastqFile, int n)throws Exception{
//		fastqParser fqp=new fastqParser(fastqFile,"");
//		for(;fqp.hasNext();){
//			char[] seq=fqp.next().getSeq().toCharArray();
//			StringBuffer kmer=new StringBuffer(" "+fqp.getSeq().substring(0, n-1));
//			for(int i=n-1;i<seq.length;++i){
//				kmer.deleteCharAt(0);
//				kmer.append(seq[i]);
//				//System.out.println(kmer+" O");
//				System.out.println(kmerFunctions.kmerToUse(kmer));
//			}
//		}
//	}
//
//	public static void kmerCount(String fastqFile, int kmerSize)throws Exception{
//		KmerMap_binary kMap= new KmerMap_binary(31);
//		fastqParser fqp=new fastqParser(new BufferedReader(new FileReader(fastqFile)),"");
//		for(int i=1;fqp.hasNext();++i){
//			if(i%100000==0){
//				System.err.println("\t"+i+"\t"+kMap.size());
//			}
//			kMap.addSequenceBuffered(fqp.next().getSeq());
//		}
//		kMap.flushBuffer();
//		System.err.println("Printing...");
//		kMap.print();
//	}
//	/**
//	public static void kmerCount(String fastqFile, int kmerSize,String tmpPrefix,String outFile)throws Exception{
//		fastqParser fqp=new fastqParser(new BufferedReader(new FileReader(fastqFile)),"");
//		KmerTable kTable= new KmerTable(kmerSize,tmpPrefix,1000000);
//		for(int i=1;fqp.hasNext();i++){
//			if(i%10000==0){
////				System.err.print("\t"+i+"\t"+kTable.size()+"\r");
//			}
//			kTable.addSequence(fqp.next().getSeq());
//		}
//		System.err.println();
//		kTable.merge(outFile);
////		kTable.printHistogram();
//	}
//	*/
//	
//	public static void kmerCount_merge(ArrayList<String> files,String outFile)throws Exception{
//		KmerTable kTable= new KmerTable(files);
//		kTable.merge(outFile);
//	}
//	
//	public static void kmerCount_histogram(String kmerFile,int bin)throws Exception{
//		ObjectInputStream in= new ObjectInputStream(new FileInputStream(kmerFile));
//		HashMap<Integer, Integer>counts= new HashMap<Integer, Integer>();
//		int max=0,curBin;
//		Kmer_count_tuple kmer;
//		try{
//			for(int i=0;true;++i){
//				if(i%10000==0){
//					System.err.print("\t"+i+"\r");
//				}
//				kmer=(Kmer_count_tuple) in.readObject();
//				curBin=kmer.getCount()/bin;
//				max=curBin>max?curBin:max;
//				if(counts.containsKey(curBin)){
//					counts.put(curBin, counts.get(curBin)+1);
//				}else{
//					counts.put(curBin, 1);
//				}
//			}
//		}catch (EOFException e) {
//		}
//		System.err.println("printing...");
//		//print
//		System.out.println("kmer count (bin size="+bin+")\t#unique kmers");
//		for(int i=1;i<=max;i++){
//			if(counts.containsKey(i)){
//				System.out.println(i*bin+"\t"+counts.get(i));
//			}else{
//				System.out.println(i*bin+"\t0");
//			}
//		}
//	}
	
//	public static void randomVelvetPairedSubset(String fastqFile,double fraction)throws Exception{
//		fastqParser fqp= new fastqParser(new BufferedReader(new FileReader(fastqFile)),"");
//		Random rg= new Random();
//		for(;fqp.hasNext();){
//			if(rg.nextDouble()<fraction){
//				System.out.println(fqp.next());
//				System.out.println(fqp.next());
//			}else{
//				fqp.next();
//				fqp.next();
//			}
//		}
//	}
	
	private static void convertNumberToPhredQual(String fastqFile)throws Exception{
		fastqParser fqp= new fastqParser(new BufferedReader(new FileReader(fastqFile)),"");
		FastqSeq fqs;
		for(;fqp.hasNext();){
			fqs=fqp.next();
			fqs.convertNumberToPhredQual();
			System.out.println(fqs.toString());
		}
	}
	/**
	* @deprecated
	*/
	public static String getPrefix(String fastqFile)throws Exception{
		return getPrefix(fastqFile, 5);
	}
	
	/**
	* @deprecated
	*/
	public static String getPrefix(String fastqFile,int prefixSize)throws Exception{
		fastqParser fqp;
		FastqSeq fqs;
		String prefix="";
		boolean done=false;
		for (int i=prefixSize;i>-1&&!done;--i){
			fqp= new fastqParser(new BufferedReader(new FileReader(fastqFile)),"");
			fqs=fqp.next();
			prefix=fqs.getQname().substring(0,i);
			done=true;
			for(int j=0;j<10&&fqp.hasNext();++j){
				fqs=fqp.next();
				if(!prefix.equals(fqs.getQname().subSequence(0, i))){
//					System.out.println(prefix+"\t"+fqs.getQname().subSequence(0, i));
					done=false;
				}
			}
		}
		return prefix;
	}
	
//	private static void removeShortVelvetPair(String fastqFile, int n,String singletFile,String pairFile)throws Exception{
//		fastqParser fqp= new fastqParser(new BufferedReader(new FileReader(fastqFile)),"");
//		FastqSeq fqs1,fqs2;
//		BufferedWriter singlet= new BufferedWriter(new FileWriter(singletFile));
//		BufferedWriter pair= new BufferedWriter(new FileWriter(pairFile));
//		for(;fqp.hasNext();){
//			fqs1=fqp.next();
//			if(fqp.hasNext()){
//				fqs2=fqp.next();
//				if(fqs1.length()<n){
//					if(fqs2.length()>=n){
//						singlet.write(fqs2.toString()+"\n");
//					}
//				}else{
//					if(fqs2.length()<n){
//						singlet.write(fqs1.toString()+"\n");
//					}else{
//						pair.write(fqs1.toString()+"\n");
//						pair.write(fqs2.toString()+"\n");
//					}
//				}
//			}else{
//				singlet.close();
//				pair.close();
//				throw new Exception("File does not contain an even number of reads... last read: "+fqs1.getQname());
//			}
//		}
//		singlet.close();
//		pair.close();
//	}
	
	private static void removeShort(String fastqFile, int n)throws Exception{
		fastqParser fqp= new fastqParser(new BufferedReader(new FileReader(fastqFile)),"");
		FastqSeq fqs;
		for(;fqp.hasNext();){
			fqs= fqp.next();
			if(fqs.length()>=n){
				System.out.println(fqs.toString());
			}
		}
	}
	
	private static void reverseComplement(BufferedReader fastq)throws Exception{
		fastqParser fqp= new fastqParser(fastq,"");
		for(;fqp.hasNext();){
			System.out.println(fqp.next().reverseComplement().toString());
		}
	}
	
//	private static void clean(String fastqFile,String fastaFile,int minLength,int wordSize, String word, int maxNs) throws Exception{
//		fastqParser fqp= new fastqParser(new BufferedReader(new FileReader(fastqFile)),"");
//		fastaParser fp= new fastaParser(new BufferedReader(new FileReader(fastaFile)));
//		FastqSeq fqs;
//		FastaSeq fs;
////		ArrayList<Pattern> toMask= new ArrayList<Pattern>();
//		ArrayList<Matcher> Masks= new ArrayList<Matcher>();
//		
//		for(int i=0;i<word.length();i++){
//			Masks.add(Pattern.compile(word.charAt(i)+"{"+wordSize+",}").matcher(""));
//		}
//		MatchResult mr;
//		ArrayList<int[]> toMask;
//		Comparator<int[]> pairSorter= new IntPairComparator();
//		int printed,start,end,nCount;
//		//assumes all reads are present in both files
//		for(;fqp.hasNext()&&fp.hasNext();){
//			fqs=fqp.next();
//			fs=fp.next();
//			if(!fs.getQname().equals(fqs.getQname())){
//				System.err.println("Unmatched names - fasta: "+fs.getQname()+", fastq: "+fqs.getQname());
//			}
//			nCount=0;
//			for (int i=0;i<fs.getSeq().length();i++) {
//				if(fs.getSeq().charAt(i)=='N'){
//					nCount++;
//				}
//			}
//			if(nCount<=maxNs){
//				toMask= new ArrayList<int[]>();
//				toMask.add(new int[]{0,0});
//				toMask.add(new int[]{fs.length(),fs.length()});
//				for (Matcher matcher : Masks) {
//					matcher.reset(fs.getSeq());
//					while(matcher.find()){
//						mr=matcher.toMatchResult();
//						toMask.add(new int[]{mr.start(),mr.end()});
//					}
//				}
//				Collections.sort(toMask, pairSorter);
//				printed=0;
//				for(int i=1;i<toMask.size();i++){
//					start=toMask.get(i-1)[1];
//					end=toMask.get(i)[0];
//					if(end-start>minLength){
//						//print
//						printed++;
//						if(printed==1){
//							//without change to name
//							System.out.println(fqs.toString(start, end));
//						}else{
//							//with added idenifier
//							System.out.println(fqs.toString(start, end, fqs.getQname()+"_"+printed));
//						}
//					}
//				}
//			}
//		}
//	}
//	
//	private static void randomsub(String fastqFile,int N)throws Exception{
////		fastqParser fqp= new fastqParser(new BufferedReader(new FileReader(fastqFile)),getPrefix(fastqFile));
//		FastqParser_nio fqp= new FastqParser_nio(fastqFile, "");
//		Random rs= new Random();
//		final int luckynumber=rs.nextInt(N);
//		for (;fqp.hasNext();){
//			if(rs.nextInt(N)==luckynumber){
//				System.out.println(fqp.next().toString());
//			}else{
//				fqp.next();
//			}
//		}
//	}
	
	private static void sub(ArrayList<BufferedReader> fastqFiles,BufferedReader subsetFile,boolean inSubset) throws Exception{
		HashSet<String> subset= fileToHashset(subsetFile);
		for (BufferedReader fastqFile : fastqFiles) {
			sub(fastqFile,subset,inSubset);	
		}
	}
	
	private static void subToFile(ArrayList<String> fastqFiles,BufferedReader subsetFile,boolean inSubset)throws Exception{
		HashSet<String> subset= fileToHashset(subsetFile);
		for (String fastqFile : fastqFiles) {
			subToFile(fastqFile,subset,inSubset,fastqFile+"sub.fq");	
		}
	}
	
	private static void subToFile(String fastqFile, HashSet<String> subset,boolean inSubset,String outFile)throws Exception{
		fastqParser fqp= new fastqParser(new BufferedReader(new FileReader(fastqFile)),"");
		FastqSeq fqs;
		BufferedWriter out= new BufferedWriter(new FileWriter(outFile));
		for(;fqp.hasNext();){
			fqs=fqp.next();
			if(subset.contains(fqs.getQname())==inSubset){
				out.write(fqs.toString()+"\n");
			}
		}
		out.close();
	}
	
//	private static void sub(String fastqFile,String subsetFile, boolean inSubset)throws Exception{
//		sub(fastqFile,fileToHashset(subsetFile),inSubset);
		
//		fastqParser fqp= new fastqParser(new BufferedReader(new FileReader(fastqFile)),getPrefix(fastqFile));
//		FastqSeq fqs;
//		
//		for(;fqp.hasNext();){
//			fqs=fqp.next();
//			if(subset.contains(fqs.getQname())==inSubset){
//				System.out.println(fqs.toString());
//			}
//		}
//	}
	
	
	private static HashSet<String> fileToHashset(BufferedReader in)throws Exception{
		HashSet<String> hashset= new HashSet<String>();
//		BufferedReader in= new BufferedReader(new FileReader(file));
		for(String s= in.readLine();s!=null;s=in.readLine()){
			hashset.add(s);
		}
		in.close();
		return hashset;
	}
	
	private static void sub(BufferedReader fastqFile, HashSet<String> subset, boolean inSubset)throws Exception{
		fastqParser fqp= new fastqParser(fastqFile,"");
		FastqSeq fqs;
		
		for(;fqp.hasNext();){
			fqs=fqp.next();
			if(subset.contains(fqs.getQname())==inSubset){
				System.out.println(fqs.toString());
			}
		}
	}
	
//	private static void subJia(String fastqFile,String subsetFile, boolean inSubset,int dropPost)throws Exception{
//		HashSet<String> subset= fileToHashset(subsetFile);
//		
//		fastqParser fqp= new fastqParser(new BufferedReader(new FileReader(fastqFile)),"");
//		FastqSeq fqs;
//		
//		for(;fqp.hasNext();){
//			fqs=fqp.next();
//			if(subset.contains(fqs.getQname().substring(0, fqs.getQname().length()-dropPost))==inSubset){
//				System.out.println(fqs.toString());
//			}
//		}
//	}
	
	private static void subReg(String fastqFile, String regExpString, boolean inSubset) throws Exception{
		Matcher regExpMatcher=Pattern.compile(regExpString).matcher("");
		fastqParser fqp= new fastqParser(new BufferedReader(new FileReader(fastqFile)),"");
		FastqSeq fqs;
		for(;fqp.hasNext();){
			fqs=fqp.next();
			regExpMatcher.reset(fqs.getQname());
			if(regExpMatcher.find()==inSubset){
				System.out.println(fqs.toString());
			}
		}
	}
	
//	private static void information(String fastqFile,int burnin,int sampleSize,int sampleDistance,int minWordLength, int maxWordLength,String outPrefix)throws Exception{
//		fastqParser fqp= new fastqParser(new BufferedReader(new FileReader(fastqFile)),"");
//		FastqSeq fqs;
//		StringContent seqInfo= new StringContent();
////		StringContent qual= new StringContent(minWordLength);
//		final int lastSample= sampleSize==0?Integer.MAX_VALUE:burnin+sampleSize*sampleDistance+1; 
//		for(int count=1;fqp.hasNext()&&count<lastSample;count++){
//			if(count>burnin && (count-burnin)%sampleDistance==0){
//				//sample
//				fqs=fqp.next();
//				final String seq=fqs.getSeq();
//				final int length=seq.length();
//				for(int wordLength=minWordLength;wordLength<maxWordLength;wordLength++){
//					for(int i=0;i<length-wordLength;i++){
//						seqInfo.add(seq.substring(i, i+wordLength));
//					}
//				}
//			}else{
//				//discard
//				fqp.next();
//			}
//		}
//		//print
//		for(int wordLength=minWordLength;wordLength<maxWordLength;wordLength++){
//			BufferedWriter out= new BufferedWriter(new FileWriter(outPrefix+"_wordLength_"+wordLength+".csv"));
//			out.write("#Total"+sep+seqInfo.getTotal(wordLength+"")+"\nWord"+sep+"Count\n");
//			HashMap<String, Integer> dist=seqInfo.getCount(wordLength+"");
//			for(String key : dist.keySet()){
//				out.write(key+sep+dist.get(key)+"\n");
//			}
//			out.close();
//		}
//	}
}

class IntPairComparator implements Comparator<int[]>{

	public int compare(int[] arg0, int[] arg1) {
		//assumes that the arrays will be of size two, but at least one.
		return arg0[0]-arg1[0];
	}
	
}

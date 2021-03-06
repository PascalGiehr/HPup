package tools.fasta;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.LinkedList;
import java.util.Random;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPOutputStream;

//import tools.blast.blastM8Alignment;
//import tools.blast.blastM8Parser;
import tools.fastq.FastqSeq;
//import tools.gff.gffLine_startComparator;
//import tools.hmmer.hmmerAlignment;
//import tools.hmmer.hmmerModelAlignment;
//import tools.hmmer.hmmerParser;
//import tools.kmer.KmerMap_binary;
//import tools.kmer.KmerSet_binary;
//import tools.kmer.kmerFunctions;
//import tools.nexus.NexusSet;
//import tools.pfam.pfamAlignment;
//import tools.pfam.pfamParser;
//import tools.phylip.PhylipSet;
import tools.rocheQual.RocheQualSeq;
import tools.rocheQual.RocheQualParser;
import tools.sequences.sequenceUtils;
import tools.utils.general;
import tools.utils.multiPatternMatcher.multipatternTree;
//import tools.vcf.bisSNP.BisSNPbedEntry;
//import tools.vcf.bisSNP.BisSNPbedEntryParser;

public class fastaUtils {


	final static String sep="\t";
	static HashSet<String> nucleotides;
	/**
	 * @param args
	 */
	public static void main(String[] args)throws Exception {
		//args=new String[]{"clean","/local/membrane_protein/","/local/slumpen/f3strict/f3_ncRNA_conservative.txt"};
		nucleotides=new HashSet<String>();
		nucleotides.add("A");
		nucleotides.add("T");
		nucleotides.add("G");
		nucleotides.add("C");
		nucleotides.add("N");
		if(args.length==0){
			System.out.println(printHelp());
		}else{
			if(args[0].equals("sub")){
				if(args.length==3)
					subSet(args[1],new BufferedReader(new FileReader(args[2])),true);
				else{
					System.err.println(printHelp());
				}
			}else if(args[0].equals("subPipe")){
				if(args.length==2)
					subSet(args[1],new BufferedReader(new InputStreamReader(System.in)),true);
				else
					System.err.println(printHelp());
			}else if(args[0].equals("sort")){
				if(args.length==3)
					sort(args[1],args[2]);
				else
					System.err.println(printHelp());
			}else if(args[0].equals("orfStat")){
				if(args.length==3)
					if(args[2].equals("forward")){
						orfStat(args[1],false);
					}else if(args[2].equals("both")){
						orfStat(args[1],true);
					}
				else
					System.err.println(printHelp());
			}else if(args[0].equals("antisub")){
				if(args.length==3)
					subSet(args[1],new BufferedReader(new FileReader(args[2])),false);
				else
					System.err.println(printHelp());
			}else if(args[0].equals("subReg")){
				if(args.length==3)
					subSetReg(args[1],args[2],true);
				else
					System.err.println(printHelp());
			}else if(args[0].equals("antisubReg")){
				if(args.length==3)
					subSetReg(args[1],args[2],false);
				else
					System.err.println(printHelp());
			}else if(args[0].equals("kmerizeBis")){
				if(args.length==3)
					kmerizeBis(args[1],Integer.parseInt(args[2]));
				else
					System.err.println(printHelp());
//			}else if(args[0].equals("kmerize")){
//				if(args.length==3)
//					kmerize(new BufferedReader(new FileReader(args[1])),Integer.parseInt(args[2]));
//				else
//					System.err.println(printHelp());
			}else if(args[0].equals("kmerizeNoRC")){
				if(args.length==3)
					kmerizeNoRC(new BufferedReader(new FileReader(args[1])),Integer.parseInt(args[2]));
				else
					System.err.println(printHelp());
//			}else if(args[0].equals("kmerizePipe")){
//				if(args.length==2)
//					kmerize(new BufferedReader(new InputStreamReader(System.in)),Integer.parseInt(args[1]));
//				else
//					System.err.println(printHelp());
//			}else if(args[0].equals("kmerizeWithTag")){
//				if(args.length==3)
//					kmerizeWithTag(args[1],Integer.parseInt(args[2]));
//				else
//					System.err.println(printHelp());
//			}else if(args[0].equals("kmerizeWithTagNoRC")){
//				if(args.length==3)
//					kmerizeWithTagNoRC(args[1],Integer.parseInt(args[2]));
//				else
//					System.err.println(printHelp());
//			}else if(args[0].equals("kmerizeAllAlt")){
//				if(args.length==5)
//					kmerizeAllAlt(args[1],args[2],args[3],Integer.parseInt(args[4]));
//				else
//					System.err.println(printHelp());
//			}else if(args[0].equals("kmerizeAllAltPosWithTag")){
//				if(args.length==4)
//					kmerizeAllAltPosWithTag(args[1],Integer.parseInt(args[2]),Integer.parseInt(args[3]));
//				else
//					System.err.println(printHelp());
			}else if(args[0].equals("allAltPos")){
				if(args.length==4)
					allAltPos(args[1],Integer.parseInt(args[2]));
				else
					System.err.println(printHelp());
//			}else if(args[0].equals("kmerStat")){
//				if(args.length==4)
//					kmerStat(args[1],args[2],Integer.parseInt(args[3]));
//				else
//					System.err.println(printHelp());
			}else if(args[0].equals("kmerizeWRevComp")){
				if(args.length==3)
					kmerizeWRevComp(args[1],Integer.parseInt(args[2]));
				else
					System.err.println(printHelp());
//			}else if(args[0].equals("domains")){
//				if(args.length==4)
//					getDomains(args[1],args[2],args[3],true);
//				else
//					System.err.println(printHelp());
//			}else if(args[0].equals("pfamDomains")){
//				if(args.length==4)
//					pfamDomains(args[1],args[2],args[3],true,Integer.MAX_VALUE);
//				else if(args.length==5)
//					pfamDomains(args[1],args[2],args[3],true,Integer.parseInt(args[4]));
//				else
//					System.err.println(printHelp());
//			}else if(args[0].equals("antidomains")){
//				if(args.length==4)
//					getDomains(args[1],args[2],args[3],false);
//				else
//					System.err.println(printHelp());
//			}else if(args[0].equals("clean")){
//				if(args.length==4)
//					domain_clean(args[1],args[2],args[3],true);
//				else
//					System.err.println(printHelp());
//			}else if(args[0].equals("anticlean")){
//				if(args.length==4)
//					domain_clean(args[1],args[2],args[3],false);
//				else
//					System.err.println(printHelp());
			}else if(args[0].equals("split")){
				if(args.length==4) {
					split(args[1],args[2],Integer.parseInt(args[3]));
				}else
					System.err.println(printHelp());
			}else if(args[0].equals("description")){
				if(args.length==3)
					description(args[1],args[2]);
				else if(args.length==2)
					description(args[1]);
				else
					System.err.println(printHelp());
			}else if(args[0].equals("cleanse")){
				if(args.length==2)
					cleanse(args[1]);
				else
					System.err.println(printHelp());
			}else if(args[0].equals("cleanseShort")){
				if(args.length==2)
					cleanseShort(args[1]);
				else
					System.err.println(printHelp());
			}else if(args[0].equals("phylipWarning")){
				if(args.length==2)
					phylipWarning(args[1]);
				else
					System.err.println(printHelp());
			}else if(args[0].equals("tag")){
				if(args.length==3)
					tag(new BufferedReader(new FileReader(args[1])),args[2]);
				else
					System.err.println(printHelp());
			}else if(args[0].equals("tagPipe")){
				if(args.length==2)
					tag(new BufferedReader(new InputStreamReader(System.in)),args[1]);
				else
					System.err.println(printHelp());
			}else if(args[0].equals("dropShort")){
				if(args.length==3)
					dropShort(args[1], Integer.parseInt(args[2]));
				else
					System.err.println(printHelp());
			}else if(args[0].equals("chance")){
				if(args.length==3)
					chance(args[1], Integer.parseInt(args[2]));
				else
					System.err.println(printHelp());
			}else if(args[0].equals("diffTag")){
				if(args.length==3)
					diffTag(args[1], args[2]);
				else
					System.err.println(printHelp());
			}else if(args[0].equals("lengths")){
				if(args.length==2)
					lengths(args[1]);
				else
					System.err.println(printHelp());
			}else if(args[0].equals("countTagged")){
				if(args.length==2)
					countTagged(args[1]);
				else
					System.err.println(printHelp());
			}else if(args[0].equals("cutColumns")){
				if(args.length==3)
					cutColumns(args[1],Double.parseDouble(args[2]));
				else
					System.err.println(printHelp());
			}else if(args[0].equals("removeDash")){
				if(args.length==2)
					removeDash(args[1]);
				else
					System.err.println(printHelp());
//			}else if(args[0].equals("toPhy")){
//				if(args.length==2)
//					toPhy(args[1]);
//				else
//					System.err.println(printHelp());
//			}else if(args[0].equals("toNexus")){
//				if(args.length==2)
//					toNexus(args[1]);
//				else
//					System.err.println(printHelp());
			}else if(args[0].equals("toStockholm")){
				if(args.length==2)
					toStockholm(args[1]);
				else
					System.err.println(printHelp());
			}else if(args[0].equals("shortenNames")){
				if(args.length==3)
					shortenNames(args[1],args[2]);
				else
					System.err.println(printHelp());
			}else if(args[0].equals("reduceTaggedByChance")){
				if(args.length==3)
					reduceTaggedByChance(args[1],Integer.parseInt(args[2]));
				else
					System.err.println(printHelp());
			}else if(args[0].equals("subFiles")){
				if(args.length==4)
					subFiles(args[1],args[2],args[3]);
				else
					System.err.println(printHelp());
			}else if(args[0].equals("giList")){
				if(args.length==2)
					giList(args[1]);
				else
					System.err.println(printHelp());
			}else if(args[0].equals("removeDesc")){
				if(args.length==2)
					removeDesc(args[1]);
				else
					System.err.println(printHelp());
			}else if(args[0].equals("getEditedRNA")){
				if(args.length==2)
					getEditedRNA(args[1]);
				else
					System.err.println(printHelp());
			}else if(args[0].equals("extractSubSeq")){
				if(args.length==3)
					extractSubSeq(args[1], general.getBufferedReader(args[2]),0);
				else if(args.length==4)
					extractSubSeq(args[1], general.getBufferedReader(args[2]), Integer.parseInt(args[3]));
				else
					System.err.println(printHelp());
			}else if(args[0].equals("cut")){
				if(args.length==3)
					cut(args[1],args[2]);
				else
					System.err.println(printHelp());
			}else if(args[0].equals("cleanXs")){
				if(args.length==5)
					cleanXs(args[1],args[2],Integer.parseInt(args[3]),args[4]);
				else
					System.err.println(printHelp());
//			}else if(args[0].equals("renameFromM8blast")){
//				if(args.length==3)
//					renameFromM8blast(args[1],args[2]);
//				else
//					System.err.println(printHelp());
			}else if(args[0].equals("genomicPrimerSite")){
				if(args.length==2)
					genomicPrimerSite(args[1]);
				else
					System.err.println(printHelp());
			}else if(args[0].equals("rename")){
				if(args.length==3)
					rename(args[1],args[2]);
				else
					System.err.println(printHelp());
			}else if(args[0].equals("oneLetterPerLine")){
				if(args.length==2)
					oneLetterPerLine(args[1]);
				else
					System.err.println(printHelp());
			}else if(args[0].equals("alignedaa2rna_shiftFrame")){
				if(args.length==3)
					alignedaa2rna_shiftFrame(args[1],args[2]);
				else
					System.err.println(printHelp());
			}else if(args[0].equals("drop2andfilter")){
				if(args.length==2)
					drop2andfilter(args[1]);
				else
					System.err.println(printHelp());
			}else if(args[0].equals("countPos")){
				if(args.length==2)
					countPos(args[1]);
				else
					System.err.println(printHelp());
			}else if(args[0].equals("culmativeXcount")){
				if(args.length==2)
					culmativeXcount(args[1]);
				else
					System.err.println(printHelp());
			}else if(args[0].equals("GCcontentPerRow")){
				if(args.length==3)
					GCcontentPerRow(args[1],Integer.parseInt(args[2]));
				else
					System.err.println(printHelp());
			}else if(args[0].equals("motifFind") && args.length==3){
				motifFind(args[1],args[2]);
			}else if(args[0].equals("motifFindReg") && args.length==3){
				motifFindReg(args[1],args[2]);
			}else if(args[0].equals("randomCut")){
				if(args.length==3)
					randomCut(args[1],Integer.parseInt(args[2]));
				else
					System.err.println(printHelp());
			}else if(args[0].equals("alignedaa2rna")){
				if(args.length==3)
					alignedaa2rna(args[1],args[2],true);
				else if(args.length==4 && args[3].toUpperCase().equals("F")){
					alignedaa2rna(args[1],args[2],false);
				}else
					System.err.println(printHelp());
//			}else if(args[0].equals("cContext")&&args.length>3){
//				if(args[2].endsWith("gz")){
//					System.err.println("No support for gzipped vcf files due to problems with bgzip");
//				}
//				ArrayList<String> context= new ArrayList<String>();
//				for(int i=3;i<args.length;++i){
//					context.add(args[i]);
//				}
//				cContext(args[1],general.getBufferedReader(args[2]),context);
			}else if(args[0].equals("reverseComplementGenome")&&args.length==2){
				reverseComplementGenome(args[1]);
//			}else if(args[0].equals("kmerContains")&&args.length>3){
//				kmerContains(args[1], args[2],Integer.parseInt(args[3]));
//			}else if(args[0].equals("kmerLinks")&&args.length>2){
//				kmerLinks(new BufferedReader(new FileReader(args[1])), Integer.parseInt(args[2]));
			}else if(args[0].equals("MSAtoVistaPlot")&&args.length>1){
				MSAtoVistaDiagram(args[1], args[4], Integer.parseInt(args[2]), true, args[3],false);
			}else if(args[0].equals("MSAtoVistaPlotF")&&args.length>1){
				MSAtoVistaDiagram(args[1], args[4], Integer.parseInt(args[2]), false, args[3],false);
			}else if(args[0].equals("MSAtoVistaPlotM")&&args.length>1){
				MSAtoVistaDiagram(args[1], args[4], Integer.parseInt(args[2]), false, args[3],true);
			}else if(args[0].equals("splitContigs")&&args.length>3){
				splitContigs(args[1],Integer.parseInt(args[2]),Integer.parseInt(args[3]));
			}else if(args[0].equals("splitVelvetContigs")&&args.length>2){
				splitVelvetContigs(args[1],Integer.parseInt(args[2]));
			}else if(args[0].equals("splitOnN")&&args.length>3){
				splitOnN(args[1],Integer.parseInt(args[2]),Integer.parseInt(args[3]));
			}else if(args[0].equals("splitNewblerContigs")&&args.length>3){
				splitNewblerContigs(args[1],args[2],Integer.parseInt(args[3]));
			}else if(args[0].equals("nCount")&&args.length>1){
				nCount(args[1]);
			}else if(args[0].equals("extractPaired")&&args.length>3){
				extractPaired(new BufferedReader(new FileReader(args[1])),new BufferedReader(new FileReader(args[2])),args[3]);
			}else if(args[0].equals("extractPairedNoOrder")&&args.length>4){
				extractPairedNoOrder(new BufferedReader(new FileReader(args[1])),new BufferedReader(new FileReader(args[2])),args[3],Integer.parseInt(args[4]));
			}else if(args[0].equals("compare")&&args.length>1){
				compare(args[1],args[2]);
			}else if(args[0].equals("count")&&args.length>1){
				count(args[1],args[2].charAt(0));
//			}else if(args[0].equals("cutAtBlastStart")&&args.length>1){
//				cutAtBlastStart(args[1],args[2]);
			}else if(args[0].equals("tails")&&args.length>2){
				tails(args[1],Integer.parseInt(args[2]));
//			}else if(args[0].equals("tailsTable")&&args.length>2){
//				tailsTable(args[1],Integer.parseInt(args[2]));
//			}else if(args[0].equals("getMismatch")&&args.length>2){
//				getMismatch(args[1],args[2],args[3]);
//			}else if(args[0].equals("kmerCoverage_extractPerfectDirect")&&args.length>3){
//				ArrayList<String> fastaFiles= new ArrayList<String>();
//				for(int i=3;i<args.length;++i){
//					fastaFiles.add(args[i]);
//				}
//				kmerCoverage_extractPerfectDirect(kmerFunctions.readKmers(args[1], Integer.parseInt(args[2])), fastaFiles);
//			}else if(args[0].equals("kmer_frequency")&&args.length>2){
//				kmer_frequency(args[1],args[2]);
//			}else if(args[0].equals("mask")&&args.length>2){
//				mask(args[1],args[2]);
			}else if(args[0].equals("minLengthPipe")&&args.length>1){
				minLength(new BufferedReader(new InputStreamReader(System.in)),Integer.parseInt(args[1]));
			}else if(args[0].equals("pad")&&args.length>2){
				pad(args[1],Integer.parseInt(args[2]));
			}else if(args[0].equals("reverseComplement")&&args.length>1){
				reverseComplementFile(args[1]);
//			}else if(args[0].equals("removeDoublets")&&args.length>1){
//				removeDoublets(new BufferedReader(new FileReader(args[1])));
//			}else if(args[0].equals("removeDoubletsPipe")&&args.length>0){
//				removeDoublets(new BufferedReader(new InputStreamReader(System.in)));
			}else if(args[0].equals("toFastq")){
				if(args.length==3)
					toFastq(new BufferedReader(new FileReader(args[1])),args[2]);
				else if(args.length==2)
					toFastq(new BufferedReader(new FileReader(args[1])));
				else
					System.err.println(printHelp());
			}else if(args[0].equals("toFastqPipe")&&args.length==1){
				toFastq(new BufferedReader(new InputStreamReader(System.in)));
			}else if(args[0].equals("toFastqRoche")&&args.length>2){
				toFastqRoche(args[1],args[2]);
			}else if(args[0].equals("vistaPairsToPlot")&&args.length>1){
				ArrayList<String> pairFiles= new ArrayList<String>();
				for(int i=1;i<args.length;i++){
					pairFiles.add(args[i]);
				}
				vistaPairsToPlot(pairFiles);
			}else{
				System.err.println("end");
				System.err.println(printHelp());
				System.err.println("end");
			}
		}
	}

	private static String printHelp(){
		String help="First argument is the method to use (sub,antisub)\n\n";
		help+="countTagged - counts the number of transcripts in each tagged group\n";
		help+="\tArguments: fastaUtils countTagged <tagged fafile>\n";
		help+="cleanXs - reads a fasta and a qual file. Searches for longer stretches of X, split and prints all which are longer than minLength\n";
		help+="\tArguments: fastaUtils cleanXs <fafile> <qual file> <minLength> <outprefix>\n";
		help+="getEditedRNA - prints the rna sequence for each file in the filelist\n";
		help+="\tArguments: fastaUtils getEditedRNA <filelist>\n";
		help+="reduceTaggedByChance - takes n transcripts from each tagged group (>tag_id) by chance\n";
		help+="\tArguments: fastaUtils reduceTaggedByChance <fafile> <n>\n";
		help+="antisub - extracts a subset, defined as those not present in a file with one entry on each line, to a new fa-file printed as output.\n";
		help+="\tArguments: fastaUtils antisub <fafile> <subSetFile>\n";
		help+="subReg - extracts a subset, defined by a regular expression, to a new fa-file printed as output.\n";
		help+="\tArguments: fastaUtils sub <fafile> <regExp>\n";
		help+="oneLetterPerLine - Prints the fasta file with one letter per line\n";
		help+="\tArguments: fastaUtils oneLetterPerLine <fafile>\n";
		help+="antisubReg - extracts a subset, defined as those not matching the regular expression, to a new fa-file printed as output.\n";
		help+="\tArguments: fastaUtils antisub <fafile> <regExp>\n";
		help+="drop2andfilter - drops the first two letters of the sequence and then all third position letters.\n";
		help+="\tArguments: fastaUtils drop2andfilter <fafile>\n";
		help+="domains - extracts the sequences for the domains in the list\n";
		help+="\tArguments: fastaUtils domains <fafile> <hmmerFile> <domainSet>\n";
		help+="pfamDomains - extracts the sequences for the domains in the list\n";
		help+="\tArguments: fastaUtils domains <fafile> <pfam_scan.pl result file> <domainSet> <nr of domains to print sorted by evalue=Infinity>\n";
		help+="antidomains - extracts the sequences for the domains not in the list\n";
		help+="\tArguments: fastaUtils antidomains <fafile> <hmmerFile> <domainSet>\n";
		help+="clean - removes all domains not in the list\n";
		help+="\tArguments: fastaUtils clean <fafile> <hmmerFile> <domainSet>\n";
		help+="anticlean - removes all domains in the list\n";
		help+="\tArguments: fastaUtils anticlean <fafile> <hmmerFile> <domainSet>\n";
		help+="split - splits file into x files with at most n transcripts in each. (Not to be applied on qual files)\n";
		help+="\tArguments: fastaUtils split <fafile> <outprefix> <n>\n";
		help+="description - prints the header of the gi:s in the gi_list\n";
		help+="\tArguments: fastaUtils description <fafile> (<gi_list>)\n";
		help+="removeDesc - removes the description from the header\n";
		help+="\tArguments: fastaUtils removeDesc <fafile> \n";
		help+="cleanse - removes doublets\n";
		help+="\tArguments: fastaUtils cleanse <fafile>\n";
		help+="cleanseShort - removes doublets, keeps the longer variant\n";
		help+="\tArguments: fastaUtils cleanseShort <fafile>\n";
		help+="tag(Pipe) - adds a tag directly after >\n";
		help+="\tArguments: fastaUtils tag <fafile (Std.in)> <tag>\n";
		help+="giList - prints a list of the qnames >\n";
		help+="\tArguments: fastaUtils giList <fafile> \n";
		help+="dropShort - removes all transcripts shorter than n \n";
		help+="\tArguments: fastaUtils dropShort <fafile> <n>\n";
		help+="diffTag -  tags a fa-file according to tagfile (id\ttag)\n";
		help+="\tArguments: fastaUtils dropShort <fafile> <tagfile>\n";
		help+="phylipWarning - prints all name groups which do not have an unique begining\n";
		help+="\tArguments: fastaUtils phylipWarning <fafile>\n";
		help+="chance - select n transcripts on chance\n";
		help+="\tArguments: fastaUtils chance <fafile> <n>\n";
		help+="removeDash - clean the file from dashes (-)\n";
		help+="\tArguments: fastaUtils removeDash <fafile>\n";
		help+="toPhy - converts the fasta-file to .phy format... inserts dashes (-) in the end if the sequences are of different length\n";
		help+="\tArguments: fastaUtils toPhy <fafile>\n";
		help+="toNexus - converts the fasta-file to .nxs format... inserts dashes (-) in the end if the sequences are of different length\n";
		help+="\tArguments: fastaUtils toNexus <fafile>\n";
		help+="toStockholm - converts the fasta-file to Stockholm format...\n";
		help+="\tArguments: fastaUtils toStockholm <fafile>\n";
		help+="cutColumns - removes columns from an (aligned) fa-file witch is represented in less than covCut (0.0-1.0:% ,>1: number) sequences \n";
		help+="\tArguments: fastaUtils cutColumns <fafile> <covCut>\n";
		help+="lengths - prints the length of each sequence\n";
		help+="\tArguments: fastaUtils lengths <fafile>\n";
		help+="shortenNames - replaces the first instance of the regExp with nothing\n";
		help+="\tArguments: fastaUtils shortenNames <fafile> <regExp>\n";
		help+="cut - prints the subsequences defined by the position table file (tname\\tstart\\tstop\\tstrand\\outname)\n";
		help+="\tArguments: fastaUtils cut <fafile> <position table>\n";
		help+="randomCut - If possible extracts a section of the given size from a random location in each sequence\n";
		help+="\tArguments: fastaUtils cut <fafile> <size>\n";
		help+="sub(Pipe) - extracts a subset, defined in a file with one entry on each line, to a new fa-file printed as output.\n";
		help+="\tArguments: fastaUtils sub <fafile> <subSetFile (Std.in)>\n";
		help+="subFiles - extracts a subsets, defined in a file with one entry together with its group on line (id\\tgroup), to a new fa-files with the given prefix\n";
		help+="\tArguments: fastaUtils subFiles <fafile> <subSetFile> <outPrefix>\n";
		help+="alignedaa2rna - takes a protein fasta file (aligned or not,cut or not) and translates it according to a rna fasta file. The names must be the same. If the third parameter is set to T (F is optional), all codons will be checked to be correct.\n";
		help+="\tArguments: fastaUtils alignedaa2rna <aa file> <rna file> <check all=T (or F)>\n";
		help+="alignedaa2rna_shiftFrame - Inserts a N in the beginning of each sequence in the giList.\n";
		help+="\tArguments: fastaUtils alignedaa2rna_shiftFrame <fafile> <giList> \n";
		help+="sort - prints the sequences according to the order in the list\n";
		help+="\tArguments: fastaUtils sort <fa file> <list> \n";
		help+="renameFromM8blast - prints the sequences according to the order in the list\n";
		help+="\tArguments: fastaUtils renameFromM8blast <fa file> <blast file> \n";
		help+="rename - renames the sequences according to the giTable file\n";
		help+="\tArguments: fastaUtils rename <fa file> <giTable> \n";
		help+="orfStat - prints the lengths of the orfs in all six reading frames in the fasta sequence\n";
		help+="\tArguments: fastaUtils orfStat <fa file> <both/forward> \n";
		help+="genomicPrimerSite - reads a fasta seq and tells if there is a genomic primer site according to Majds prerequsite\n";
		help+="\tArguments: fastaUtils genomicPrimerSite <fa file> \n";
		help+="cContext - Summarizes the average coverage and methylation for all reads in a specific context. A context is defined as a IUPAC string s and a 1 based number n for which C is in question (s,n or GCG,2)\n";
		help+="\tArguments: fastaUtils cContext <fa file> <bisSNP vcf file> <context1> ... <context n> \n";
		help+="countPos - counts the number of different letters with respect to position in a fasta file\n";
		help+="\tArguments: fastaUtils countPos <fa file> \n";
		help+="culmativeXcount - take a file with screened reads and outputs a tabbed table with the format rownumber,id,isscreened (1 or 0),culmative count of isscreened\n";
		help+="\tArguments: fastaUtils culmativeXcount <fa file> \n";
		help+="GCcontentPerRow - calculates the GC-content per row in a file with only one sequence. assumes fixed row length. prints fixed wiggle format\n";
		help+="\tArguments: fastaUtils GCcontentPerRow <fa file> <row length>\n";
		help+="filterIlluminaPairsPre - Identifies sequences that shares the n first base pairs\n";
		help+="\tArguments: fastaUtils filterIlluminaPairsPre <fa file> <motif file>\n";
		help+="compare - compares two fasta files and prints mutated positions. Expects point mutations and similar lines\n";
		help+="\t<input> = <fastaFile1> <fastaFile2> \n";
		help+="count - counts the number of the given char X in a sequence \n";
		help+="\t<input> = <fastaFile> <char X> \n";
		help+="cutAtBlastStart - takes a fasta file and a blast file, cuts at the starting point of each blasthit. Meant to cut Sanger sequences at an adapter \n";
		help+="\t<input> = <fastaFile> <blastFile> \n";
		help+="extractSubSeq - prints the sub sequences defined by the position file (qname\\tstart\\tstop). Each region will optionally be flanked by x bp\n";
		help+="\tArguments: fastaUtils extractSubSeqExt <fa file> <position file>\n";
		help+="extractPaired - given a gi-list of the original order as std.in, reads present in both files will be extracted";
		help+="\tArguments: fastaUtils extractPaired <fa file 1> <fa file 2> <outPrefix>\n";
		help+="extractPairedNoOrder - given a buffer size, reads present in both files will be extracted";
		help+="\tArguments: fastaUtils extractPaired <fa file 1> <fa file 2> <outPrefix> <bufferSize>\n";
		help+="getMismatch - Extends the blast-file with the first mismatch\n";
		help+="\t<input> = <blastFile> <fastaFile1 (query)> <fastaFile2 (target)> \n";
		help+="kmerContains - prints all ids for sequencing containing the given kmers \n";
		help+="\t<input> = <fastaFile> <kmerFile> <kmerPos>\n";
		help+="kmerCoverage_extractPerfecDirect - prints all sequences with perfect kmerCoverage given the kmerCount file \n";
		help+="\t<input> = <kmerFile> <kmerPos> <fasta file 1> ... <fasta file n>\n";
		help+="kmerize(NoRC) - prints all kmers of length n in the sequences \n";
		help+="\t<input> = <fastaFile> <n>\n";
		help+="kmerizeBis - prints all kmers of length n in the sequences together with possible C to T changes. Forward and Reverse \n";
		help+="\t<input> = <fastaFile> <n>\n";
		help+="kmerizeWithTag(NoRC) - prints all kmers of length n in the sequences together with the qname it comes from \n";
		help+="\t<input> = <fastaFile> <n>\n";
		help+="kmerizeAllAlt - prints all kmers of length n in the sequences together with all posible variants for the mutation\n";
		help+="\t<input> = <blastFile> <fastaFile1 (query)> <fastaFile2 (target)> <n>\n";
		help+="kmerizeAllAltPosWithTag - prints all kmers of length n in the sequences together with all posible variants for a mutation at pos p\n";
		help+="\t<input> = <fastaFile1> <n> <p>\n";
		help+="allAltPos - prints all sequences with all possible nucleotides at pos p\n";
		help+="\t<input> = <fastaFile1> <p>\n";
		help+="kmerStat - take a kmer file and prints the total kmer-coverage for each column together with the length\n";
		help+="\t<input> = <fastafile> <kmerfile> <kmerPos>\n";
		help+="kmerizeWRevComp - prints all kmers of length n in the sequences, together with their reverse complements \n";
		help+="\t<input> = <fastaFile> <n>\n";
		help+="kmerLinks - prints all pairs of neighboring kmers given the sequence data \n";
		help+="\t<input> = <fastaFile> <n>\n";
		help+="kmer_frequency - takes a reduced kmer distribution (count\\tkmer) and prints a tab separated line for each sequence with the frequency of the occuring kmers in the distribution \n";
		help+="\t<input> = <fastaFile> <(reduced) kmer dist>\n";
		help+="mask - Replaces the regions represented in the gff file by N's\n";
		help+="\t<input> = <fastaFile> <gff file>\n";
		help+="minLengthPipe - removes sequences shorter than n\n";
		help+="\t<input> = <n>\n";
		help+="motifFind - Identifies the presence of the motifs given in the motif file among the sequences in the fasta file\n";
		help+="\tArguments: fastaUtils motifFind <fa file> <motif file>\n";
		help+="motifFindReg - Identifies the presence of the motifs given in the regular expression among the sequences in the fasta file (only on the forward strand)\n";
		help+="\tArguments: fastaUtils motifFind <fa file> <regExp>\n";
		help+="MSAtoVistaPlot - calculates a windowed count for the conservation of a multiple alignment\n";
		help+="\t<input> = <fastaFile> <windowSize> <gapSymbols> <outfile>\n";
		help+="nCount - counts the number of N in a sequence \n";
		help+="\t<input> = <fastaFile> \n";
		help+="pad - adds n N before and after each sequence \n";
		help+="\t<input> = <fastaFile> <N> \n";
		help+="removeDoublets(Pipe) - removes exact doublets\n";
		help+="\t<input> = <fastaFile>\n";
		help+="reverseComplement - reverse complements the sequences\n";
		help+="\t<input> = <fastaFile>\n";
		help+="reverseComplementGenome - reverse complements the lines and shifts the headers. The file needs to be piped through tac or reversed line-by-line in some other fashion\n";
		help+="\t<input> = <fastaFile>\n";
		help+="splitContigs - Takes a fasta-file. Splits the sequences into n bp big parts, with the given coverage\n";
		help+="\t<input> = <fastaFile> <size> <coverage>\n";
		help+="splitNewblerContigs - Takes a file with Newbler contigs and the 454ContigGraph.txt file... Splits the contigs into n bp big parts, with the given coverage if it is smaller than 100 (discards contigs with coverage higher than 100), but at least two if the contig is longer than n\n";
		help+="\t<input> = <fastaFile> <454ContigGrapth.txt> <size>\n";
		help+="splitOnN - Splits the sequences on all stretches of N longer or equeal to Length. Only prints sequences longer than minSeqLength\n";
		help+="\t<input> = <fastaFile> <length> <minSeqLength>\n";
		help+="splitVelvetContigs - Takes a file with velvet contigs... assumes the last part of the qname is the coverage. Splits the contigs into n bp big parts, with the given coverage if it is smaller than 100 (discards contigs with coverage higher than 100), but at least two if the contig is longer than n\n";
		help+="\t<input> = <fastaFile> <size>\n";
		help+="tails - prints n bp from each end and rev.comp the 3' end \n";
		help+="\t<input> = <fastaFile> <n>\n";
		help+="tailsTable - prints the representative n-mer bp from each end\n";
		help+="\t<input> = <fastaFile> <n>\n";
		help+="toFastq - combines a fasta and a qual file to a fastq file\n";
		help+="\tArguments: fastaUtils toFastq <fa file> <qual file>\n";
		help+="toFastq(Pipe) - combines a fasta and a qual file to a fastq file with quality 30 for all positions\n";
		help+="\tArguments: fastaUtils toFastq <fa file>\n";
		help+="toFastqRoche - combines a fasta and a Roche qual file to a fastq file\n";
		help+="\tArguments: fastaUtils toFastq <fa file> <Roche qual file>\n";
		help+="vistaPairsToPlot - takes all the pair-wise alignments for a reference in mVista and \n";
		help+="\t<input> = <pair fastaFile 1> ... <pair fastaFile n> \n";
		
		


		return help;
	}
	
	private static void pad(final String faFile, final int N)throws Exception{
		StringBuilder padding= new StringBuilder(256);
		while(padding.length()<N){
			padding.append('N');
		}
		final String Npad= padding.toString();
		
		fastaParser fp= new fastaParser(faFile);
		FastaSeq fs;
		
		for(;fp.hasNext();){
			fs=fp.next();
			System.out.println(fs.getHeader());
			System.out.println(Npad+fs.getSeq()+Npad);
		}
	}
	
	private static void kmerizeBis(String faFile, final int k)throws Exception{
		fastaKmerParser fkp= new fastaKmerParser(general.getBufferedReader(faFile), k);
		String kmer;
		
		for(;fkp.hasNext();){
			kmer=fkp.next();
			kmerizeBis_printCombinations(new StringBuilder(kmer), 0);
			kmerizeBis_printCombinations(new StringBuilder(sequenceUtils.reverseComplement(kmer)), 0);
		}
	}
	
	private static void kmerizeBis_printCombinations(StringBuilder kmer, int pos){
		if(pos==kmer.length()){
			System.out.println(kmer.toString());
		}else{
			kmerizeBis_printCombinations(kmer, pos+1);
			char curChar=kmer.charAt(pos);
			if(curChar=='C' || curChar=='c'){
				StringBuilder kmer2= new StringBuilder(kmer.toString());
				kmer2.setCharAt(pos, 'T');
				kmerizeBis_printCombinations(kmer2, pos);
			}
		}
	}
	
//	public static void cContext(String faFile, BufferedReader vcf,ArrayList<String> rawContexts)throws Exception{
//		final int contextSize=5,cPos=2;
//		multipatternTree pattern= new multipatternTree();
//		for (String rawContext : rawContexts) {
//			String[] l= rawContext.split(",");
//			String tag=l[0];
//			int pos= Integer.parseInt(l[1]);
//			while(pos<=cPos){
//				l[0]="N"+l[0];
//				++pos;
//			}
//			while(l[0].length()<contextSize){
//				l[0]=l[0]+"N";
//			}
//			pattern.consume(l[0], tag);
//		}
//		
//		RandomAccessFastaParser rafp= new RandomAccessFastaParser(faFile);
////		rafp.getAllNames();
//		FastaSeq fs= new FastaSeq(">NotAsequence", new StringBuffer());
//		BisSNPbedEntryParser bsbep= new BisSNPbedEntryParser(vcf);
//		BisSNPbedEntry bsbe;
//		String curKmer;
//		int count=0,snpPosCount=0;
//		int start,stop;
//		HashMap<String, Long> cCount= new HashMap<String, Long>(),
//				tCount= new HashMap<String, Long>(),
//				nrOfLoci= new HashMap<String, Long>();
//		
//		for(;bsbep.hasNext();++count){
//			bsbe=bsbep.next();
//			if(count%100000==0){
//				System.err.print(count+" "+bsbe.getChr()+"\r");
//			}
//			if(!bsbe.getChr().equals(fs.getQname())){
//				fs= rafp.get(bsbe.getChr());
//			}
//			if(bsbe.isPassed() && bsbe.getCoverage()>0){
//				start=bsbe.getStart()-cPos;
//				stop=start+contextSize;
//				if(start>0 && stop<fs.length()){
//					curKmer=fs.getSeq(start, stop);
//					if(curKmer.charAt(cPos)=='G'){
//						curKmer=sequenceUtils.reverseComplement(curKmer);
//					}
//					
//					if(curKmer.charAt(cPos)=='C'){
//						for(String context : pattern.get(curKmer)){
//							if(nrOfLoci.containsKey(context)){
//								nrOfLoci.put(context, nrOfLoci.get(context) + 1);
//								cCount.put(context, cCount.get(context) + bsbe.getcCount());
//								tCount.put(context, tCount.get(context) + bsbe.gettCount());
//							}else{
//								nrOfLoci.put(context, 1L);
//								cCount.put(context, (long) bsbe.getcCount());
//								tCount.put(context, (long) bsbe.gettCount());
//							}
//						}
//					}else{
//						//the reference doesn't contain a C at this position
//						++snpPosCount;
//					}
//				}
//			}
//		}
//		System.err.print(count+"\n");
//		System.err.println("Skipped "+snpPosCount+" positions as they weren't supported by the reference. They were SNPs and this is too simple of a calculation");
//		System.out.println("context\tavg(cov)\tavg(meth)\tnrOfLoci");
//		for(String context : nrOfLoci.keySet()){
//			double c=cCount.get(context),
//					t=tCount.get(context);
//			long nr= nrOfLoci.get(context);
//			System.out.println(context+sep+((c+t)/nr)+sep+(100*c/(c+t))+sep+nr);
//		}
//	}
	
	
	
	public static void extractPairedNoOrder(BufferedReader fa1, BufferedReader fa2, String outPrefix, int bufferSize)throws Exception{
		fastaParser fp1= new fastaParser(fa1),
				fp2= new fastaParser(fa2);
		
		BufferedWriter out1=new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outPrefix+"_1.fa.gz")))),
				out2=new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outPrefix+"_2.fa.gz")))),
				out1s=new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outPrefix+"singlet_1.fa.gz")))),
				out2s=new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outPrefix+"singlet_2.fa.gz"))));

		HashMap<String, FastaSeq> buffer1= new HashMap<String, FastaSeq>(bufferSize*5/4),
				buffer2= new HashMap<String, FastaSeq>(bufferSize*5/4);
		
		LinkedList<String> queue1= new LinkedList<String>(),
				queue2= new LinkedList<String>();
		
		FastaSeq fs1,fs2;
		String qname1,qname2;
		
		int count=0,paired=0,s1=0,s2=0;
		
		while(fp1.hasNext()||fp2.hasNext()){
			++count;
			if(count%10000==0){
				System.err.print(count+"###"+paired+"###"+s1+"###"+s2+"\r");
			}
			if(fp1.hasNext()){
				fs1=fp1.next();
				qname1=fs1.getQname();
				if(buffer2.containsKey(qname1)){
					out1.write(fs1.toString()+'\n');
					out2.write(buffer2.remove(qname1).toString()+'\n');
					++paired;
				}else{
					buffer1.put(qname1, fs1);
					queue1.add(qname1);
				}
			}
			
			if(fp2.hasNext()){
				fs2=fp2.next();
				qname2=fs2.getQname();
				if(buffer1.containsKey(qname2)){
					out2.write(fs2.toString()+'\n');
					out1.write(buffer1.remove(qname2).toString()+'\n');
					++paired;
				}else{
					buffer1.put(qname2, fs2);
					queue1.add(qname2);
				}
			}
			
			while(buffer1.size()>bufferSize){
				qname1=queue1.pop();
				if(buffer1.containsKey(qname1)){
					fs1=buffer1.remove(qname1);
					if(buffer2.containsKey(qname1)){
						out1.write(fs1.toString()+'\n');
						out2.write(buffer2.remove(qname1).toString()+'\n');
						++paired;
					}else{
						out1s.write(fs1.toString()+'\n');
						++s1;
					}
				}
			}
			
			while(buffer2.size()>bufferSize){
				qname2=queue2.pop();
				if(buffer2.containsKey(qname2)){
					fs2=buffer2.remove(qname2);
					if(buffer1.containsKey(qname2)){
						out2.write(fs2.toString()+'\n');
						out1.write(buffer1.remove(qname2).toString()+'\n');
						++paired;
					}else{
						out2s.write(fs2.toString()+'\n');
						++s2;
					}
				}
			}
		}
		while(buffer1.size()>0){
			qname1=queue1.pop();
			if(buffer1.containsKey(qname1)){
				fs1=buffer1.remove(qname1);
				if(buffer2.containsKey(qname1)){
					out1.write(fs1.toString()+'\n');
					out2.write(buffer2.remove(qname1).toString()+'\n');
					++paired;
				}else{
					out1s.write(fs1.toString()+'\n');
					++s1;
				}
			}
		}
		
		while(buffer2.size()>0){
			qname2=queue2.pop();
			if(buffer2.containsKey(qname2)){
				fs2=buffer2.remove(qname2);
				if(buffer1.containsKey(qname2)){
					out2.write(fs2.toString()+'\n');
					out1.write(buffer1.remove(qname2).toString()+'\n');
					++paired;
				}else{
					out2s.write(fs2.toString()+'\n');
					++s2;
				}
			}
		}
		System.err.println(count+"###"+paired+"###"+s1+"###"+s2);
		
		out1.close();
		out2.close();
		out1s.close();
		out2s.close();
		
		
	}
	
	public static void extractPaired(BufferedReader fa1, BufferedReader fa2, String outPrefix)throws Exception{
		fastaParser fp1= new fastaParser(fa1),
				fp2= new fastaParser(fa2);
		BufferedReader giOrder= new BufferedReader(new InputStreamReader(System.in));
		
		FastaSeq fs1=fp1.next(),
				fs2=fp2.next();
		
		String gi1=fs1.getQname()
				,gi2=fs2.getQname();
		
		BufferedWriter out1=new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outPrefix+"_1.fa.gz")))),
				out2=new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outPrefix+"_2.fa.gz")))),
				out1s=new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outPrefix+"singlet_1.fa.gz")))),
				out2s=new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outPrefix+"singlet_2.fa.gz"))));
		
		
//		out1=new BufferedWriter(new OutputStreamWriter(System.out));
		 
		int count=1,s1=0,s2=0,p=0;
		boolean keepOn=true;
		
		for(String gi=giOrder.readLine();gi!=null && keepOn;gi=giOrder.readLine(),++count){
			if(count%10000==0){
				System.err.print(count+"###"+p+"###"+s1+"###"+s2+"\r");
				System.err.println(gi1);
				System.err.println(gi2);
				break;
			}
			if(gi.equals(gi1)){
				if(gi.equals(gi2)){
					++p;
//					System.err.println("print");
					out1.write(fs1+"\n");
					out2.write(fs2+"\n");
					if(fp2.hasNext()){
						fs2=fp2.next();
						gi2=fs2.getQname();
					}else{
						keepOn=false;
					}
				}else{
					++s1;
					out1s.write(fs1+"\n");
				}
				if(fp1.hasNext()){
//					System.err.println("read1");
					fs1=fp1.next();
					gi1=fs1.getQname();
				}else{
					keepOn=false;
				}
			}else{
				if(gi.equals(gi2)){
//					System.err.println("read2");
					++s2;
					out2s.write(fs2+"\n");
					if(fp2.hasNext()){
						fs2=fp2.next();
						gi2=fs2.getQname();
					}else{
						keepOn=false;
					}
				}
			}
		}
		System.err.println(count+"###"+p+"###"+s1+"###"+s2+"\r");
		out1.close();
		out2.close();
		out1s.close();
		out2s.close();
	}
	
	public static void allAltPos(String fa, int pos)throws Exception{
		fastaParser fp= new fastaParser(new BufferedReader(new FileReader(fa)));
		FastaSeq fs;
		StringBuffer seq;
		String[] nucs= new String[]{"T","C","G"};
		for(;fp.hasNext();){
			fs=fp.next();
			seq= new StringBuffer(fs.getSeq(0,pos));
			seq.append("A");
			seq.append(fs.getSeq().substring(pos+1));
			System.out.println(">"+fs.getQname()+"_A");
			System.out.println(seq.toString());
			for(String nuc :nucs){
				seq.replace(pos, pos+1, nuc);
				System.out.println(">"+fs.getQname()+"_"+nuc);
				System.out.println(seq.toString());
			}
		}
		
	}
	
//	public static void kmerizeAllAltPosWithTag(String fa, int n,int pos)throws Exception{
//		fastaParser fp= new fastaParser(new BufferedReader(new FileReader(fa)));
//		FastaSeq fs;
//		StringBuffer seq;
//		String[] nucs= new String[]{"T","C","G"};
//		for(;fp.hasNext();){
//			fs=fp.next();
//			seq= new StringBuffer(fs.getSeq(0,pos));
//			seq.append("A");
//			seq.append(fs.getSeq().substring(pos+1));
////			System.err.println(seq.toString());
//			for(String k : kmerFunctions.kmerToUse(seq, n)){
//				System.out.println(k+sep+fs.getQname()+"_A");
//			}
//			for(String nuc :nucs){
//				seq.replace(pos, pos+1, nuc);
////				System.err.println(seq.toString());
//				for(String k : kmerFunctions.kmerToUse(seq, n)){
//					System.out.println(k+sep+fs.getQname()+"_"+nuc);
//				}
//			}
//		}
//	}
	
	public static void compare(String fa1,String fa2)throws Exception{
		BufferedReader in1= new BufferedReader(new FileReader(fa1));
		BufferedReader in2= new BufferedReader(new FileReader(fa2));
		String chrname="";
		int pos=1;
		char[] c1,c2;
		for(String s1=in1.readLine(),s2=in2.readLine();s1!=null&&s2!=null;s1=in1.readLine(),s2=in2.readLine()){
			if(s1.length()!=s2.length()){
				throw new Exception("files differently formatted");
			}
			if(s1.startsWith(">")){
				chrname=s1.split(" ")[0].substring(1);
				pos=1;
			}else{
				if(!s1.equals(s2)){
					c1=s1.toCharArray();
					c2=s2.toCharArray();
					for(int i=0;i<c1.length;++i){
						if(c1[i]!=c2[i]){
							System.out.println(chrname+sep+(pos+i)+sep+c1[i]+sep+c2[i]);
						}
					}
				}
				pos+=s1.length();
			}
		}
	}
	
//	public static void kmerCoverage_extractPerfectDirect(KmerSet_binary kmers,ArrayList<String> fastaFiles)throws Exception{
//		fastaParser fp;
//		FastaSeq fs;
//		BufferedWriter out;
//		int count;
//		System.err.println("parsing fastq files");
//		for (String fileName : fastaFiles) {
//			count=0;
//			out= new BufferedWriter(new FileWriter(fileName+".perfect.fa"));
//			fp=new fastaParser(new BufferedReader(new FileReader(fileName)));
//			for(;fp.hasNext();++count){
//				if(count%10000==0){
//					System.err.print(fileName+":  "+count+"\r");
//				}
//				fs=fp.next();
//				if(kmers.isPerfectND(fs.getSeq())){
//					fs.write(out);
//				}
//			}
//			System.err.println(fileName+":  "+count);
//			out.close();
//		}
//	}

	
	public static void minLength(BufferedReader in,int n)throws Exception{
		fastaParser fp=new fastaParser(in);
		FastaSeq fs;
		for(;fp.hasNext();){
			fs=fp.next();
			if(fs.length()>=n){
				System.out.println(fs);
			}
		}
	}
	
//	public static void kmerLinks(BufferedReader in, int kmerSize)throws Exception{
//		fastaParser fp=new fastaParser(in);
//		ArrayList<String>kmers;
//		String last;
//		int i=0;
//		for(;fp.hasNext();++i){
//			if(i%10000==0){
//				System.err.print("      "+i+"\r");
//			}
//			kmers=kmerFunctions.kmerToUse(fp.next().getSeq(), kmerSize);
//			last=kmers.remove(0);
//			for (String kmer : kmers) {
//				System.out.println(last+sep+kmer);
//				last=kmer;
//			}
//		}
//		System.err.println("      "+i);
//	}
	
//	public static void removeDoublets(BufferedReader in)throws Exception{
//		fastaParser fp=new fastaParser(in);
//		FastaSeq fs;
//		String rep;
//		HashSet<String> printed=new HashSet<String>();
//		int i=0;
//		for(;fp.hasNext();++i){
//			if(i%10000==0){
//				System.err.print("      "+i+"\r");
//			}
//			fs=fp.next();
//			rep=kmerFunctions.kmerToUse(fs.getSeq());
//			if(!printed.contains(rep)){
//				System.out.println(fs);
//				printed.add(rep);
//			}
//		}
//		System.err.println("      "+i);
//	}
	
//	public static void kmerStat(String faFile, String kmerFile, int kmerPos)throws Exception{
//		HashMap<String, int[]> kmers= new HashMap<String, int[]>(10000000);
//		BufferedReader in= new BufferedReader(new FileReader(kmerFile));
//		String[] l;
//		int[] counts,totCounts;
//		for(String s=in.readLine();s!=null;s=in.readLine()){
//			l=s.split("\t");
//			if (l.length>1) {
//				counts = new int[l.length - 1];
//				for (int i = 0, j = 0; i < l.length; ++i) {
//					if (i != kmerPos) {
//						counts[j] = Integer.parseInt(l[i]);
//						++j;
//					}
//				}
//				kmers.put(l[kmerPos], counts);
//			}
//		}
//		fastaParser fp= new fastaParser(new BufferedReader(new FileReader(faFile)));
//		FastaSeq fs;
//		StringBuffer seq;
//		final int kmerSize=kmers.keySet().iterator().next().length();
//		System.err.println("kmerSize: "+kmerSize);
//		final int countSize= kmers.values().iterator().next().length;
//		System.err.println("countSize: "+countSize);
//		for(;fp.hasNext();){
//			fs=fp.next();
//			totCounts=new int[countSize];
//			for(int i=0;i<countSize;++i){
//				totCounts[i]=0;
//			}
//			seq= new StringBuffer(fs.getSeq());
//			for(String kmer : kmerFunctions.kmerToUse(seq, kmerSize)){
//				if (kmers.containsKey(kmer)) {
//					counts = kmers.get(kmer);
//					for (int i = 0; i < countSize; ++i) {
//						totCounts[i] += counts[i];
//					}
//				}
//			}
//			System.out.print(fs.getQname()+sep+seq.length());
//			for(int i=0;i<countSize;++i){
//				System.out.print(sep+totCounts[i]);
//			}
//			System.out.println();
//		}
//	}
	
//	public static void cutAtBlastStart(String faFile, String blastFile)throws Exception{
//		blastM8Parser bp= new blastM8Parser(new BufferedReader(new FileReader(blastFile)));
//		blastM8Alignment ba;
//		HashMap<String, ArrayList<Integer>> cuts= new HashMap<String, ArrayList<Integer>>();
//		fastaParser fp= new fastaParser(new BufferedReader(new FileReader(faFile)));
//		FastaSeq fs;
//		int pos,start,end;
//		ArrayList<Integer> cutList;
//		
//		for(;bp.hasMore();){
//			ba= bp.nextHit();
//			pos=ba.strand[0]=='+'?ba.tstart:ba.tend;
//			if(cuts.containsKey(ba.tname)){
//				cuts.get(ba.tname).add(pos);
//			}else{
//				cutList= new ArrayList<Integer>();
//				cutList.add(pos);
//				cuts.put(ba.tname, cutList);
//			}
//		}
//		
//		for(;fp.hasNext();){
//			fs=fp.next();
//			if(cuts.containsKey(fs.getQname())){
//				cutList=cuts.get(fs.getQname());
//				Collections.sort(cutList);
//				start=1;
//				for (int i=0;i<cutList.size();++i){
//					end=cutList.get(i);
////					System.out.println(">"+fs.getQname()+"_"+i+" "+start+" "+end+"\n"+fs.getSeq(start, end));
//					System.out.println(fs.getQname()+"\t"+i+"\t"+start+"\t"+end);
//					start=end;
//				}
//				end=fs.length();
////				System.out.println(">"+fs.getQname()+"_"+cutList.size()+" "+start+" "+end+"\n"+fs.getSeq(start, end));
//				System.out.println(fs.getQname()+"\t"+cutList.size()+"\t"+start+"\t"+end);
//				
//			}else{
////				System.out.println(fs);
//				System.out.println(fs.getQname()+"\t0\t1\t"+fs.length());
//			}
//		}
//	}
	
	public static void splitOnN(String faFile,int length,int minSeqLength)throws Exception{
		fastaParser fp= new fastaParser(new BufferedReader(new FileReader(faFile)));
		FastaSeq fs;
		String[] seqs;
		for(;fp.hasNext();){
			fs=fp.next();
			seqs=fs.getSeq().split("N{"+length+",}");
			if(seqs.length>1){
				for(int i=1;i<=seqs.length;++i){
					if(seqs[i-1].length()>=minSeqLength){
						System.out.println(">"+fs.getQname()+"_"+i+"\n"+seqs[i-1]);
					}
				}
			}else{
				System.out.println(fs.toString());
			}
		}
	}
	
//	public static void mask(String faFile, String gffFile)throws Exception{
//		BufferedReader in= new BufferedReader(new FileReader(gffFile));
//		gffLine gl;
//		HashMap<String, ArrayList<gffLine>> maskStore= new HashMap<String, ArrayList<gffLine>>();
//		ArrayList<gffLine> masks;
//		Comparator<gffLine> c= new gffLine_startComparator();
//		fastaParser fp= new fastaParser(new BufferedReader(new FileReader(faFile)));
//		FastaSeq fs;
//		StringBuffer seq;
//		int lastend;
//		
//		for(String s=in.readLine();s!=null;s=in.readLine()){
//			if(s.length()>0&& s.charAt(0)!='#'){
//				gl=new gffLine(s);
//				if(!maskStore.containsKey(gl.getChr())){
//					maskStore.put(gl.getChr(), new ArrayList<gffLine>());
//				}
//				maskStore.get(gl.getChr()).add(gl);
//			}
//		}
//		
//
//		for(;fp.hasNext();){
//			fs=fp.next();
//			if(maskStore.containsKey(fs.getQname())){
//				masks=maskStore.get(fs.getQname());
//				Collections.sort(masks, c);
//				System.out.println(fs.getHeader());
//				lastend=0;
//				seq=new StringBuffer(fs.getSeq());
//				for(int i=0;i<masks.size();++i){
//					gl=masks.get(i);
//					System.out.print(seq.substring(lastend, gl.getStart()-1));
//					System.out.print(mask_NStretch(gl.getEnd()-gl.getStart()+1));
//					lastend=gl.getEnd();
//				}
//				System.out.println(seq.substring(lastend));
//			}else{
//				System.out.println(fs.toString());
//			}
//		}
//	}
	
	private static String mask_NStretch(int i){
		StringBuffer r= new StringBuffer(i);
		while(r.length()<i){
			r.append('N');
		}
		return r.toString();
	}
	
//	public static void getMismatch(String blastFile,String faFile1,String faFile2)throws Exception{
//		blastM8Parser bmp= new blastM8Parser(new BufferedReader(new FileReader(blastFile)));
//		blastM8Alignment bma;
//		ArrayList<blastM8Alignment> bmas= new ArrayList<blastM8Alignment>();
//		HashMap<String,String> first=new HashMap<String, String>(), second=new HashMap<String, String>();
//		fastaParser fp;
//		FastaSeq fs;
//		StringBuffer seq1,seq2;
//		@SuppressWarnings("unused")
//		int modPos;
//		
//		//Parse blastfile
//		for(;bmp.hasMore();){
//			bma=bmp.nextHit();
//			first.put(bma.getQname(), "");
//			second.put(bma.getTname(), "");
//			bmas.add(bma);
//		}
//		//parse first fastaFile
//		fp= new fastaParser(new BufferedReader(new FileReader(faFile1)));
//		for(;fp.hasNext();){
//			fs=fp.next();
//			if(first.containsKey(fs.getQname())){
//				first.put(fs.getQname(), fs.getSeq());
//			}
//		}
//		//parse second fastaFile
//		fp= new fastaParser(new BufferedReader(new FileReader(faFile2)));
//		for(;fp.hasNext();){
//			fs=fp.next();
//			if(second.containsKey(fs.getQname())){
//				second.put(fs.getQname(), fs.getSeq());
//			}
//		}
//		
//		//merge information and generate kmers
//		for (blastM8Alignment bmatmp : bmas) {
//			seq1=new StringBuffer(first.get(bmatmp.getQname()).substring(bmatmp.qstart-1, bmatmp.qend));
//			if(bmatmp.strand[0]=='-'){
//				seq2=new StringBuffer(sequenceUtils.reverseComplement(second.get(bmatmp.getTname()).substring(bmatmp.tstart-1, bmatmp.tend)));
//			}else{
//				seq2=new StringBuffer(second.get(bmatmp.getTname()).substring(bmatmp.tstart-1, bmatmp.tend));
//			}
//			//identify changed position
//			modPos=-1;
////			System.err.println();
////			System.err.println(seq1);
//			for(int i=0;i<seq1.length();++i){
//				if(seq1.charAt(i)==seq2.charAt(i)){
////					System.err.print(" ");
//				}else{
//					modPos=i;
//					System.out.println(bmatmp.toString()+sep+seq1.charAt(i)+">"+seq2.charAt(i));
////					System.err.println("|");
//					break;
//				}
//			}
////			System.err.println(seq2);
//			
//		}
//	}
	
//	public static void kmerizeAllAlt(String blastFile,String faFile1,String faFile2,int kmerSize)throws Exception{
//		blastM8Parser bmp= new blastM8Parser(new BufferedReader(new FileReader(blastFile)));
//		blastM8Alignment bma;
//		ArrayList<blastM8Alignment> bmas= new ArrayList<blastM8Alignment>();
//		HashMap<String,String> first=new HashMap<String, String>(), second=new HashMap<String, String>();
//		fastaParser fp;
//		FastaSeq fs;
//		StringBuffer seq1,seq2;
//		int modPos;
//		
//		//Parse blastfile
//		for(;bmp.hasMore();){
//			bma=bmp.nextHit();
//			first.put(bma.getQname(), "");
//			second.put(bma.getTname(), "");
//			bmas.add(bma);
//		}
//		//parse first fastaFile
//		fp= new fastaParser(new BufferedReader(new FileReader(faFile1)));
//		for(;fp.hasNext();){
//			fs=fp.next();
//			if(first.containsKey(fs.getQname())){
//				first.put(fs.getQname(), fs.getSeq());
//			}
//		}
//		//parse second fastaFile
//		fp= new fastaParser(new BufferedReader(new FileReader(faFile2)));
//		for(;fp.hasNext();){
//			fs=fp.next();
//			if(second.containsKey(fs.getQname())){
//				second.put(fs.getQname(), fs.getSeq());
//			}
//		}
//		
//		//merge information and generate kmers
//		for (blastM8Alignment bmatmp : bmas) {
//			seq1=new StringBuffer(first.get(bmatmp.getQname()).substring(bmatmp.qstart-1, bmatmp.qend));
//			if(bmatmp.strand[0]=='-'){
//				seq2=new StringBuffer(sequenceUtils.reverseComplement(second.get(bmatmp.getTname()).substring(bmatmp.tstart-1, bmatmp.tend)));
//			}else{
//				seq2=new StringBuffer(second.get(bmatmp.getTname()).substring(bmatmp.tstart-1, bmatmp.tend));
//			}
//			//identify changed position
//			modPos=-1;
////			System.err.println();
////			System.err.println(seq1);
//			for(int i=0;i<seq1.length();++i){
//				if(seq1.charAt(i)==seq2.charAt(i)){
////					System.err.print(" ");
//				}else{
//					modPos=i;
////					System.err.println("|");
//					break;
//				}
//			}
////			System.err.println(seq2);
//			for (String nuc : new String[]{"A","T","G","C"}) {
//				seq1.replace(modPos, modPos+1, nuc);
////				System.err.println(seq1);
//				for (String kmer : kmerFunctions.kmerToUse(seq1, kmerSize)) {
//					System.out.println(kmer);
//				}
//			}
//		}
//	}
	
//	public static void reverseComplementFile(String faFile)throws Exception{
//		fastaParser fp= new fastaParser(new BufferedReader(new FileReader(faFile)));
//		FastaSeq fs;
//		int count;
//		for(count=0;fp.hasNext();count++){
//			if(count%1000==0){
//				System.err.print("          "+count+"\r");
//			}
//			fs=fp.next();
//			System.out.println(reverseComplement(fs));
//		}
//		System.err.println("          "+count);
//			
//	}
	
	public static void reverseComplementFile(String faFile)throws Exception{
		BufferedReader in= new BufferedReader(new FileReader(faFile));
		StringBuffer curSeq= new StringBuffer("");
		long count=0;
		int seq=0;
		for(String s= in.readLine();s!=null;s= in.readLine()){
			++count;
			if(count%100==0){
				System.err.print("lines read: "+count+"\r");
			}
			if(s.startsWith(">")){
				if(seq>0){
					System.out.print(curSeq);
					curSeq.setLength(0);
					System.out.print("\n"+s);
				}else{
					System.out.print(s);
				}
				++seq;
			}else{
				curSeq.insert(0, sequenceUtils.reverseComplement(s));
				curSeq.insert(0, '\n');
			}
		}
		in.close();
		System.err.println("lines read: "+count);
		System.out.println(curSeq);
	}
	
	public static void reverseComplementGenome(String genomeFile)throws Exception{
		BufferedReader in= new BufferedReader(new FileReader(genomeFile));
		String curHeader="";
		int seq=0;
		long count=0;
		for(String s= in.readLine();s!=null;s= in.readLine()){
			++count;
			if(count%10000==0){
				System.err.print("lines read: "+count+"\r");
			}
			if(s.startsWith(">")){
				if(seq>0){
					System.out.println(curHeader);
				}
				++seq;
				curHeader=s;
			}else{
				System.out.println(sequenceUtils.reverseComplement(s));
			}
		}
		in.close();
		System.err.println("lines read: "+count);
		System.out.println(curHeader);
	}
	
//	public static void kmerContains(String faFile,String kmerFile, final int kmerPos) throws Exception{
//		BufferedReader in= new BufferedReader(new FileReader(kmerFile));
//		String s=in.readLine();
//		String[] l=s.split("\t");
//		final int kmerLength=l[kmerPos].length();
//		KmerSet_binary kmers= new KmerSet_binary(kmerLength,false);
//		kmers.addKmer(l[kmerPos]);
//		for(s=in.readLine();s!=null;s=in.readLine()){
//			l=s.split("\t");
//			if(l.length>kmerPos){
//				kmers.addKmer(l[kmerPos]);
//			}
//		}
//		fastaParser fp= new fastaParser(new BufferedReader(new FileReader(faFile)));
//		FastaSeq fs;
//		StringBuffer kmer;
//		for(;fp.hasNext();){
//			fs=fp.next();
//			if(fs.getSeq().length()>=kmerLength){
//				kmer=new StringBuffer(" "+fs.getSeq(0, kmerLength-1));
//				final char[] seq=fs.getSeq().toCharArray();
//				for(int i=kmerLength-1;i<seq.length;++i){
//					kmer.deleteCharAt(0);
//					kmer.append(seq[i]);
//					if(kmers.exists(kmerFunctions.kmerToUse(kmer))){
//						System.out.println(fs.getQname());
//						break;
//					}
//				}
//			}
//
//		}
//	}
	
//	public static void kmer_frequency(String faFile, String kmerDistFile) throws Exception{
//		BufferedReader in= new BufferedReader(new FileReader(kmerDistFile));
//		String s=in.readLine();
//		String[] l=s.split("\t");
//		int count=0;
//		final int kmerLength=l[1].length();
//		System.err.println("\n\nReading the kmer distribution....");
//		KmerMap_binary kmers= new KmerMap_binary(kmerLength);
//		kmers.addKmer(l[1], Integer.parseInt(l[0]));
//		for(s=in.readLine();s!=null;s=in.readLine(),count++){
//			if(count%100000==0){
//				System.err.print("          "+count+"\r");
//			}
//			l= s.split("\t");
//			if(l.length>1){
//				kmers.addKmer(l[1], Integer.parseInt(l[0]));
//			}	
//		}
//		System.err.println("          "+count);
//		System.err.println("Size: "+kmers.size());
//		in.close();
//		
//		fastaParser fp= new fastaParser(new BufferedReader(new FileReader(faFile)));
//		FastaSeq fs;
//		StringBuffer kmer;
//		System.err.println("\n\nProcessing the sequences..."); 
//		for(count=0;fp.hasNext();count++){
//			if(count%10==0){
//				System.err.print("          "+count+"\r");
//			}
//			fs=fp.next();
//			System.out.print(fs.getQname());
//			if(fs.getSeq().length()>=kmerLength){
//				kmer=new StringBuffer(" "+fs.getSeq(0, kmerLength-1));
//				final char[] seq=fs.getSeq().toCharArray();
//				for(int i=kmerLength-1;i<seq.length;++i){
//					kmer.deleteCharAt(0);
//					kmer.append(seq[i]);
//					System.out.print(sep+kmers.get(kmer)); //this function already converts to the right kmer
//				}
//			}
//			System.out.println();
//		}
//		System.err.println("          "+count);
//	}
	
	public static void tails(String faFile, int n)throws Exception{
		fastaParser fp= new fastaParser(new BufferedReader(new FileReader(faFile)));
		FastaSeq fs;
		
		for(;fp.hasNext();){
			fs=fp.next();
			System.out.println(">"+fs.getQname()+"_5p "+n);
			System.out.println(fs.getSeq(0, n));
			System.out.println(">"+fs.getQname()+"_3p "+n+" rev.comp.");
			System.out.println(sequenceUtils.reverseComplement(fs.getSeq(fs.length()-n, fs.length())));
		}
	}
	
//	public static void tailsTable(String faFile, int n)throws Exception{
//		fastaParser fp= new fastaParser(new BufferedReader(new FileReader(faFile)));
//		FastaSeq fs;
//		
//		for(;fp.hasNext();){
//			fs=fp.next();
//			System.out.println(kmerFunctions.kmerToUse(fs.getSeq(0, n))+sep+fs.getQname()+sep+"1");
//			System.out.println(kmerFunctions.kmerToUse(fs.getSeq(fs.length()-n, fs.length()))+sep+fs.getQname()+sep+"2");
//		}
//	}
	
	public static void kmerizeWRevComp(String faFile,int length)throws Exception{
		fastaKmerParser fkp=new fastaKmerParser(new BufferedReader(new FileReader(faFile)), length);
		String kmer;
		for(;fkp.hasNext();){
			kmer=fkp.next();
			System.out.println(kmer);
			System.out.println(sequenceUtils.reverseComplement(kmer));
		}
	}
	
	public static void kmerizeNoRC(BufferedReader faFile,int length)throws Exception{
		fastaKmerParser fkp=new fastaKmerParser(faFile, length);
		for(;fkp.hasNext();){
			System.out.println(fkp.next());
		}
	}
	
//	public static void kmerize(BufferedReader faFile,int length)throws Exception{
//		fastaKmerParser fkp=new fastaKmerParser(faFile, length);
//		for(;fkp.hasNext();){
//			System.out.println(kmerFunctions.kmerToUse(fkp.next()));
//		}
//	}
	
//	public static void kmerizeWithTag(String faFile, int length) throws Exception{
//		fastaParser fp= new fastaParser(new BufferedReader(new FileReader(faFile)));
//		FastaSeq fs;
//		String kmer;
//		int cutoff;
//		for (;fp.hasNext();){
//			fs=fp.next();
//			if(fs.length()>length){
//				cutoff=(fs.length()+length-1)/2;
//				kmer=fs.getSeq(0, length);
//				System.out.println(kmerFunctions.kmerToUse(kmer)+sep+fs.getQname()+sep+"L");
//				for(int i=length;i<fs.length();++i){
//					kmer=kmer.substring(1)+fs.getSeq().charAt(i);
//					if(i<=cutoff){
//						System.out.println(kmerFunctions.kmerToUse(kmer)+sep+fs.getQname()+sep+"L");
//					}else{
//						System.out.println(kmerFunctions.kmerToUse(kmer)+sep+fs.getQname()+sep+"R");
//					}
//				}
//			}
//		}
//	}
	
//	public static void kmerizeWithTagNoRC(String faFile, int length) throws Exception{
//		fastaParser fp= new fastaParser(new BufferedReader(new FileReader(faFile)));
//		FastaSeq fs;
//		String kmer;
//		int cutoff;
//		for (;fp.hasNext();){
//			fs=fp.next();
//			if(fs.length()>length){
//				cutoff=(fs.length()+length-1)/2;
//				kmer=fs.getSeq(0, length);
//				System.out.println(kmerFunctions.kmerToUse(kmer)+sep+fs.getQname()+sep+"L");
//				for(int i=length;i<fs.length();++i){
//					kmer=kmer.substring(1)+fs.getSeq().charAt(i);
//					if(i<=cutoff){
//						System.out.println(kmer+sep+fs.getQname()+sep+"L");
//					}else{
//						System.out.println(kmer+sep+fs.getQname()+sep+"R");
//					}
//				}
//			}
//		}
//	}

	public static void count(String faFile,char X)throws Exception{
		fastaParser fp= new fastaParser(new BufferedReader(new FileReader(faFile)));
		FastaSeq fs;
		int count;
		for(;fp.hasNext();){
			fs=fp.next();
			count=0;
			for (char N : fs.getSeq().toCharArray()) {
				if(N==X){
					count++;
				}
			}
			System.out.println(fs.getQname()+sep+count);
		}
	}
	
	public static void nCount(String faFile)throws Exception{
		fastaParser fp= new fastaParser(new BufferedReader(new FileReader(faFile)));
		FastaSeq fs;
		int count,pos;
		System.out.println("qname"+sep+"#N"+sep+"end"+sep+"length");
		for(;fp.hasNext();){
			fs=fp.next();
			count=0;
			pos=0;
			for (char c : fs.getSeq().toCharArray()) {
				pos++;
				if(c=='N'){
					count++;
				}else if(count>0){
					System.out.println(fs.getQname()+sep+count+sep+pos+sep+fs.length());
					count=0;
				}
			}
			if(count>0){
				System.out.println(fs.getQname()+sep+count+sep+pos+sep+fs.length());
			}
		}
	}
	
	public static void vistaPairsToPlot(ArrayList<String> pairFiles) throws Exception{
		ArrayList<ArrayList<FastaSeq>> pairs= new ArrayList<ArrayList<FastaSeq>>();
		fastaParser fp;
		ArrayList<FastaSeq> cur;
		for (String pairFile : pairFiles) {
			cur= new ArrayList<FastaSeq>();
			fp= new fastaParser(new BufferedReader(new FileReader(pairFile)));
			cur.add(fp.next());
			cur.add(fp.next());
			pairs.add(cur);
		}
		int[] pairPos=new int[pairs.size()];
		for(int i=0;i<pairs.size();i++){
			pairPos[i]=-1;
		}
		Hashtable<String, Integer> majorityCount;
		vistaPairsToPlot_pos pos;
		@SuppressWarnings("unused")
		String ref,maxnuc;
		int max;
		System.out.println("Pos"+sep+"Ref"+sep+"count");
		for(int curPos=1;;curPos++){
			pos=vistaPairsToPlot_nextPos(pairs.get(0), pairPos[0]);
			if(pos.getNewPos()==-616){
				//flush and finish
				
				break;
			}else{
				ref=pos.getOne();
				maxnuc=ref;
				majorityCount= new Hashtable<String, Integer>();
				if(ref.equals(pos.getTwo())){
					majorityCount.put(ref, 2);
					max=2;
				}else{
					majorityCount.put(ref, 1);
					max=1;
					majorityCount.put(pos.getTwo(), 1);
				}
				pairPos[0]=pos.getNewPos();
				for(int i=1;i<pairs.size();i++){
					pos=vistaPairsToPlot_nextPos(pairs.get(i), pairPos[i]);
					if(!ref.equals(pos.getOne())){
						System.err.println("out of phase");
					}
					if(majorityCount.containsKey(pos.getTwo())){
						if(majorityCount.get(pos.getTwo())==max){
							max++;
							maxnuc=pos.getTwo();
						}
						majorityCount.put(pos.getTwo(), majorityCount.get(pos.getTwo())+1);
					}else{
						majorityCount.put(pos.getTwo(), 1);
					}
					pairPos[i]=pos.getNewPos();
				}
				System.out.println(curPos+sep+ref+sep+majorityCount.get(ref));
			}
		}
	}
	
	public static vistaPairsToPlot_pos vistaPairsToPlot_nextPos(ArrayList<FastaSeq> seqs, int startPos){
		int pos=startPos+1;
		final String ref=seqs.get(0).getSeq();
		for(;pos<ref.length();pos++){
			if(nucleotides.contains(ref.charAt(pos)+"")){
				break;
			}else if (ref.charAt(pos)!='-'){
				System.err.print(ref.charAt(pos)+"");
			}
		}
		if(pos<ref.length()){
			return new vistaPairsToPlot_pos(ref.charAt(pos)+"", seqs.get(1).getSeq().charAt(pos)+"", pos);
		}else{
			return new vistaPairsToPlot_pos("", "", -616);
		}
	}
	
	public static void splitContigs(String faFile, int outSize, int coverage) throws Exception{
		fastaParser fp= new fastaParser(new BufferedReader(new FileReader(faFile)));
		FastaSeq fs;
		final int maxGap=2;
		String[] seqs;
		for(;fp.hasNext();){
			fs=fp.next();
			seqs=fs.getSeq().split("N{"+maxGap+",}");
			for(int i=0;i<seqs.length;i++){
				splitContigs_sub(new FastaSeq(">"+fs.getQname()+"_"+i,seqs[i]), coverage, outSize);
			}
		}
	}
	
	public static void splitVelvetContigs(String faFile, int outSize) throws Exception{
		fastaParser fp= new fastaParser(new BufferedReader(new FileReader(faFile)));
		FastaSeq fs;
		double coverage;
		int printcoverage;
		final int maxGap=2;
		String[] seqs;
		for(;fp.hasNext();){
			fs=fp.next();
			String[] l=fs.getQname().split("_");
			coverage=Double.parseDouble(l[l.length-1]);
			printcoverage=(int)Math.round(coverage);
			printcoverage=printcoverage<1?1:printcoverage;
//			printcoverage=printcoverage>10?10:printcoverage;
//			splitContigs_sub(fs, printcoverage, outSize);
			if(printcoverage<100){
				seqs=fs.getSeq().split("N{"+maxGap+",}");
				for(int i=0;i<seqs.length;i++){
					splitContigs_sub(new FastaSeq(">"+fs.getQname()+"_"+i,seqs[i]), printcoverage, outSize);
				}
			}
		}
	}
	
	public static void splitNewblerContigs(String faFile,String graphFile, int outSize) throws Exception{
		//parse the graph file for coverage
		HashMap<String, Integer> coverages= new HashMap<String, Integer>();
		BufferedReader in= new BufferedReader(new FileReader(graphFile));
		String[] l;
		for(String s=in.readLine();s!=null;s=in.readLine()){
			l=s.split("\t");
			if(l.length>3){
				if(l[1].length()>6 &&l[1].substring(0, 6).equals("contig")){
					coverages.put(l[1], (int)Math.round(Double.parseDouble(l[3])));
				}else{
					break;
				}
			}
		}
		in.close();
		//parse fastafile and split
		// HOW ABOUT THE QUAL FILE??? SHOULD IT ALSO BE SPLIT??
		fastaParser fp= new fastaParser(new BufferedReader(new FileReader(faFile)));
		FastaSeq fs;
		int coverage;
		for(;fp.hasNext();){
			fs=fp.next();
			coverage=coverages.get(fs.getQname());
			if(coverage<100){
				splitContigs_sub(fs, coverage, outSize);
			}
		}
	}
	
	private static void splitContigs_sub(FastaSeq fs, int coverage, int outSize) throws Exception{
		final String qbase=fs.getQname();
		final int minLength=75;
		int printStart,printStop;
		String seq;
		if(fs.length()>outSize){
			coverage=coverage==1?2:coverage;
			for(int start=-(coverage-1)*outSize/coverage,i=0;start<fs.length();start+=outSize/coverage,i++){
				printStart=start>0?start:0;
				printStop=start+outSize<fs.length()?start+outSize:fs.length();
//				if(printStop-printStart>minLength){
//					System.out.println(">"+qbase+"_"+i+"\n"+fs.getSeq(printStart, printStop));
//				}
				seq=splitContigs_trimN(fs.getSeq(printStart, printStop));
				if(seq.length()>minLength){
					System.out.println(">"+qbase+"_"+i+"\n"+seq);
				}
			}
		}else{
			//just print the sequence 'coverage' times
//			if(fs.length()>minLength){
//				System.out.println(fs.toString());
//			}
			seq=splitContigs_trimN(fs.getSeq());
			if(seq.length()>minLength){
				for(int i=0;i<coverage;i++){
					System.out.println(">"+qbase+"_"+i+"\n"+seq);
				}
			}
		}
	}
	
	private static String splitContigs_trimN(String seq){
		final int count=seq.length();
		int len=seq.length();
		int st=0;
		
		final char[] val =seq.toCharArray();
		while((st<len) && val[st]=='N'){
			st++;
		}
		while((st<len)&&val[len-1]=='N'){
			len--;
		}
		return((st>0)||(len<count))?seq.substring(st,len):seq;
	}
	

//	private static boolean splitContigs_nCheck(String seq,int minLength){
//		int count=0,nCount=0,min=Integer.MAX_VALUE;
//		for(int i=0;i<seq.length();i++){
//			if(seq.charAt(i)=='N'){
//				min=min<count?min:count>0?count:min;
//				count=0;
//				nCount++;
//			}else{
//				count++;
//			}
//		}
//		min=min<count?min:count;
//		return min>minLength && (double)nCount/seq.length()<0.6;
//	}
	
	public static void MSAtoVistaDiagram(String msaFile,String outfile, int windowSize, boolean onlyRef, String gapSymbols,boolean majority)throws Exception{
		HashSet<String> gaps= new HashSet<String>();
		for (char gap : gapSymbols.toCharArray()) {
			gaps.add(gap+"");
		}
		fastaParser fp= new fastaParser(new BufferedReader(new FileReader(msaFile)));
		FastaSeq fs;
		ArrayList<FastaSeq> seqs= new ArrayList<FastaSeq>();
		int length=0;
		//Check reference
		if(fp.hasNext()){
			fs=fp.next();
			length=fs.length();
			seqs.add(fs);
		}
		for(;fp.hasNext();){
			fs=fp.next();
			if(fs.length()!=length){
				throw new Exception("10 The sequences in the input file has different lengths and is not a MSA. The reference is "+length+" bp, and sequence: "+fs.getQname()+" is "+fs.length()+" bp.");
			}
			seqs.add(fs);
		}
		final int nrOfSeqs=seqs.size();
		if(onlyRef){
			//clean sequences
			final String ref=seqs.get(0).getSeq();
			boolean[] bad=new boolean[length];
			for(int i=0;i<length;i++){
				bad[i]=gaps.contains(ref.charAt(i)+"");
			}
			for (FastaSeq fastaSeq : seqs) {
				String newSeq="";
				final String oldSeq=fastaSeq.getSeq();
				for(int i=0;i<length;i++){
					if(!bad[i]){
						newSeq+=oldSeq.charAt(i);
					}
				}
				fastaSeq.setSeq(newSeq);
			}
			length=seqs.get(0).getSeq().length();
		}
		//prep
		final char[][] charSeq=new char[nrOfSeqs][length];
		String gapLine="";
		for(int i=0;i<seqs.size();i++){
			charSeq[i]=seqs.get(i).getSeq().toCharArray();
			gapLine+="-";
		}
		//count
		//Should be altered if this should work on really large alignments
		int[] counts= new int[length];
		String[] lines= new String[length];
		int count,sum=0,curCount;
		final int shift=windowSize/2;
		String outString,line;
		BufferedWriter out= new BufferedWriter(new FileWriter(outfile));
		HashMap<String, Integer> majorityCount;
		for(int i=0;i<length;i++){
			count=1;
			line=""+charSeq[0][i];
			majorityCount= new HashMap<String, Integer>();
			majorityCount.put(""+charSeq[0][i], 1);
			for(int j=1;j<nrOfSeqs;j++){
				final String key=""+charSeq[j][i];
				line+=key;
				if(majorityCount.containsKey(key)){
					curCount=majorityCount.get(key)+1;
					if(curCount>count){
						count=curCount;
					}
					majorityCount.put(key,curCount);
				}else{
					majorityCount.put(key, 1);
				}
			}
			if(!majority){
				count=majorityCount.get(""+charSeq[0][i]);
			}
			lines[i]=line;
			counts[i]=count;
			sum+=count;
			if(i>=windowSize){
				//subtract
				sum-=counts[i-windowSize];
			}
			if(i-shift>=0){
				line=lines[i-shift];
			}else{
				line=gapLine;
			}
			outString=(i-shift+1)+"\t"+line+"\t"+((double)sum)/(nrOfSeqs)/windowSize;
			out.write(outString+"\n");
			System.out.println(outString);
		}
		out.close();
	}
	
	public static void randomCut(String faFile,int size)throws Exception{
		fastaParser fp= new fastaParser(new BufferedReader(new FileReader(faFile)));
		FastaSeq fs;
		Random randomGen= new Random();
		int start;
		for(;fp.hasNext();){
			fs=fp.next();
			if(fs.length()>=size){
				start=randomGen.nextInt(fs.length()-size+1);
				System.out.println(">"+fs.getQname()+" "+start+" "+(start+size)+"\n"+fs.getSeq(start, start+size));
			}
		}
	}
	
	public static void toFastqRoche(String faFile, String qualFile)throws Exception{
		fastaParser fp= new fastaParser(new BufferedReader(new FileReader(faFile)));
		FastaSeq fs,qs;
		RocheQualParser rqp= new RocheQualParser(new BufferedReader(new FileReader(qualFile)));
		RocheQualSeq rqs;
		FastqSeq fqs;
		//assumes the sequences are ordered in the same way
		for(;fp.hasNext()&&rqp.hasNext();){
			fs=fp.next();
			rqs=rqp.next();
			if(!fs.getQname().equals(rqs.getQname())){
				throw new Exception("Unmatched files: "+fs.getHeader()+" is not equal to "+rqs.getHeader());
			}
			qs=rqs.toSangerQual();
			if(fs.length()!=qs.length()){
				throw new Exception("The fasta and quality sequences are of different lengths for "+fs.getQname()+" ("+fs.length()+" bp and "+qs.length()+" bp, respectively)\n"+fs+"\n"+qs+"\n"+rqs);
			}
			fqs= new FastqSeq("@"+fs.getHeader().substring(1), fs.getSeq(), qs.getSeq());
			System.out.println(fqs);
		}
	}
	
	public static void toFastq(BufferedReader faFile, String qualFile)throws Exception{
		fastaParser fp= new fastaParser(faFile);
		FastaSeq fs;
		fastaParser qp = new fastaParser(new BufferedReader(new FileReader(qualFile)));
		FastaSeq qs;
		String header;
		//assumes the sequences are ordered in the same way
		for(;fp.hasNext()&&qp.hasNext();){
			fs=fp.next();
			qs=qp.next();
			if(!fs.getQname().equals(qs.getQname())){
				throw new Exception("Unmatched files: "+fs.getHeader()+" is not equal to "+qs.getHeader());
			}
			header= fs.getHeader().substring(1);
			if(fs.length()!=qs.length()){
				throw new Exception("The fasta and quality sequences are of different lengths for "+fs.getQname()+" ("+fs.length()+" bp and "+qs.length()+" bp, respectively)");
			}
			System.out.println("@"+header);
			System.out.println(fs.getSeq());
			System.out.println("+");
			System.out.println(qs.getSeq());
		}
	}

	public static void toFastq(BufferedReader faFile)throws Exception {
		fastaParser fp= new fastaParser(faFile);
		for(;fp.hasNext();){
			System.out.println((new FastqSeq(fp.next())).toString());
		}
	}
	
	public static void motifFindReg(String faFile, String regExp)throws Exception{
		fastaParser fp;
		FastaSeq fs;
		Pattern motif= Pattern.compile(regExp);
		
		if(faFile.equals("-")){
			fp=new fastaParser(new BufferedReader(new InputStreamReader(System.in)));
		}else{
			fp= new fastaParser(new BufferedReader(new FileReader(faFile)));
		}
		for(int nr=0;fp.hasNext();nr++){
			if(nr%1000==0){
				System.err.println(nr);
			}
			fs=fp.next();
			Matcher seq= motif.matcher(fs.getSeq());
			final String qname=fs.getQname();
			while (seq.find()){
				System.out.println(qname+sep+seq.start()+sep+seq.end());
			}
		}
		fp.close();
	}
	
	
	public static void motifFind(String faFile, String motifFile)throws Exception{
		fastaParser fp;
		
		
		fp= new fastaParser(new BufferedReader(new FileReader(motifFile)));
		FastaSeq fs,rev;
		ArrayList<FastaSeq> forward= new ArrayList<FastaSeq>();
		ArrayList<FastaSeq> reverse= new ArrayList<FastaSeq>();
		for(;fp.hasNext();){
			//Gather motifs
			fs=fp.next();
			forward.add(fs);
			rev=reverseComplement(fs);
			if(!fs.getSeq().equals(rev.getSeq())){
				reverse.add(rev);
			}
		}
		
		if(faFile.equals("-")){
			fp=new fastaParser(new BufferedReader(new InputStreamReader(System.in)));
		}else{
			fp= new fastaParser(new BufferedReader(new FileReader(faFile)));
		}
		for(int nr=0;fp.hasNext();nr++){
			if(nr%1000==0){
				System.err.println(nr);
			}
			fs=fp.next();
			final String seq=fs.getSeq();
			final String qname=fs.getQname();
			//Search forward motifs
			for (FastaSeq motif : forward) {
				final String motifSeq=motif.getSeq();
				final String motifQname=motif.getQname();
				for (int lastStart=0;lastStart<seq.length();){
					int start=seq.indexOf(motifSeq, lastStart);
					if(start==-1){
						lastStart=seq.length();
					}else{
						System.out.println(qname+sep+motifQname+sep+"+"+sep+start+sep+(start+motifSeq.length()));
						lastStart=start+1;
					}
				}
			}
			//Search reverse motifs
			for (FastaSeq motif : reverse) {
				final String motifSeq=motif.getSeq();
				final String motifQname=motif.getQname();
				for (int lastStart=0;lastStart<seq.length();){
					int start=seq.indexOf(motifSeq, lastStart);
					if(start==-1){
						lastStart=seq.length();
					}else{
						System.out.println(qname+sep+motifQname+sep+"-"+sep+start+sep+(start+motifSeq.length()));
						lastStart=start+1;
					}
				}
			}
		}
	}
	
//	private static String reverseComplement(String in)throws Exception{
//		String seqOut="";
//		for (char c : in.toCharArray()) {
//			seqOut=complement(c)+seqOut;
//		}
//		return seqOut;
//	}
	
	private static FastaSeq reverseComplement(FastaSeq in)throws Exception{
		return new FastaSeq(in.getHeader(), sequenceUtils.reverseComplement(in.getSeq()));
	}
	
	public static void GCcontentPerRow(String faFile, final int rowLength)throws Exception{
		BufferedReader in= new BufferedReader(new FileReader(faFile));
		int sum;
		String s=in.readLine();
		while(!s.startsWith(">")){
			s=in.readLine();
		}
		System.out.println("fixedStep chrom="+s.substring(1).split(" ")[0]+" start=0 step="+rowLength);
		for(s=in.readLine();s!=null;s=in.readLine()){
			if(s.length()>0){
				sum=0;
				if(s.length()==rowLength){
					for(int i=0;i<rowLength;i++){
						final char cur=s.charAt(i);
						if(cur=='G'||cur=='C'){
							sum++;
						}
					}
				}else{
					System.err.println("Row of length: "+s.length()+"\n"+s+"\n");
					for(int i=0;i<s.length();i++){
						final char cur=s.charAt(i);
						if(cur=='G'||cur=='C'){
							sum++;
						}
					}
				}
				System.out.println(sum);
			}
		}
	}
	
	public static void culmativeXcount(String faFile) throws Exception{
		fastaParser fp= new fastaParser(new BufferedReader(new FileReader(faFile)));
		FastaSeq fs;
		int row=1;
		int sum=0;
		for (;fp.hasNext();){
			fs=fp.next();
			if(fs.getSeq().indexOf('X')==-1){
				System.out.println((row++)+sep+fs.getQname()+sep+"0"+sep+sum);
			}else{
				System.out.println((row++)+sep+fs.getQname()+sep+"1"+sep+(++sum));
			}
			
			
		}
	}
	
	public static void countPos(String faFile) throws Exception{
		fastaParser fp= new fastaParser(new BufferedReader(new FileReader(faFile)));
		ArrayList<String> alphabet= new ArrayList<String>();
		int maxLength=0;
		ArrayList< HashMap<String, Integer>> count= new ArrayList< HashMap<String,Integer>>();
		FastaSeq fs;
		for (;fp.hasNext();){
			fs=fp.next();
			final String seq= fs.getSeq();
			for(int i=0;i<seq.length();i++){
				final char cur= seq.charAt(i);
				//very redundant solution
				if(!alphabet.contains(cur+"")){
					alphabet.add(cur+"");
				}
				if(i==maxLength){
					count.add( new HashMap<String, Integer>());
					maxLength++;
				}
				if(count.get(i).containsKey(cur+"")){
					count.get(i).put(cur+"",count.get(i).get(cur+"")+1);
				}else{
					count.get(i).put(cur+"",1);
				}
			}
		}
		//print
		System.out.print("pos");
		for(String letter: alphabet){
			System.out.print("\t"+letter);
		}
		System.out.println();
		for(int i=0;i<maxLength;i++){
			System.out.print(i+1);
			final HashMap<String, Integer> tmp= count.get(i);
			for(String letter: alphabet){
				if(tmp.containsKey(letter)){
					System.out.print("\t"+tmp.get(letter));
				}else{
					System.out.print("\t0");
				}
			}
			System.out.println();
		}
	}
	
	public static void cleanXs(String faFile,String qualFile,int minLength,String outPrefix)throws Exception{
		fastaParser fp= new fastaParser(new BufferedReader(new FileReader(faFile)));
		RocheQualParser qp= new RocheQualParser(new BufferedReader(new FileReader(qualFile)));
		BufferedWriter fastaOut= new BufferedWriter(new FileWriter(outPrefix+".fas"));
		BufferedWriter qualOut= new BufferedWriter(new FileWriter(outPrefix+".fas.qual"));
		FastaSeq fs;
		RocheQualSeq qs;
		boolean X;
		int start=-1;
		ArrayList<int[]> toPrint;
		for(;fp.hasNext()&&qp.hasNext();){
			fs=fp.next();
			qs=qp.next();
			final String qname=fs.getQname();
			if(!qname.equals(qs.getQname())){
				System.err.println("Unmatched qnames: "+qname+"\t"+qs.getQname());
			}
			final String seq=fs.getSeq();
			toPrint= new ArrayList<int[]>();
			//parse for regions without X
			if(seq.charAt(0)=='X'){
				X=true;
			}else{
				X=false;
				start=0;
			}
			for(int i=1;i<seq.length();i++){
				if(seq.charAt(i)=='X'){
					if(!X){
						if(i-start>minLength){
							toPrint.add(new int[] {start,i});
						}
						X=true;
					}
				}else{
					if(X){
						start=i;
						X=false;
					}
				}
			}
			if(!X&&seq.length()-start>minLength){
				toPrint.add(new int[] {start,seq.length()});
			}
			int count=1;
			for(int[] pair : toPrint){
				fastaOut.write(">"+qname+"_"+count+"\n"+fs.getSeq(pair[0],pair[1])+"\n");
				qualOut.write(">"+qname+"_"+count+"\n"+qs.getSeq(pair[0], pair[1])+"\n");
				count++;
			}
		}
		fastaOut.close();
		qualOut.close();
	}
	
	public static void drop2andfilter(String faFile)throws Exception{
		fastaParser fp= new fastaParser(new BufferedReader(new FileReader(faFile)));
		FastaSeq fs;
		String nSeq;
		for(;fp.hasNext();){
			fs=fp.next();
			System.out.println(fs.getHeader());
			nSeq="";
			for(int i=2;i<fs.getSeq().length();i++){
				if((i-1)%3!=0){
					nSeq+=fs.getSeq().charAt(i);
				}
			}
			System.out.println(nSeq);
		}
	}
	
	public static void oneLetterPerLine(String faFile)throws Exception{
		fastaParser fp= new fastaParser(new BufferedReader(new FileReader(faFile)));
		FastaSeq fs;
		for(;fp.hasNext();){
			fs=fp.next();
			System.out.println(fs.getHeader());
			for(int i=0;i<fs.length();i++){
				System.out.println(fs.getSeq().charAt(i));
			}
		}
	}
	
	public static void cut(String faFile,String posFile)throws Exception{
		HashMap<String, ArrayList<String>> limits= new HashMap<String, ArrayList<String>>();
		BufferedReader in = new BufferedReader(new FileReader(posFile));
		int start,stop;
		String[] l;
		FastaSeq fs;
		for(String s=in.readLine();s!=null;s=in.readLine()){
			l=s.split("\t");
			if(l.length>4){
				if(!limits.containsKey(l[0])){
					limits.put(l[0], new ArrayList<String>());
				}
				limits.get(l[0]).add(s);
			}
		}
		fastaParser fp= new fastaParser(new BufferedReader(new FileReader(faFile)));
		for(;fp.hasNext();){
			fs=fp.next();
			if(limits.containsKey(fs.getQname())){
				final ArrayList<String> tmp= limits.get(fs.getQname());
				for (String s : tmp) {
					l=s.split("\t");
					
					System.out.println(">"+l[4]);
					start=Integer.parseInt(l[1]);
					stop=Integer.parseInt(l[2]);
					if(start<0){
						start=0;
					}
					if(stop>fs.length()){
						stop=fs.length();
					}
					if(l[3].equals("-")){
						//reverse complement
						System.out.println(sequenceUtils.reverseComplement(fs.getSeq(start,stop)));
					}else if(l[3].equals("+")){
						System.out.println(fs.getSeq(start,stop));
					}else{
						throw new Exception("unknown strand: "+l[3]+" (Accepts + or -)");
					}
				}
			}
		}
	}
	
	public static void genomicPrimerSite(String faFile)throws Exception{
		boolean found;
		fastaParser fp= new fastaParser(new BufferedReader(new FileReader(faFile)));
		FastaSeq fs;
		String seq;
		for(;fp.hasNext();){
			fs=fp.next();
			seq=fs.getSeq();
			//Check
			found=false;
			for(int i=0;i<seq.length()&&i<43;i++){
				if(checkPrimeSite(seq.substring(i, seq.length()<50?seq.length():50))){
					found=true;
					break;
				}
			}
			if(found){
				System.out.println(fs.getQname()+sep+"1"+sep+seq);
			}else{
				System.out.println(fs.getQname()+sep+"0"+sep+seq);
			}
		}
	}
	private static boolean checkPrimeSite(String s){
		int min=8,count=0;
		for(int i=0;i<s.length();i++){
			count+=s.charAt(i)=='A'?1:-1;
			if(count>min){
				break;
			}
		}
		if(count>min){
			return true;
		}else{
			return false;
		}
	}
	
	public static void orfStat(String faFile, boolean both)throws Exception{
		BufferedReader in = new BufferedReader(new FileReader(faFile));
		String s,qname="";
		for(s= in.readLine();s!=null;s=in.readLine())
			if(s.length()>0)
				if(s.charAt(0)=='>')
					break;
		
		ArrayList<ArrayList<Integer>> orfs=new ArrayList<ArrayList<Integer>>();
		int[] rStop= new int[3];
		String[] codons= new String[3];
		
		codons[2]=" ";
		
		int pos=0,frame=0;
		char c;
		for(;s!=null;s=in.readLine()){
			if(s.length()>0){
				if(s.charAt(0)=='>'){
					//print old 
					for (ArrayList<Integer> rf : orfs) {
						for (Integer start : rf) {
							System.out.println(qname+sep+((pos-start)/3)+sep+start+sep+pos);
						}
					}
					//initiate new
					qname=s.substring(1).split(" ")[0];
					orfs= new ArrayList<ArrayList<Integer>>();
					orfs.add(new ArrayList<Integer>());
					orfs.add(new ArrayList<Integer>());
					orfs.add(new ArrayList<Integer>());
					rStop= new int[] {0,0,0};
					codons= new String[3];
					codons[2]=" ";
					pos=0;
					frame=0;
				}else{
					int i=0;
					if(pos<2){
						for(c=s.charAt(i);pos<2&&i<s.length()-1;c=s.charAt(++i)){
							if(c!=' '){
								codons[2]=codons[2]+c;
								pos++;
							}
						}
					}
					for(;i<s.length();i++){
						c=s.charAt(i);
						if(c!=' '){
							if(frame==0){
								codons[0]=codons[2].substring(1, 3)+c;
							}else{
								codons[frame]=codons[frame-1].substring(1, 3)+c;
							}
							//forward
							if(codons[frame].equals("ATG")){
								//add new start
								orfs.get(frame).add(pos);
							}
							if(codons[frame].toUpperCase().equals("TAA")||codons[frame].toUpperCase().equals("TAG")||codons[frame].toUpperCase().equals("TGA")){
								//terminate all orfs in reading frame
								for (Integer start : orfs.get(frame)) {
									System.out.println(qname+sep+((pos-start)/3)+sep+start+sep+pos);
								}
								orfs.set(frame, new ArrayList<Integer>());
							}
							if(both){
								if(codons[frame].equals("CAT")){
									//print a new orf
									System.out.println(qname+sep+((pos-rStop[frame])/3)+sep+pos+sep+rStop[frame]);
								}
								if(codons[frame].toUpperCase().equals("TTA")||codons[frame].toUpperCase().equals("CTA")||codons[frame].toUpperCase().equals("TCA")){
									rStop[frame]=pos;
								}
							}
							if(frame==2){
								frame=0;
							}else{
								frame++;
							}
							pos++;
						}
					}
				}
			}
		}
		for (ArrayList<Integer> rf : orfs) {
			for (Integer start : rf) {
				System.out.println(qname+sep+((pos-start)/3)+sep+start+sep+pos);
			}
		}
	}
	
	public static void orfStat2(String faFile,boolean both)throws Exception{
		fastaParser fp= new fastaParser(new BufferedReader(new FileReader(faFile)));
		FastaSeq fs;
		for(;fp.hasNext();){
			fs=fp.next();
			String[] codons= new String[3];
			codons[2]=" "+fs.getSeq().charAt(0)+fs.getSeq().charAt(1);
			ArrayList<ArrayList<Integer>> orfs= new ArrayList<ArrayList<Integer>>();
			orfs.add(new ArrayList<Integer>());
			orfs.add(new ArrayList<Integer>());
			orfs.add(new ArrayList<Integer>());
			
			for(int i=3, j=0;i<fs.getSeq().length();i++,j=i%3){
				//forward
				if(j==0){
					codons[j]=codons[2].substring(1, 3)+fs.getSeq().charAt(i);
				}else{
					codons[j]=codons[j-1].substring(1, 3)+fs.getSeq().charAt(i);
				}				
				if(codons[j].toUpperCase().equals("ATG")){
					//start new orf
					orfs.get(j).add(i);
				}
				if(codons[j].toUpperCase().equals("TAA")||codons[j].toUpperCase().equals("TAG")||codons[j].toUpperCase().equals("TGA")){
					//end all orfs in reading frame
					for (Integer start : orfs.get(j)) {
						System.out.println(fs.getQname()+sep+((i-start)/3)+sep+start+sep+i);
					}
					orfs.set(j, new ArrayList<Integer>());
				}
			}
			//empty unfinished orfs
			for (ArrayList<Integer> rf : orfs) {
				int i=fs.getSeq().length();
				for (Integer start : rf) {
					System.out.println(fs.getQname()+sep+((i-start)/3)+sep+start+sep+i);
				}
			}
			if(both){
				codons= new String[3];
				codons[2]=" "+sequenceUtils.complement(fs.getSeq().charAt(fs.getSeq().length()-1))+sequenceUtils.complement(fs.getSeq().charAt(fs.getSeq().length()-2));
				orfs= new ArrayList<ArrayList<Integer>>();
				orfs.add(new ArrayList<Integer>());
				orfs.add(new ArrayList<Integer>());
				orfs.add(new ArrayList<Integer>());
				for(int i=fs.getSeq().length()-3,j=0;i>=0;i--){
					if(j==0){
						codons[j]=codons[2].substring(1, 3)+sequenceUtils.complement(fs.getSeq().charAt(i));
					}else{
						codons[j]=codons[j-1].substring(1, 3)+sequenceUtils.complement(fs.getSeq().charAt(i));
					}
					if(codons[j].toUpperCase().equals("ATG")){
						//start new orf
						orfs.get(j).add(i);
					}
					if(codons[j].toUpperCase().equals("TAA")||codons[j].toUpperCase().equals("TAG")||codons[j].toUpperCase().equals("TGA")){
						//terminate all orfs in reading frame
						for (Integer start : orfs.get(j)) {
							System.out.println(fs.getQname()+sep+((start-i)/3)+sep+start+sep+i);
						}
						orfs.set(j, new ArrayList<Integer>());
					}
					if(j==2){
						j=0;
					}else{
						j++;
					}
				}
				//empty unfinished orfs
				for (ArrayList<Integer> rf : orfs) {
					int i=0;
					for (Integer start : rf) {
						System.out.println(fs.getQname()+sep+((i-start)/3)+sep+start+sep+i);
					}
				}
			}
		}
	}
	
//	private static char complement(char c)throws Exception{
//		switch (c) {
//		case 't':
//			return 'a';
//		case 'T':
//			return 'A';
//		case 'a':
//			return 't';
//		case 'A':
//			return 'T';
//		case 'g':
//			return 'c';
//		case 'G':
//			return 'C';
//		case 'c':
//			return 'g';
//		case 'C':
//			return 'G';
//		case 'N':
//			return 'N';
//		case 'n':
//			return 'n';
//			default:
//				return c;
////			throw new Exception("Unknown character: "+c);
//		}
//	}
	
	public static void rename(String faFile,String giTableFile)throws Exception{
		BufferedReader in= new BufferedReader(new FileReader(giTableFile));
		HashMap<String, String> code= new HashMap<String, String>();
		String[] l;
		for(String s=in.readLine();s!=null;s=in.readLine()){
			l=s.split("\t");
			if(l.length>1){
				code.put(l[0], l[1]);
			}
		}
		fastaParser fp= new fastaParser(new BufferedReader(new FileReader(faFile)));
		FastaSeq fs;
		for(;fp.hasNext();){
			fs=fp.next();
			if(code.containsKey(fs.getQname())){
				System.out.println(">"+code.get(fs.getQname())+" "+fs.getQname());
				System.out.println(fs.getSeq());
			}else{
				System.out.println(fs.toString());
			}
		}
	}
	
	public static void extractSubSeq(String faFile, BufferedReader in,int flank)throws Exception{
		RandomAccessFastaParser rafp= new RandomAccessFastaParser(faFile);
		FastaSeq fs= new FastaSeq(">notAsequence", "ATGC");
		String[] l;
		int start,stop;
		for(String s=in.readLine();s!=null;s=in.readLine()){
			l=s.split("\t");
			if(l.length>=3){
				if(fs.getQname()!=l[0]){
					fs=rafp.get(l[0]);
				}
				System.out.println(">"+fs.getQname()+"_"+l[1]+":"+l[2]);
				start=Integer.parseInt(l[1]);
				stop=Integer.parseInt(l[2]);
				System.out.println(fs.getSeq(start-flank>0?start-flank:0, stop+flank<fs.length()?stop+flank:fs.length()));
			}
		}
		in.close();
	}
	
	public static void alignedaa2rna_shiftFrame(String fafile,String giList)throws Exception{
		BufferedReader in =new BufferedReader(new FileReader(giList));
		ArrayList<String> gis= new ArrayList<String>();
		fastaParser fp = new fastaParser(new BufferedReader(new FileReader(fafile)));
		FastaSeq fs;
		for(String s=in.readLine();s!=null&&s.length()>0;s=in.readLine()){
			gis.add(s);
		}
		for(;fp.hasNext();){
			fs=fp.next();
			if(gis.contains(fs.getQname())){
				System.out.println(fs.getHeader());
				System.out.println("N"+fs.getSeq());
			}else{
				System.out.println(fs.toString());
			}
		}
	}
	
	public static void giList(String faFile)throws Exception{
		fastaParser fp= new fastaParser(new BufferedReader(new FileReader(faFile)));
		for(;fp.hasNext();){
			System.out.println(fp.next().getQname());
		}
	}
	
//	public static void renameFromM8blast(String faFile, String blastFile)throws Exception{
//		HashMap<String, ArrayList<String>> map= new HashMap<String, ArrayList<String>>();
//		blastM8Parser bp= new blastM8Parser(new BufferedReader(new FileReader(blastFile)));
//		blastM8Alignment ba;
//		for(;bp.hasMore();){
//			ba=bp.nextHit();
//			if(!map.containsKey(ba.tname)){
//				map.put(ba.tname, new ArrayList<String>());
//			}
//			map.get(ba.tname).add(ba.qname);
//		}
//		fastaParser fp= new fastaParser(new BufferedReader(new FileReader(faFile)));
//		FastaSeq fs;
//		for(;fp.hasNext();){
//			fs=fp.next();
//			if(map.containsKey(fs.getQname())){
//				for (String qname : map.get(fs.getQname())) {
//					System.out.println(">"+qname+" "+fs.getQname());
//					System.out.println(fs.getSeq());
//				}
//			}
//		}
//	}
	
	public static void toStockholm(String faFile)throws Exception{
		//int namelength=19;
		System.out.println("# STOCKHOLM 1.0");
		fastaParser fp=new fastaParser(new BufferedReader(new FileReader(faFile)));
		FastaSeq fs;
		for(;fp.hasNext();){
			fs=fp.next();
//			qname=fs.getHeader().split(" ")[0].substring(1);
//			for(;qname.length()<namelength;){
//				qname+=" ";
//			}
//			if(qname.length()>namelength){
//				qname=qname.substring(0,namelength);
//			}
			System.out.println(fs.getQname()+"\t"+fs.getSeq());
		}
		System.out.println("//");
	}
	
//	public static void pfamDomains(String faFile, String resultFile, String domainSet, boolean inSet,int top)throws Exception{
//		//Get domains to get
//		ArrayList<String> domainsToUse= new ArrayList<String>();
//		BufferedReader in= new BufferedReader(new FileReader(domainSet));
//		for(String s=in.readLine();s!=null;s=in.readLine()){
//			if(s.length()>0){
//				domainsToUse.add(s);
//			}
//		}
//		//Get domain hit results
//		HashMap<String, ArrayList<pfamAlignment>> domains= new HashMap<String, ArrayList<pfamAlignment>>();
//		pfamParser pp= new pfamParser(resultFile);
//		pfamAlignment tmp;
//		for(;pp.hasMore();){
//			tmp=pp.next();
//			if(domainsToUse.contains(tmp.getHmmname())==inSet){
//				if(!domains.containsKey(tmp.getQname())){
//					domains.put(tmp.getQname(), new ArrayList<pfamAlignment>());
//				}
//				//order domains
//				int p=domains.get(tmp.getQname()).size();
//				if(p>0){
//					if(domains.get(tmp.getQname()).get(0).getEvalue()>tmp.getEvalue()){
//						p=0;
//					}else{
//						for(;p>1&&domains.get(tmp.getQname()).get(p-1).getEvalue()>tmp.getEvalue();p--);
//						
//					}
//				}
//				domains.get(tmp.getQname()).add(p, tmp);
//			}
//		}
//		//parse fasta file and extract domains
//		fastaParser fp= new fastaParser(new BufferedReader(new FileReader(faFile)));
//		FastaSeq fs;
//		String qname;
//		int i;
//		for(;fp.hasNext();){
//			fs=fp.next();
//			qname=fs.getHeader().split(" ")[0].substring(1);
//			if(domains.containsKey(qname)){
//				i=0;
//				for (pfamAlignment domain : domains.get(qname)) {
//					i++;
//					if(i>top){
//						System.err.println(qname+sep+domains.get(qname).size()+" extra domains");
//						break;
//					}
//					System.out.println(fs.getHeader()+" "+domain.getHmmname()+" "+i);
//					System.out.println(fs.getSeq().substring(domain.getQstart()-1, domain.getQend()));
//				}
//			}else{
////				System.err.println("no domain for: "+qname );
//			}
//		}
//	}

	public static void sort(String faFile, String listFile)throws Exception{
		ArrayList<String> list= new ArrayList<String>();
		BufferedReader in= new BufferedReader(new FileReader(listFile));
		for(String s=in.readLine();s!=null;s=in.readLine()){
			if(s.length()>0){
				list.add(s.split("\t")[0]);
			}
		}

		//extract fasta-sequences
		HashMap<String, FastaSeq> faSet=new HashMap<String, FastaSeq>();
		fastaParser fp=new fastaParser(new BufferedReader(new FileReader(faFile)));
		FastaSeq tmp;
		for(;fp.hasNext();){
			tmp=fp.next();
			if(list.contains(tmp.getHeader().substring(1).split(" ")[0])){
				faSet.put(tmp.getHeader().substring(1).split(" ")[0], tmp);
			}
		}

		//print sequences
		for (String s : list) {
			System.out.println(faSet.get(s).toString());
		}
	}

	public static void getEditedRNA(String filelist)throws Exception{
		BufferedReader in= new BufferedReader(new FileReader(filelist));
		fastaParser fp;
		FastaSeq fs;
		String seq;
		for(String s= in.readLine();s!=null;s= in.readLine()){
			try{
				fp= new fastaParser(new BufferedReader(new FileReader(s)));
				System.out.println(">"+s.substring(0, s.lastIndexOf(".")));
				seq="";
				for(;fp.hasNext();){
					fs=fp.next();
					if(fs.getHeader().startsWith(">exon")&& !fs.getHeader().endsWith("trans")){
						seq+=fs.getSeq();
					}
				}
				System.out.println(seq);
			}catch (Exception e) {
				System.err.println(s);
				e.printStackTrace(System.err);
			}
		}
	}

	public static void alignedaa2rna(String aaFile, String rnaFile, boolean checkAll)throws Exception{
		char gap='-',c;
		BufferedWriter out=new BufferedWriter(new OutputStreamWriter(System.err));
		boolean error=false;
		HashMap<String, String[]> rna= new HashMap<String, String[]>();
		String[] l,m;
		String seq,tmp;
//		int frame;
		int start,nr;
		FastaSeq fs;
		fastaParser fp= new fastaParser(new BufferedReader(new FileReader(rnaFile)));
		for(;fp.hasNext();){
			fs= fp.next();
			l=new String[] {fs.getSeq(),translateRNA(fs.getSeq(),1),translateRNA(fs.getSeq(),2),translateRNA(fs.getSeq(),3)};
//			System.err.println(test=fs.getHeader().substring(1).split(" ")[0]);
			rna.put(fs.getHeader().substring(1).split(" ")[0], l);
		}
		fp= new fastaParser(new BufferedReader(new FileReader(aaFile)));
		for(;fp.hasNext();){
			fs=fp.next();
			tmp=fs.getHeader();
			if(rna.containsKey(tmp.substring(1).split(" ")[0])){
				l=rna.get(tmp.substring(1).split(" ")[0]);
				System.out.println(tmp);
				seq=fs.getSeq();
				//find start
//				test=new String[] {"","",""};
				tmp="";
				start=-1;
//				frame=1;
				for(int i=0;i<seq.length();i++){
					nr=0;
//					for(int j=1;j<4;j++){
					c=seq.charAt(i);
					if(c!=gap){
						tmp+=c;
//						test[j-1]+=c;
						if(tmp.length()>2){
							for(int j=1;j<4;j++){
								m=l[j].split(tmp);
								if(m.length==2){
									start=m[0].length()*3+j-1;
//									frame=j;
									nr++;
								}
								if(m.length>2){
//									System.err.println("tjo"+j+"z"+i+"nr"+nr);
									nr+=2;
								}
							}


							if(nr==1){
//								System.out.println(frame+"\t"+tmp+"\t"+start);
								break;
							}else{
								start=-1;
							}
						}
					}
				}
				
				
//				for(;frame<4;frame++){
//					tmp="";
//					for(int i=0;i<seq.length();i++){
//						c=seq.charAt(i);
//						if(c!=gap){
//							tmp+=c;
//							if(tmp.length()>2){
//								m=l[frame].split(tmp);
//								if(m.length==2){
//									start=m[0].length()*3;
//									break;
//								}
//							}
//						}
//					}
//					if(start!=-1){
//						break;
//					}
//				}
				if(start==-1){
					//could not find the start of the given subsequence
					System.err.println("could not find the start in the rna of the prot sequence:"+fs.getHeader());
					if(!error){
						out=new BufferedWriter(new FileWriter("err.out"));
						error=true;
					}
					out.write(fs.getQname()+"\n");
				}else{
					//print rna sequence
					for(int i=0,j=0;i<seq.length();i++){
						c=seq.charAt(i);
						if(c==gap){
							System.out.print(""+gap+gap+gap);
						}else{
							tmp=l[0].substring(start+j, start+(j+=3));
							if(checkAll){
								if(!translateRNA(tmp).equals(c+"")){
									//error
									System.err.println("error in the translation of "+fs.getHeader()+": "+translateRNA(tmp)+"!="+c);
									if(!error){
										out=new BufferedWriter(new FileWriter("err.out"));
										error=true;
									}
									out.write(fs.getQname()+"\n");
								}
							}
							System.out.print(tmp);
						}
					}
				}
				System.out.println();
			}else{
				System.err.println("Missing among the rna-sequences: \n"+tmp);
			}
		}
		out.close();
	}

	public static String translateRNA(String rna){
		return translateRNA(rna,1);
	}
	public static String translateRNA(String rna,int frame){
		//create dictionary
		HashMap<String, String> codon=standardCodon();
		String prot="";
		String c;
		if(frame>0){
			for(int i=frame-1;i<rna.length();i+=3){
				if(i+3==rna.length()){
					c=rna.substring(i);
				}else if(i+3>rna.length()){
					c="END";
				}else {
					c=rna.substring(i,i+3);
				}
				if(codon.containsKey(c.toUpperCase())){
					prot+=codon.get(c.toUpperCase());
				}else{
					prot+="X";
				}
			}
		}else{
			for(int i=rna.length()+frame;i>=0;i-=3){
				c=""+rna.charAt(i+2)+rna.charAt(i+1)+rna.charAt(i);
				if(codon.containsKey(c.toUpperCase())){
					prot+=codon.get(c.toUpperCase());
				}else{
					prot+="X";
				}
			}
		}
		return prot;
	}

	private static HashMap<String, String> standardCodon(){
		HashMap<String, String> codon=new HashMap<String, String>();
		codon.put("TTT", "F");
		codon.put("TTC", "F");
		codon.put("TTA", "L");
		codon.put("TTG", "L");
		codon.put("TCT", "S");
		codon.put("TCC", "S");
		codon.put("TCA", "S");
		codon.put("TCG", "S");
		codon.put("TAT", "Y");
		codon.put("TAC", "Y");
		codon.put("TAA", "*");
		codon.put("TAG", "*");
		codon.put("TGT", "C");
		codon.put("TGC", "C");
		codon.put("TGA", "*");
		codon.put("TGG", "W");
		codon.put("CTT", "L");
		codon.put("CTC", "L");
		codon.put("CTA", "L");
		codon.put("CTG", "L");
		codon.put("CCT", "P");
		codon.put("CCC", "P");
		codon.put("CCA", "P");
		codon.put("CCG", "P");
		codon.put("CAT", "H");
		codon.put("CAC", "H");
		codon.put("CAA", "Q");
		codon.put("CAG", "Q");
		codon.put("CGT", "R");
		codon.put("CGC", "R");
		codon.put("CGA", "R");
		codon.put("CGG", "R");
		codon.put("ATT", "I");
		codon.put("ATC", "I");
		codon.put("ATA", "I");
		codon.put("ATG", "M");
		codon.put("ACT", "T");
		codon.put("ACC", "T");
		codon.put("ACA", "T");
		codon.put("ACG", "T");
		codon.put("AAT", "N");
		codon.put("AAC", "N");
		codon.put("AAA", "K");
		codon.put("AAG", "K");
		codon.put("AGT", "S");
		codon.put("AGC", "S");
		codon.put("AGA", "R");
		codon.put("AGG", "R");
		codon.put("GTT", "V");
		codon.put("GTC", "V");
		codon.put("GTA", "V");
		codon.put("GTG", "V");
		codon.put("GCT", "A");
		codon.put("GCC", "A");
		codon.put("GCA", "A");
		codon.put("GCG", "A");
		codon.put("GAT", "D");
		codon.put("GAC", "D");
		codon.put("GAA", "E");
		codon.put("GAG", "E");
		codon.put("GGT", "G");
		codon.put("GGC", "G");
		codon.put("GGA", "G");
		codon.put("GGG", "G");
		codon.put("END", "");
		return codon;
	}

	private static void removeDesc(String fafile)throws Exception{
		fastaParser fa=new fastaParser(new BufferedReader(new FileReader(fafile)));
		FastaSeq tmp;
		for(;fa.hasNext();){
			tmp=fa.next();
			System.out.println(tmp.getHeader().split(" ")[0]+"\n"+tmp.getSeq());
		}
	}

	private static void cutColumns(String fafile,double covCut)throws Exception{
		ArrayList<FastaSeq> seqs =new ArrayList<FastaSeq>();
		String[] newSeq;
		FastaSeq tmp;
		int maxLength=0,count;
		fastaParser fa=new fastaParser(new BufferedReader(new FileReader(fafile)));
		for(;fa.hasNext();){
			tmp=fa.next();
			maxLength=tmp.getSeq().length()>maxLength?tmp.getSeq().length():maxLength;
			seqs.add(tmp);
		}
		newSeq=new String[seqs.size()];
		for(int i=0;i<newSeq.length;i++){
			newSeq[i]="";
		}
		if(covCut<=1){
			covCut*=seqs.size();
		}
		covCut=seqs.size()-covCut; //the maximum number of dashes allowed
		for(int i=0;i<maxLength;i++){
			count=0;
			for (FastaSeq seq : seqs) {
				if(seq.getSeq().length()<i){
					count++;
				}else if(seq.getSeq().charAt(i)=='-'||seq.getSeq().charAt(i)==' '){
					count++;
				}
			}
			if(count<covCut){
				for(int j=0;j<seqs.size();j++){
					tmp=seqs.get(j);
					if(tmp.getSeq().length()>i){
						newSeq[j]+=""+tmp.getSeq().charAt(i);
					}else{
						newSeq[j]+="-";
					}
				}
			}
		}
		for(int i=0;i<seqs.size();i++){
			System.out.println(seqs.get(i).getHeader()+"\n"+newSeq[i]);
		}
	}

//	private static void toPhy(String fafile)throws Exception{
//		fastaParser fa=new fastaParser(new BufferedReader(new FileReader(fafile)));
//		ArrayList<FastaSeq> seqs=new ArrayList<FastaSeq>();
//		for(;fa.hasNext();){
//			seqs.add(fa.next());
//		}
//		PhylipSet tmp=new PhylipSet(seqs);
//		tmp.print2Screen();
//	}
//
//	private static void toNexus(String fafile)throws Exception{
//		fastaParser fa=new fastaParser(new BufferedReader(new FileReader(fafile)));
//		ArrayList<FastaSeq> seqs=new ArrayList<FastaSeq>();
//		for(;fa.hasNext();){
//			seqs.add(fa.next());
//		}
//		NexusSet tmp=new NexusSet(seqs,true);
//		tmp.print2Screen();
//	}
	

	private static void subFiles(String fafile,String subSetFile,String outPrefix)throws Exception{
		HashMap<String, ArrayList<String>> groups= new HashMap<String, ArrayList<String>>();
		HashMap<String, BufferedWriter> outs=new HashMap<String, BufferedWriter>();
		BufferedReader in=new BufferedReader(new FileReader(subSetFile));
		String[] l;
		String header,id;
		for(String s=in.readLine();s!=null;s=in.readLine()){
			l=s.split("\t");
			if(l.length>=2){
				if(!groups.containsKey(l[0])){
					groups.put(l[0], new ArrayList<String>());
				}
				if(!groups.get(l[0]).contains(l[1])){
					if(!outs.containsKey(l[1])){
						outs.put(l[1], new BufferedWriter(new FileWriter(outPrefix+"_"+l[1]+".fa")));
					}
					groups.get(l[0]).add(l[1]);
				}
			}
		}
		fastaParser fa=new fastaParser(new BufferedReader(new FileReader(fafile)));
		for(;fa.hasNext();){
			header=fa.nextHit();
			id=header.substring(1).split(" ")[0];
			if(groups.containsKey(id)){
				for ( String group : groups.get(id)) {
					outs.get(group).write(header+"\n"+fa.getSeq()+"\n");
				}
			}
		}
		for (BufferedWriter out: outs.values()) {
			out.close();
		}
	}

	private static void countTagged(String fafile)throws Exception{
		fastaParser fa=new fastaParser(new BufferedReader(new FileReader(fafile)));
		HashMap<String, Integer> nr=new HashMap<String, Integer>();
		String tag;
		for(;fa.hasNext();){
			tag=fa.nextHit().split("_")[0].substring(1);
			if(nr.containsKey(tag)){
				nr.put(tag, nr.get(tag)+1);
			}else{
				nr.put(tag, 1);
			}
		}
		//print result
		for(String s:nr.keySet()){
			System.out.println(s+"\t"+nr.get(s));
		}
	}

	public static boolean phylipWarning(ArrayList<String> list){
		HashMap<String, ArrayList<String>> groups =new HashMap<String, ArrayList<String>>();
		boolean duplicates=false;
		
		String phy;
		for(String header: list){
			phy=header.length()>10?header.substring(0, 11):header;
			if(!groups.containsKey(phy)){
				groups.put(phy, new ArrayList<String>());
			}
			groups.get(phy).add(header);
		}
		//print result
		for (ArrayList<String> group : groups.values()) {
			if(group.size()>1){
				duplicates=true;
				for (String s : group) {
					System.err.println(s);
				}
				System.err.println();
			}
		}
		return duplicates;
	}
	
	
	
	private static void phylipWarning(String fafile)throws Exception{
		fastaParser fa=new fastaParser(new BufferedReader(new FileReader(fafile)));
		HashMap<String, ArrayList<String>> groups =new HashMap<String, ArrayList<String>>();
		String header,phy;
		for(;fa.hasNext();){
			header=fa.nextHit();
			phy=header.length()>10?header.substring(0, 11):header;
			if(!groups.containsKey(phy)){
				groups.put(phy, new ArrayList<String>());
			}
			groups.get(phy).add(header);
		}
		//print result
		for (ArrayList<String> group : groups.values()) {
			if(group.size()>1){
				for (String s : group) {
					System.out.println(s);
				}
				System.out.println();
			}
		}
	}

	private static void shortenNames(String fafile, String regExp) throws Exception{
		fastaParser fa=new fastaParser(new BufferedReader(new FileReader(fafile)));
		String header;
		for(;fa.hasNext();){
			header=fa.nextHit();
			System.out.println(header.replaceFirst(regExp, "")+"\n"+fa.getSeq());
		}
	}

	private static void reduceTaggedByChance(String fafile, int N) throws Exception{
		fastaParser fa=new fastaParser(new BufferedReader(new FileReader(fafile)));
		Hashtable<String, ArrayList<String[]>> seqs=new Hashtable<String, ArrayList<String[]>>();
		String header,tag;
		String[] seq;
		for(;fa.hasNext();){
			header=fa.nextHit();
			tag=header.substring(1).split("_")[0];
			seq=new String[] {header,fa.getSeq()};
			if(!seqs.containsKey(tag)){
				seqs.put(tag, new ArrayList<String[]>());
			}
			seqs.get(tag).add(seq);
		}
		Random slump=new Random();
		ArrayList<String> used;
		for (ArrayList<String[]> group:seqs.values()) {
			double chance=0.1/group.size();
			used=new ArrayList<String>();
			int n=group.size()<N?group.size():N;
			while(used.size()<n){
				for (String[] s : group) {
					if(!used.contains(s[0])){
						if(slump.nextDouble()<=chance){
							System.out.println(s[0]);
							System.out.println(s[1]);
							used.add(s[0]);
							if(used.size()>=n){
								break;
							}
						}
					}
				}
			}
		}
	}

	public static void lengths(String fafile)throws Exception{
		fastaParser fa=new fastaParser(new BufferedReader(new FileReader(fafile)));
		for(;fa.hasNext();){
			System.out.println(fa.nextHit().substring(1).split(" ")[0]+"\t"+fa.getSeq().length());
		}
	}

	public static void removeDash(String fafile)throws Exception{
		fastaParser fa=new fastaParser(new BufferedReader(new FileReader(fafile)));
		String header, seq, seqclean;
		for(;fa.hasNext();){
			header=fa.nextHit();
			seq=fa.getSeq();
			seqclean="";
			for(int i=0;i<seq.length();i++){
				if(!(seq.charAt(i)=='-')){
					seqclean+=seq.charAt(i);
				}
			}
			System.out.println(header+"\n"+seqclean);
		}
	}
	public static void chance(String fafile,int N)throws Exception{
		double chance=1.0/(N*10.0);
		Random slump=new Random();
		fastaParser fa; 
		ArrayList<String> used=new ArrayList<String>();
		for(int n=0,i=0;n<N && i<100;i++){
			fa = new fastaParser(new BufferedReader(new FileReader(fafile)));
			for(;fa.hasNext();){
				String s=fa.nextHit();
				if(!used.contains(s)){
					if(slump.nextDouble()<=chance){
						System.out.println(s+"\n"+fa.getSeq());
						used.add(s);
						n++;
					}
				}
				if(n>=N){
					break;
				}
			}
		}
	}

	public static void diffTag(String fafile,String tagfile)throws Exception{
		String leftoverTag="919919";
		BufferedReader in=new BufferedReader(new FileReader(tagfile));
		Hashtable<String, String> tags=new Hashtable<String, String>();
		for(String s=in.readLine();s!=null;s=in.readLine()){
			String[] l=s.split("\t");
			if(l.length>=2){
				tags.put(l[0], l[1]);
			}
		}
		fastaParser fa=new fastaParser(new BufferedReader(new FileReader(fafile)));
		for(;fa.hasNext();){
			String s=fa.nextHit();
			String id=s.substring(1).split(" ")[0];//which part of the header is the id? change here
			if(tags.containsKey(id)){
				System.out.println(">"+tags.get(id)+"_"+s.substring(1)+"\n"+fa.getSeq());
			}else{
				System.out.println(">"+leftoverTag+"_"+s.substring(1)+"\n"+fa.getSeq());
			}
		}
	}

	public static void dropShort(String fafile, int n)throws Exception{
		fastaParser fa=new fastaParser(new BufferedReader(new FileReader(fafile)));
		for(;fa.hasNext();){
			String s=fa.nextHit();
			String seq=fa.getSeq();
			if(seq.length()>=n){
				System.out.println(s+"\n"+seq);
			}
		}
	}
	public static void tag(BufferedReader fafile, String tag)throws Exception{
		fastaParser fa =new fastaParser(fafile);
		for(;fa.hasNext();){
			String s=fa.nextHit();
			System.out.println(">"+tag+"_"+s.substring(1)+"\n"+fa.getSeq());
		}
	}

	public static void cleanse(String fafile)throws Exception{
		HashSet<String> v=new HashSet<String>();
		fastaParser fa =new fastaParser(new BufferedReader(new FileReader(fafile)));
		for(;fa.hasNext();){
			String s=fa.nextHit();
			if(!v.contains(s)){
				System.out.println(s+"\n"+fa.getSeq());
				v.add(s);
			}
		}
	}
	
	public static void cleanseShort(String fafile)throws Exception{
		HashMap<String, FastaSeq> seqs= new HashMap<String, FastaSeq>();
		fastaParser fa =new fastaParser(new BufferedReader(new FileReader(fafile)));
		FastaSeq fs;
		for(;fa.hasNext();){
			fs=fa.next();
			if(seqs.containsKey(fs.getQname())){
				if(fs.length()>seqs.get(fs.getQname()).length()){
					seqs.put(fs.getQname(), fs);
				}
			}else{
				seqs.put(fs.getQname(), fs);
			}
		}
		for (FastaSeq seq : seqs.values()) {
			System.out.println(seq.toString());
		}
	}

	public static void description(String fafile, String gi_list)throws Exception{
		Hashtable<String, Integer> ht=new Hashtable<String, Integer>();
		BufferedReader in =new BufferedReader(new FileReader(gi_list));
		for(String s=in.readLine();s!=null;s=in.readLine()){
			if(ht.containsKey(s))
				ht.put(s,new Integer(((Integer)ht.get(s)).intValue()+1));
			else
				ht.put(s,new Integer(1));
		}
		fastaParser fa=new fastaParser(new BufferedReader(new FileReader(fafile)));
		for(;fa.hasNext();){
			String s=fa.nextHit();
			String gi=s.substring(1).split(" ")[0];
			if(ht.containsKey(gi)){
				int nr=((Integer)ht.get(gi)).intValue();
				for(int i=0;i<nr;i++)
					System.out.println(gi+"\t"+s.substring(s.indexOf(" ")+1));
			}
		}
	}

	public static void description(String fafile)throws Exception{
		fastaParser fa=new fastaParser(new BufferedReader(new FileReader(fafile)));
		for(;fa.hasNext();){
			String s=fa.nextHit();
			String gi=s.substring(1).split(" ")[0];
			System.out.println(gi+"\t"+s.substring(s.indexOf(" ")+1));
		}
	}

	public static void split(String fafile,String outprefix, int n)throws Exception{
		fastaParser fa=new fastaParser(new BufferedReader(new FileReader(fafile)));
		String postfix="fa";//fafile.split(".")[fafile.split(".").length-1];
		BufferedWriter out=new BufferedWriter(new FileWriter(outprefix+"_00001."+postfix));
		String nr;
		for(int printed=0,filenr=1;fa.hasNext();++printed){
			if(printed==n){
				out.close();
				++filenr;
				nr=""+filenr;
				while(nr.length()<5){
					nr="0"+nr;
				}
				out=new BufferedWriter(new FileWriter(outprefix+"_"+nr+"."+postfix));
				printed=0;
			}
			out.write(fa.nextHit()+"\n");
			out.write(fa.getSeq()+"\n");
		}
		out.close();
	}
	public static void subSet(String fafile,BufferedReader subSet,boolean inSubSet)throws Exception{
		HashSet<String> set=new HashSet<String>();
//		BufferedReader in=new BufferedReader(new FileReader(subFile));
		for(String s=subSet.readLine();s!=null;s=subSet.readLine())
			set.add(s.split("\t")[0]);
		fastaParser fa=new fastaParser(new BufferedReader(new FileReader(fafile)));
		for(;fa.hasNext();){
			String header=fa.nextHit();
//			System.out.println(header);
//			System.out.println(header.split(" ")[0].substring(1).replace('|', '#').split("#")[2]);
//			if(set.contains(header.split(" ")[0].substring(1).replace('|', '#').split("#")[2])==inSubSet){ //change here to decide which part of the header that is in the subList
			if(set.contains(header.split(" ")[0].substring(1))==inSubSet){ //change here to decide which part of the header that is in the subList
				System.out.println(header+"\n"+fa.getSeq()/*+"\n"*/);
			}
		}
	}
	public static void subSetReg(String fafile,String regExp,boolean inSubset)throws Exception{
		fastaParser fa=new fastaParser(new BufferedReader(new FileReader(fafile)));
		for(;fa.hasNext();){
			String header=fa.nextHit();
			if(header.matches(regExp)==inSubset){
				System.out.println(header+"\n"+fa.getSeq());
			}
		}
	}
	public static void subSet(String fafile, ArrayList<String> subList, boolean inSubSet,String outFile)throws Exception{
		fastaParser fa=new fastaParser(new BufferedReader(new FileReader(fafile)));
		BufferedWriter out=new BufferedWriter(new FileWriter(outFile));
		String header;
		for(;fa.hasNext();){
			header=fa.nextHit();
			if (subList.contains(header.split(" ")[0].substring(1))==inSubSet){ //change here to decide which part of the header that is in the subList
				out.write(header+"\n"+fa.getSeq()+"\n");
			}
		}
		out.close();
	}
//	public static void getDomains(String fafile,String hmmerFile,String domainSet,boolean inSubSet)throws Exception{
//		fastaParser fa=new fastaParser(new BufferedReader(new FileReader(fafile)));
//		Hashtable<String, String[]> seq=new Hashtable<String, String[]>();
//		for(;fa.hasNext();){
//			String header=fa.nextHit();
//			String key=header.split(" ")[0].substring(1); //change here to decide which part of the header that should represent the sequence
//			seq.put(key,new String[] {header,fa.getSeq()});
//		}
//		BufferedReader in=new BufferedReader(new FileReader(domainSet));
//		Vector<String> set=new Vector<String>();
//		for(String s=in.readLine();s!=null;s=in.readLine())
//			set.add(s);
//		hmmerParser hp=new hmmerParser(new BufferedReader(new FileReader(hmmerFile)));
//		for(;hp.hasMore();){
//			hmmerAlignment ha=hp.nextHit();
//			if(seq.containsKey(ha.qname)){
//				String[] currSeq=(String[])seq.get(ha.qname);
//				hmmerModelAlignment[] hma=ha.getTopDomainsByScore();
//				Hashtable<String, Integer> nr=new Hashtable<String, Integer>();
//				for(int i=0;i<hma.length;i++){
//					if(set.contains(hma[i].hmmname)==inSubSet){
//						if(nr.containsKey(hma[i].hmmname))
//							nr.put(hma[i].hmmname,new Integer(((Integer)nr.get(hma[i].hmmname)).intValue()+1));
//						else
//							nr.put(hma[i].hmmname,new Integer(1));
//						System.out.println(">"+hma[i].hmmname+"|"+ha.qname+"|"+nr.get(hma[i].hmmname));
//						System.out.println(currSeq[1].substring(hma[i].tstart-1,hma[i].tend-1));
//					}
//				}
//			}
//		}
//	}
//	public static void domain_clean(String fafile,String hmmerFile,String domainSet,boolean inSubSet)throws Exception{
//		fastaParser fa=new fastaParser(new BufferedReader(new FileReader(fafile)));
//		Hashtable<String, String[]> seq=new Hashtable<String, String[]>();
//		for(;fa.hasNext();){
//			String key=fa.nextHit(); 
//			seq.put(key.split(" ")[0].substring(1),new String[] {key,fa.getSeq()});//change here to decide which part of the header that should represent the sequence
//		}
//		BufferedReader in=new BufferedReader(new FileReader(domainSet));
//		Vector<String> set=new Vector<String>();
//		for(String s=in.readLine();s!=null;s=in.readLine())
//			set.add(s);
//		hmmerParser hp=new hmmerParser(new BufferedReader(new FileReader(hmmerFile)));
//		for(;hp.hasMore();){
//			hmmerAlignment ha=hp.nextHit();
//			if(seq.containsKey(ha.qname)){
//				String[] currSeq=(String[])seq.get(ha.qname);
//				hmmerModelAlignment[] hma=ha.getTopDomainsByScore();
//				int pos=0;
//				String newSeq="";
//				for(int i=0;i<hma.length;i++){
//					newSeq+=currSeq[1].substring(pos,hma[i].tstart-1);
//					if(set.contains(hma[i].hmmname)==inSubSet)
//						newSeq+=currSeq[1].substring(hma[i].tstart-1,hma[i].tend);
//					else
//						newSeq+=nNs(hma[i].tend-hma[i].tstart);
//					pos=hma[i].tend;
//				}
//				newSeq+=currSeq[1].substring(pos);
//				System.out.println(currSeq[0]);
//				System.out.println(newSeq);
//			}
//		}
//	}

	private static String nNs(int n){
		String N="";
		for(int i=0;i<n;i++)
			N+="N";
		return N;
	}
}

class vistaPairsToPlot_pos{
	private String one,two;
	private int newPos;
	
	public vistaPairsToPlot_pos(String one, String two, int newPos) {
		super();
		this.one = one;
		this.two = two;
		this.newPos = newPos;
	}
	
	public String getOne() {
		return one;
	}
	public String getTwo() {
		return two;
	}
	public int getNewPos() {
		return newPos;
	}
	
	
}

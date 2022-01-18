package tools.bam;

import htsjdk.samtools.DefaultSAMRecordFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileHeader;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import tools.fasta.FastaSeq;
import tools.fasta.fastaParser;
import tools.fastq.FastqSeq;
import tools.fastq.fastqParser;
import tools.sequences.sequenceUtils;
import tools.utils.general;

public class bamUtils{
	
	final static String sep= "\t";
	private static HashMap<String, String> methodsHelp= new HashMap<String, String>();
	
	public static void main(String[] args) throws Exception{
		methodsHelp.put("footprinting", "footprinting - Reads bam file, and annotates given positions as small or large given size cutoff\n\targs = <bam file> <bed file> <size cutoff>\n");
		methodsHelp.put("hpSeqPileup", "hpSeqPileup - Generates a raw pileup by replacing the sequence with that in the fastq file, no filtering...\n\targs = <fastq file> <bam file>\n");
		methodsHelp.put("hpSeqPileupCall", "hpSeqPileupCall - Making methylation calls (in CpG context) based on a pileup file from hpSeqPileup\n\targs = <pileup file> \n");
		methodsHelp.put("toFastqRev2", "toFastqRev2 - Printing the reads in fastq format. The second read will be reverse complemented\n\targs = <bam file> \n");
		methodsHelp.put("bamToBiQtab", "bamToBiQtab - Converts a bam alignment to biQtab format. Given an output folder, it outputs references and table. A bed file can be used to specify \n\targs = <bam file> <reference fasta file> <output folder> \n");
		methodsHelp.put("startPoints", "startPoints - prints the start point of each read\n\targs = <bam file> \n");
		methodsHelp.put("hpSeqBamAnalyze", "hpSeqBamAnalyze - Analyze the bam file with wobble sequences = <fastq file> <bam file>\n");
		methodsHelp.put("hpReplaceSeq", "hpReplaceSeq - Replace sequence from bam file with methylated sequence of fastq file= <fastq file> <bam file>\n");
		methodsHelp.put("hpSeqPileup2", "hpSeqPileup2 - Generates a raw pileup, no filtering...\n\targs = <bam file>\n");
		
		if(args.length>0){
			if(args[0].equals("footprinting")&&args.length==3){
				footprinting(args[1],args[2]);
			}else if(args[0].equals("hpSeqPileup")&&args.length==3){
				hpSeqPileup(args[2],args[1]);
			}else if(args[0].equals("hpSeqPileup2")&&args.length==2){
					hpSeqPileup2(args[1]);
			}else if(args[0].equals("hpSeqPileupCall")&&args.length==2){
				hpSeqPileupCall(args[1]);
			}else if(args[0].equals("toFastqRev2")&&args.length==2){
				toFastqRev2(args[1]);
			}else if(args[0].equals("bamToBiQtab")&&args.length==4){
				bamToBiQtab(args[1],args[2],args[3],null);
			}else if(args[0].equals("bamToBiQtab")&&args.length==5){
				bamToBiQtab(args[1],args[2],args[3],args[4]);
			}else if(args[0].equals("startPoints")&&args.length==2){
				startPoints(args[1]);
			}else if(args[0].equals("hpSeqBamAnalyze")&&args.length==3){
				hpSeqBamAnalyze(args[1],args[2]);
			}else if(args[0].equals("hpReplaceSeq")&&args.length==3){
				hpReplaceSeq(args[1],args[2]);
			}else if(args[0].equals("nextMethod")&&args.length==5){

			}else{
				System.err.println(printHelp(args[0]));
				System.exit(616);
			}
		}else{
			System.err.println(printHelp());
			System.exit(616);
		}

	}
	
	private static String printHelp(String cmd){
		if(methodsHelp.containsKey(cmd)){
			return methodsHelp.get(cmd);
		}else{
			return printHelp();
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
	
	private static void startPoints(String bamFileName)throws Exception{
		final SamReader samReader= SamReaderFactory.make()
				.enable(SamReaderFactory.Option.DONT_MEMORY_MAP_INDEX)
				.validationStringency(ValidationStringency.SILENT)
				.samRecordFactory(DefaultSAMRecordFactory.getInstance())
				.open(new File(bamFileName));
		
		for ( final SAMRecord read : samReader){
			if(!read.getNotPrimaryAlignmentFlag() && !read.getReadUnmappedFlag()){
				final String chr= read.getReferenceName();
				final int start= read.getAlignmentStart();
				final int end= read.getAlignmentEnd();
				if(read.getReadNegativeStrandFlag()){
					System.out.println(chr + sep + end + sep + "-");
				}else{
					System.out.println(chr + sep + start + sep + "+");
				}
				 
			}
		}
		
	}

	
	private static void footprinting(String bamFileName, String bedFileName)throws Exception{
		final SamReader samReader= SamReaderFactory.make()
				.enable(SamReaderFactory.Option.DONT_MEMORY_MAP_INDEX)
				.validationStringency(ValidationStringency.SILENT)
				.samRecordFactory(DefaultSAMRecordFactory.getInstance())
				.open(new File(bamFileName));
		
		BufferedReader bedFile= general.getBufferedReader(bedFileName);
		
		BufferedWriter out;
		
		HashMap<String, SortedSet<Integer>> refPos= new HashMap<String, SortedSet<Integer>>();
		SortedSet<Integer> cytosines;
		SortedSet<Integer> overlapped;
		
		for(String s=bedFile.readLine();s!=null;s=bedFile.readLine()){
			final String l[]=s.split("\t");
			if(refPos.containsKey(l[0])){
				cytosines=refPos.get(l[0]);
			}else{
				cytosines= new TreeSet<Integer>();
				refPos.put(l[0], cytosines);
			}
//			does bed files need to be shifted to be compatible with bam files?
			cytosines.add(Integer.parseInt(l[1]));
		}
		
//		TODO: this should be a completely streamed process
		for ( final SAMRecord read : samReader){
			if (!read.getNotPrimaryAlignmentFlag() && read.getMappingQuality()>=20) {
				final String chr = read.getReferenceName();
				if(!refPos.containsKey(chr)){
					refPos.put(chr, new TreeSet<Integer>());
				}
				cytosines= refPos.get(chr);
				overlapped= new TreeSet<Integer>(cytosines.subSet(read.getAlignmentStart(), read.getAlignmentEnd()));
				if (overlapped.size()>0) {
					final String methPattern = (new biQtabLine(read,
							overlapped, "default")).getMethPattern();
					
					char type = methPattern.charAt(0);
					int start = -60000, i = 0;
					System.out.println(chr + sep + read.getAlignmentStart() + sep
							+ read.getAlignmentEnd() + sep + "Read" + sep + read.getReadName());
					for (Integer pos : overlapped) {
						char curType = methPattern.charAt(i);
						i += 1;
						if (type != curType && curType != 'x') {
							if (start > 0) {
								System.out.println(chr + sep + start + sep
										+ (pos - 1) + sep + type + sep + read.getReadName());
							}
							start = pos + 1;
							type = curType;
						}
					}
				}
				
			}
		}
	}
	
	private static void bamToBiQtab(String bamFileName, String faFile, String outputFolder,String bedFileName)throws Exception{
		fastaParser fp= new fastaParser(general.getBufferedReader(faFile));
		FastaSeq fs;
		final String compSuffix="_comp";
		
		final SamReader samReader= SamReaderFactory.make()
				.enable(SamReaderFactory.Option.DONT_MEMORY_MAP_INDEX)
				.validationStringency(ValidationStringency.SILENT)
				.samRecordFactory(DefaultSAMRecordFactory.getInstance())
				.open(new File(bamFileName));
		
		biQtabLine biqLine= new biQtabLine();
		
		
		
		BufferedWriter out;
		
		HashMap<String, SortedSet<Integer>> refPos= new HashMap<String, SortedSet<Integer>>();
		
		int count;
		
		SortedSet<Integer> cytosines;
		
		
		
		
		System.err.println("LOGG: Start reading reference");
		count=0;
		for(;fp.hasNext();++count){
			fs=fp.next();
//			forward strand
			cytosines= new TreeSet<Integer>();
			String qname=fs.getQname();
			String seq=fs.getSeq();
			
			int toAdd= fs.getSeq().indexOf('C');
			while(toAdd!=-1){
				cytosines.add(toAdd);
				toAdd= fs.getSeq().indexOf('C', toAdd+1);
			}
			refPos.put(fs.getQname(), cytosines);
//			Drop forward and reverse Ns (in order to fit into pipeline)
			seq= seq.replaceAll("^N+", "").replaceAll("N+$", "");
			out= new BufferedWriter(new FileWriter(outputFolder+"/"+fs.getQname()+".fasta"));
			out.write(">"+qname+"\n"+seq+"\n");
			out.close();
			
//			complementary strand. No reverse
			qname=fs.getQname()+compSuffix;
			seq=sequenceUtils.Complement(fs.getSeq());
			cytosines= new TreeSet<Integer>();
			toAdd= seq.indexOf('C');
			while(toAdd!=-1){
				cytosines.add(toAdd);
				toAdd= seq.indexOf('C', toAdd+1);
			}
			
			refPos.put(qname, cytosines);
//			Drop forward and reverse Ns (in order to fit into pipeline)
			seq= seq.replaceAll("^N+", "").replaceAll("N+$", "");
			out= new BufferedWriter(new FileWriter(outputFolder+"/"+qname+".fasta"));
			out.write(">"+qname+"\n"+seq+"\n");
			out.close();
		}
		System.err.println("LOGG: Read "+count+" reference files");
		
		
		if(bedFileName!=null){
			System.err.println("LOGG: Reading positions from bedFile. This overrides the fasta file");
			
			BufferedReader bedFile= general.getBufferedReader(bedFileName);
			refPos= new HashMap<String, SortedSet<Integer>>();
			
			for(String s=bedFile.readLine();s!=null;s=bedFile.readLine()){
				final String[] l= s.split("\t");
				
				if(refPos.containsKey(l[0])){
					cytosines= refPos.get(l[0]);
				}else{
					cytosines= new TreeSet<Integer>();
					refPos.put(l[0], cytosines);
				}
				cytosines.add(Integer.parseInt(l[1]));
			}
		}
		
		
		out= new BufferedWriter(new FileWriter(outputFolder+"/BiQtable.tsv"));
		out.write(biqLine.getHeader()+"\n");
		
		System.err.println("LOGG: Parsing BAM file");
		count=0;
		for ( final SAMRecord read : samReader){
			if (!read.getNotPrimaryAlignmentFlag()) {
				
				if (!read.getReadNegativeStrandFlag()) {
					cytosines= refPos.get(read.getReferenceName());
					biqLine= new biQtabLine(read, cytosines, "default");
					out.write(biqLine.toString()+"\n");
				}else{
					cytosines= refPos.get(read.getReferenceName()+compSuffix);
					biqLine= new biQtabLine(read, cytosines, "default");
					biqLine.setRefSeqName(read.getReferenceName()+compSuffix);
					out.write(biqLine.toString()+"\n");
				}
			}
			++count;
		}
		out.close();
		System.err.println("LOGG: Read "+count+" alignments");
	}
	
	private static void toFastqRev2(String bamFileName)throws Exception{
		final SamReader samReader= SamReaderFactory.make()
				.enable(SamReaderFactory.Option.DONT_MEMORY_MAP_INDEX)
				.validationStringency(ValidationStringency.SILENT)
				.samRecordFactory(DefaultSAMRecordFactory.getInstance())
				.open(new File(bamFileName));
		
		
		for ( final SAMRecord read : samReader){
			if (!read.getNotPrimaryAlignmentFlag()) {
				StringBuilder out = new StringBuilder(512);
				out.append("@");
				out.append(read.getReadName());
				out.append("\n");
				if (!read.getReadPairedFlag() || read.getFirstOfPairFlag()) {
					//				read on first strand... just output
					out.append(read.getReadString());
					out.append("\n+\n");
					out.append(read.getBaseQualityString());
				} else {
					//				second in pair... reverse complement
					out.append(sequenceUtils.reverseComplement(read
							.getReadString()));
					out.append("\n+\n");
					out.append(new StringBuilder(read.getBaseQualityString())
							.reverse());
				}
				System.out.println(out.toString());
			}
		}
	}
	
	private static void hpSeqBamAnalyze(String fqFileName, String bamFileName)throws Exception{
		// read bam file
		List<SAMRecord> samFiles= new ArrayList<>();
		Map<Integer,List<SAMRecord>> samByStartPos = new HashMap<>();
		final SamReader samReader= SamReaderFactory.make()
				.enable(SamReaderFactory.Option.DONT_MEMORY_MAP_INDEX)
				.validationStringency(ValidationStringency.SILENT)
				.samRecordFactory(DefaultSAMRecordFactory.getInstance())
				.open(new File(bamFileName));

		//read fastq file 
		int count=0,
			position=0;
		String chr= "",curChr="",seq;
		fastqParser fqp= new fastqParser(general.getBufferedReader(fqFileName),"");
		FastqSeq fq;
		HashMap<String, String> seqStore= new HashMap<String, String>(5500000);
		List<SAMRecord> storage= new ArrayList<>();
		
		System.err.println("Storing sequences...");
		for(;fqp.hasNext();++count){
			if(count%10000==0){
				System.err.print(count+"\r");
			}
			fq=fqp.next();
			seqStore.put(fq.getQname(), fq.getSeq());
		}
		System.err.println(count+"");
		SortedSet<Integer> starts =  new TreeSet<Integer>();
		Map<String,Integer> meta_count = new HashMap<String,Integer>();
		List<String> metakeys = new ArrayList<String>();
		Map<Integer,List<String>> alignStart = new HashMap<Integer,List<String>>();

		for (final SAMRecord samRecord : samReader){
				//replace sequence
				if(seqStore.containsKey(samRecord.getReadName())){
					seq=seqStore.get(samRecord.getReadName());
				}else{
					seq=samRecord.getReadString();
				}
		
				if(samRecord.getReadNegativeStrandFlag()){
					seq=new StringBuilder(seq).reverse().toString();
				}
				if(!samRecord.getReadUnmappedFlag()){
					//for each start group
					
					ArrayList<MyStruct> structlist = new ArrayList<MyStruct>();
					String seqlink = hpSeqLink(seq);
					ArrayList<Integer> seqpos = getCharPosition(seq);
					String metakey = "";
					chr = samRecord.getReferenceName();
					
//					TODO: fix condition to work properly when the chromosome changes
					while(starts.size()>0  && (!chr.equals(curChr)  || (Integer) starts.first() < samRecord.getAlignmentStart()) ){
//						flush and reset data structures
//						final Integer pos = starts.first();
						metakeys= alignStart.get(starts.first());
						for(String mkey : metakeys){
							System.out.println(mkey + "\t" + meta_count.get(mkey));
							meta_count.remove(mkey);
						}
						alignStart.remove(starts.first());
						starts.remove(starts.first());
					}

					curChr=chr;
					double size = seqlink.length()-1;
					for(int i = 0;i<seqpos.size()-1;++i){
						int start = samRecord.getAlignmentStart();
						int pos1=0;
						int pos2=0;
//						TODO: account for indels
						if(samRecord.getReadNegativeStrandFlag()){
							pos1=start-seqpos.get(i);
							pos2=start-seqpos.get(i+1);
						}else{
							pos1=start+seqpos.get(i);
							pos2=start+seqpos.get(i+1);
						}
						char c1=seqlink.charAt(i);
						char c2=seqlink.charAt(i+1);
						metakey=chr+ "\t" +start+ "\t" +pos1+ "\t" +pos2 + "\t" +c1+ "\t" +c2;
						if(meta_count.containsKey(metakey)){
							meta_count.put(metakey, meta_count.get(metakey) + 1);
						}else{
							meta_count.put(metakey, 1);
							if(starts.add(start)){
								//System.out.println(starts.size());
								metakeys= new ArrayList<String>();
								alignStart.put(start,metakeys);
							}else{
								metakeys= alignStart.get(start);
							}
							metakeys.add(metakey);
						}

					}
//					curChr=chr;
//					while((Integer) starts.first() < samRecord.getAlignmentStart() || (!chr.equals(curChr) && starts.size()>0 )){
//					while(starts.size()>0 && (Integer) starts.first() < samRecord.getAlignmentStart() && chr.equals(curChr)){
//						flush and reset data structures
//						final Integer pos = starts.first();
//						metakeys= alignStart.get(starts.first());
//						for(String mkey : metakeys){
//							System.out.println(mkey + "\t" + meta_count.get(mkey));
//							meta_count.remove(mkey);
//						}
//						alignStart.remove(starts.first());
//						starts.remove(starts.first());	
//					}
				}

			}
		}
	
	/**
	 * Replace the sequence from the bam files with the methylated sequence from the fastq files
	 * Also mark duplicates reads by considering the replaced sequence
	 * TODO:   see System.out.print()
	 * No bam output possible, because the lower case letters will be transform by samRecord to capital letters
	 * AND by samtools the non iupac letters will be replaced with 'N'
	 */
	private static void hpReplaceSeq(String fqFileName, String bamFileName)throws Exception{
		final SamReader samReader= SamReaderFactory.make()
				.enable(SamReaderFactory.Option.DONT_MEMORY_MAP_INDEX)
				.validationStringency(ValidationStringency.SILENT)
				.samRecordFactory(DefaultSAMRecordFactory.getInstance())
				.open(new File(bamFileName));

		String chr="",seq,reference;
		int position=0,count=0,duplicate=0;
		int blockreads=0,groupblocks=0,totalblocks=0;
		fastqParser fqp= new fastqParser(general.getBufferedReader(fqFileName),"");
		FastqSeq fq;
		HashMap<String, String> seqStore= new HashMap<String, String>(5500000);
		List<SAMRecord> storage= new ArrayList<>();
		// prepare fastq sequence
		System.err.println("Storing sequences...");
		for(;fqp.hasNext();++count){
			if(count%10000==0){
				System.err.print(count+"\r");
			}
			fq=fqp.next();
			seqStore.put(fq.getQname(), fq.getSeq());
		}
		System.err.println(count+"");
		for (final SAMRecord samRecord : samReader){
			//replace bam sequence with methylated sequence
			if(seqStore.containsKey(samRecord.getReadName())){
				seq=seqStore.get(samRecord.getReadName());
			}else{
				seq=samRecord.getReadString();
			}	
			if(samRecord.getReadNegativeStrandFlag()){
				seq=new StringBuilder(seq).reverse().toString();
			}
			samRecord.setReadString(seq);
			if(!samRecord.getReadUnmappedFlag()){
				reference= samRecord.getReferenceName()+" dna";
				samRecord.setReferenceName(reference);
				//for each start group
				if(!(chr==samRecord.getReferenceName())||!(position==samRecord.getAlignmentStart())){
					//process all read aligning to currentposition
					if(storage.size()>1){
						hpMarkDuplicates(storage);
						blockreads+=storage.size();
						groupblocks++;						
					}
					for(SAMRecord s:storage){
						System.out.print(s.getSAMString());
						if(s.getFlags()==1024||s.getFlags()==1040){
							duplicate++;
						}
					}
					totalblocks++;
					storage.clear();
					storage.add(samRecord);
					chr= samRecord.getReferenceName();
					position= samRecord.getAlignmentStart();
				}else{
					storage.add(samRecord);
				}
			}else{
				System.out.print(samRecord.getSAMString());
			}	
		}
//		System.out.println("duplicate: "+duplicate);
//		System.out.println("Number of reads in grouped blocks: "+blockreads);
//		System.out.println("Number of grouped blocks: "+groupblocks);
//		System.out.println("Number of total blocks: "+totalblocks);
	}
	
	/**
	 * Get SAMRecords blockwise to mark possible duplicates by comparing the CpG positions and the wobble sequence
	 * Wobble sequence is a known sequence in the hairpin linker, which can be seen as index/marker
	 */
	private static void hpMarkDuplicates(List<SAMRecord> samlist){
		int length1=0,length2=0;
		int duplicate=0;
		for(int i=0;i<samlist.size()-1;++i){
			SAMRecord current = samlist.get(i);
			length1= current.getReadString().length();
			String sequence1 = current.getReadString();
			for(int j=i+1;j<samlist.size()-1;++j){
				SAMRecord next= samlist.get(j);
				length2= next.getReadString().length();
				String sequence2 = next.getReadString();
				if(current.getAttribute("ZW")==null || next.getAttribute("ZW")==null){
					if(compareMethylation(sequence1,sequence2)){
						if(current.getFlags()==16){
							if(length1<length2){
								current.setFlags(1040);
							}else{
								next.setFlags(1040);
							}
						}else if(current.getFlags()==0){
							if(length1<length2){
								current.setFlags(1024);
							}else{
								next.setFlags(1024);
							}
						}
					}
				}else if(current.getAttribute("ZW")!=null && next.getAttribute("ZW")!=null){
					String wobble1 = current.getAttribute("ZW").toString();
					String wobble2 = next.getAttribute("ZW").toString();
					if(compareWobble(wobble1,wobble2) && compareMethylation(sequence1,sequence2)){
						if(current.getFlags()==16){
							if(length1<length2){
								current.setFlags(1040);
							}else{
								next.setFlags(1040);
							}
						}else if(current.getFlags()==0){
							if(length1<length2){
								current.setFlags(1024);
							}else{
								next.setFlags(1024);
							}
						}
					}
				}
			}
		}
	}
	
	/**
	 * Compare thw wobble sequence of two reads
	 */
	private static boolean compareWobble(String wobble1, String wobble2){
		boolean iswobdup = false;
		// check only possible if both are not null
		if(wobble1!=null && wobble2!=null){
			if(wobble1==wobble2){
				iswobdup = true;
			}
		}
		// otherwise no check is possible
		else if(wobble1==null || wobble2==null){
			iswobdup = true;
		}
		return iswobdup;
	}
	
	/**
	 * Compare the sequence of two reads only by methylated positions
	 * Trims the sequences to equal length 
	 */
	private static boolean compareMethylation(String sequence1, String sequence2){
		int min1 = sequence1.length();
		int min2 = sequence2.length();
		char c1,c2;
		//trim to equal sequence length
		if(min1<min2){
			sequence2=sequence2.substring(0,min1);
		}else if(min2<min1){
			sequence1=sequence1.substring(0,min2);
		}
		//compare CpG chars 
		for(int i=0;i<sequence1.length()-1;++i){
			c1=sequence1.charAt(i);
			c2=sequence2.charAt(i);
			if(isCpG(c1,c2)){
				if(c1!=c2){
					return false;
				}
			}
		}
		return true;
	}
	
	/**
	 * Ckecks for CpG position
	 */
	private static boolean isCpG(char c1,char c2){
		boolean isCpG = false;
		if(c1=='F'||c1=='P'||c1=='M'||c1=='U'||c1=='f'||c1=='p'||c1=='m'||c1=='u'){
			if(c2=='F'||c2=='P'||c2=='M'||c2=='U'||c2=='f'||c2=='p'||c2=='m'||c2=='u'){
				isCpG=true;
			}
		}
		return isCpG;
	}
	
	private static boolean isCpG(char c1){
		boolean isCpG = false;
		if(c1=='F'||c1=='P'||c1=='M'||c1=='U'){
			isCpG=true;
		}
		return isCpG;
	}
	
	/**
	 * Write linkage statistics
	 * @param seq - string sequence of SAMRecord
	 */
	public static String hpSeqLink(String seq){
		char c = ' ';
		StringBuffer methSeq = new StringBuffer();
		for(int i = 0; i < seq.length();++i){
			c = seq.charAt(i);
			if(isCpG(c)){
				methSeq.append(c);
			}
		}
		return methSeq.toString();
	}
	
	/**
	 * Write statistic for string of all CpG positions
	 * Output will be show like following format:
	 * BASE		BASE	OCCURENCE	FREQUENCY
	 * example input sequence:	UUUUMF
	 * U	U	3	0.6
	 * U	M	1	0.2
	 * M	F	1	0.2
	 * @param seq - string of all CpG positions
	 */
	public static Map<String, Integer> hpMatch(String seq){
		char c1,c2;
		int count =0;
		String substring=new String();
		Map<String,Integer> map = new HashMap<String,Integer>();
		for(int i=0;i<seq.length()-1;++i){
			substring=seq.substring(i, i+2);
			count=countMatches(seq,substring);
			map.put(substring, count);
		}
		return map;
	}
	
	/**
	 * Counts occurency of substring in input string
	 * @param str - input string
	 * @param sub - substring of str of size 2
	 * @return
	 */
	public static int countMatches(String str, String sub) {
		if(isEmpty(str) || isEmpty(sub)) {
			return 0;
			}
		int count = 0;
		int idx = 0;
		while ((idx = str.indexOf(sub, idx)) != -1) {
			count++;
			idx += sub.length()-1;
		}
		return count;
	}
	
	/**
	 * Checks if input string is empty
	 * @param input - input string
	 * @return
	 */
	public static boolean isEmpty(String input){
        if(input != null && input.length() == 0){
            return true;
        }
        return false;
    }
	
	/**
	 * Get positions of char in 
	 * @param str - input string
	 * @param mychar - searched char
	 * @return
	 */
	public static ArrayList<Integer> getCharPosition(String str) {
        ArrayList<Integer> positions = new ArrayList<Integer>();

        if (str.length() == 0)
            return null;

        for (int i = 0; i < str.length(); i ++) {
            if (isCpG(str.charAt(i))) {
                positions.add(i);
            }
        }
        return positions;
	}
	
	public static ArrayList<Integer> getCharPosition(String str, char mychar) {
        ArrayList<Integer> positions = new ArrayList<Integer>();

        if (str.length() == 0)
            return null;

        for (int i = 0; i < str.length(); i ++) {
            if (str.charAt(i) == mychar) {
                positions.add(i);
            }
        }

        return positions;
	}
	
	private static void hpSeqPileupCall(String pileupFileName)throws Exception{
		BufferedReader pileupFile= general.getBufferedReader(pileupFileName);
		hpSeqPileupCounterPair hspcp= new hpSeqPileupCounterPair();
		boolean first=true;
		int pos;
		long count=0;
		String[] l;
		String chr,r1,r2;
		
		for(String s=pileupFile.readLine();s!=null;s=pileupFile.readLine(),++count){
			if(count%10000==0){
				System.err.print(count+"\r");
			}
			l= s.split("\t");
			if(l.length>3){
				chr=l[0];
				pos=Integer.parseInt(l[1]);
				r1=l[2];
				r2=l[3];
				hspcp.propagate(chr, pos, r1, r2);
				
				if(first){
					first=false;
				}else if(hspcp.couldBeMethylated()){
					System.out.println(hspcp);
				}
			}
		}
		++count;
		hspcp.propagate("", -717, "", "");
		if(hspcp.couldBeMethylated()){
			System.out.println(hspcp);
		}
		System.err.println(count);
	}
	
	private static void hpSeqPileup(String bamFileName, String fqFileName)throws Exception{
		final SamReader samReader= SamReaderFactory.make()
				.enable(SamReaderFactory.Option.DONT_MEMORY_MAP_INDEX)
				.validationStringency(ValidationStringency.SILENT)
				.samRecordFactory(DefaultSAMRecordFactory.getInstance())
				.open(new File(bamFileName));
		
//		TODO: add filters
		String curRef= "",seq;
		int curRefIndex= -616,
				refPos,
				count=0;
		
//		is this needed? is there to speed up the flush process
		SortedSet<Integer> activePositions= new TreeSet<Integer>();
		
		HashMap<Integer, StringBuffer> plus= new HashMap<Integer, StringBuffer>(),
				minus= new HashMap<Integer, StringBuffer>(),
				current;
		
		fastqParser fqp= new fastqParser(general.getBufferedReader(fqFileName),"");
		FastqSeq fq;
		HashMap<String, String> seqStore= new HashMap<String, String>(5500000);
		
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
		
		
		for (final SAMRecord samRecord : samReader){
			if(!samRecord.getReadUnmappedFlag()){
				if(samRecord.getReferenceIndex() != curRefIndex){
//					flush data and reset structures
					for(;activePositions.size()>0;){
						final Integer pos=activePositions.first();
						final String p=plus.containsKey(pos)?plus.remove(pos).toString():"",
								m=minus.containsKey(pos)?minus.remove(pos).toString():"";
						System.out.println(curRef+sep+activePositions.first()+sep+p+sep+m);
						activePositions.remove(pos);
					}
					
					plus.clear();
					minus.clear();
					
					curRef=samRecord.getReferenceName();
					curRefIndex=samRecord.getReferenceIndex();
				}

				for(;activePositions.size()>0 && activePositions.first()<samRecord.getStart();){
					final Integer pos=activePositions.first();
					final String p=plus.containsKey(pos)?plus.remove(pos).toString():"",
							m=minus.containsKey(pos)?minus.remove(pos).toString():"";
					System.out.println(curRef+sep+activePositions.first()+sep+p+sep+m);
					activePositions.remove(pos);
				}
				
				
				current=samRecord.getReadNegativeStrandFlag()?minus:plus;
				
				if(seqStore.containsKey(samRecord.getReadName())){
					seq=seqStore.get(samRecord.getReadName());
				}else{
					seq=samRecord.getReadString();
				}
				
				if(samRecord.getReadNegativeStrandFlag()){
					seq=new StringBuilder(seq).reverse().toString();
				}
				

				for(int i=0;i<seq.length();++i){
					if((refPos= samRecord.getReferencePositionAtReadPosition(i+1))!=0){
						activePositions.add(refPos);
						if(!current.containsKey(refPos)){
							current.put(new Integer(refPos), new StringBuffer(""));
						}
						current.get(refPos).append(seq.charAt(i));
					}
				}
			}
		}
		
		for(;activePositions.size()>0;){
			final Integer pos=activePositions.first();
			final String p=plus.containsKey(pos)?plus.remove(pos).toString():"",
					m=minus.containsKey(pos)?minus.remove(pos).toString():"";
			System.out.println(curRef+sep+activePositions.first()+sep+p+sep+m);
			activePositions.remove(pos);
		}
		
		
	}
	private static void hpSeqPileup2(String bamFileName)throws Exception{
		final SamReader samReader= SamReaderFactory.make()
				.enable(SamReaderFactory.Option.DONT_MEMORY_MAP_INDEX)
				.validationStringency(ValidationStringency.SILENT)
				.samRecordFactory(DefaultSAMRecordFactory.getInstance())
				.open(new File(bamFileName));
		
//		TODO: add filters
		String curRef= "",seq;
		int curRefIndex= -616,
				refPos,
				count=0;
		
//		is this needed? is there to speed up the flush process
		SortedSet<Integer> activePositions= new TreeSet<Integer>();
		
		HashMap<Integer, StringBuffer> plus= new HashMap<Integer, StringBuffer>(),
				minus= new HashMap<Integer, StringBuffer>(),
				current;
	
		for (final SAMRecord samRecord : samReader){
			if(!samRecord.getReadUnmappedFlag()){
				if(samRecord.getReferenceIndex() != curRefIndex){
//					flush data and reset structures
					for(;activePositions.size()>0;){
						final Integer pos=activePositions.first();
						final String p=plus.containsKey(pos)?plus.remove(pos).toString():"",
								m=minus.containsKey(pos)?minus.remove(pos).toString():"";
						System.out.println(curRef+sep+activePositions.first()+sep+p+sep+m);
						activePositions.remove(pos);
					}
					
					plus.clear();
					minus.clear();
					
					curRef=samRecord.getReferenceName();
					curRefIndex=samRecord.getReferenceIndex();
				}

				for(;activePositions.size()>0 && activePositions.first()<samRecord.getStart();){
					final Integer pos=activePositions.first();
					final String p=plus.containsKey(pos)?plus.remove(pos).toString():"",
							m=minus.containsKey(pos)?minus.remove(pos).toString():"";
					System.out.println(curRef+sep+activePositions.first()+sep+p+sep+m);
					activePositions.remove(pos);
				}
				
				current=samRecord.getReadNegativeStrandFlag()?minus:plus;
				
				for(int i=0;i<samRecord.getReadString().length();++i){
					if((refPos= samRecord.getReferencePositionAtReadPosition(i+1))!=0){
						activePositions.add(refPos);
						if(!current.containsKey(refPos)){
							current.put(new Integer(refPos), new StringBuffer(""));
						}
						current.get(refPos).append(samRecord.getReadString().charAt(i));
					}
				}
			}
		}
		for(;activePositions.size()>0;){
			final Integer pos=activePositions.first();
			final String p=plus.containsKey(pos)?plus.remove(pos).toString():"",
					m=minus.containsKey(pos)?minus.remove(pos).toString():"";
			System.out.println(curRef+sep+activePositions.first()+sep+p+sep+m);
			activePositions.remove(pos);
		}
	}
}

class hpSeqPileupCounter{
	
	private char consensus;
	private int meth, unmeth, cov, pos, 
	F1, F2, P1, P2, M1, M2, U1, U2, 
	f1, f2, p1, p2, m1, m2, u1, u2,
	countA, countT, countC, countG, countN;
	private String chr;
	
	
	protected hpSeqPileupCounter(String chr, int pos, String r1, String r2){
//		all this counting back and forth could surely be optimized...
		
		
		meth=unmeth=cov=F1=F2=P1=P2=M1=M2=U1=U2=f1=f2=p1=p2=m1=m2=u1=u2=countA=countT=countC=countG=countN=0;
		
		this.pos=pos;
		this.chr=chr;
		
		addSeqs(r1,r2);	
	}
	
	protected hpSeqPileupCounter(){
		this("",-616,"", "");
	}
	
	private void addSeqs(String r1, String r2){
//		CGfull('F'), CGplus('P'), CGminus('M'), CGnone('U'),
//		GCfull('f'), GCplus('p'), GCminus('m'), GCnone('u'),
//		Cmeth('C'), Cnone('Y'), CnoneLast('y'),
//		Gmeth('G'), Gnone('R'), GnoneLast('r'),
//		Anone('A'), Tnone('T'),
//		unknown('N');
		
		cov+=r1.length()+r2.length();
		
		for(int i=0, n=r1.length(); i<n; ++i){
			final char c= r1.charAt(i);
			
			switch (c) {
			case 'A':
				++countA;
				break;
			case 'T':
				++countT;
				break;
			case 'Y':
			case 'y':
				++countC;
				++unmeth;
				break;
			case 'R':
			case 'r':
				++countG;
				++unmeth;
				break;
			case 'C':
				++countC;
				++meth;
				break;
			case 'G':
				++countG;
				++meth;
				break;
			case 'f':
				++countG;
				++meth;
				++f1;
				break;
			case 'u':
				++countG;
				++unmeth;
				++u1;
				break;
			case 'F':
				++countC;
				++meth;
				++F1;
				break;
			case 'U':
				++countC;
				++unmeth;
				++U1;
				break;
			case 'p':
				++countG;
				++meth;
				++p1;
				break;
			case 'm':
				++countG;
				++unmeth;
				++m1;
				break;
			case 'P':
				++countC;
				++meth;
				++P1;
				break;
			case 'M':
				++countC;
				++unmeth;
				++M1;
				break;
			case 'N':
				++countN;
				break;
			default:
				break;
			}
		}
		
		for(int i=0, n=r2.length(); i<n; ++i){
			final char c= r2.charAt(i);
			
			switch (c) {
			case 'A':
				++countT;
				break;
			case 'T':
				++countA;
				break;
			case 'Y':
			case 'y':
				++countG;
				++unmeth;
				break;
			case 'R':
			case 'r':
				++countC;
				++unmeth;
				break;
			case 'C':
				++countG;
				++meth;
				break;
			case 'G':
				++countC;
				++meth;
				break;
			case 'f':
				++countC;
				++meth;
				++f2;
				break;
			case 'u':
				++countC;
				++unmeth;
				++u2;
				break;
			case 'F':
				++countG;
				++meth;
				++F2;
				break;
			case 'U':
				++countG;
				++unmeth;
				++U2;
				break;
			case 'p':
				++countC;
				++meth;
				++p2;
				break;
			case 'm':
				++countC;
				++unmeth;
				++m2;
				break;
			case 'P':
				++countG;
				++meth;
				++P2;
				break;
			case 'M':
				++countG;
				++unmeth;
				++M2;
				break;
			case 'N':
				++countN;
				break;
			default:
				break;
			}
		}
		consensus= calculateConsensus();
	}

	
	private char calculateConsensus(){
		int max=countA,
				sum=countA;
		char nuc='A';
		
		if(max>cov-sum){
			return nuc;
		}
		
		if(countT>max){
			max=countT;
			nuc='T';
		}else if(countT==max){
			nuc='N';
		}
		sum+=countT;
		
		if(max>cov-sum){
			return nuc;
		}
		
		
		if(countC>max){
			max=countC;
			nuc='C';
		}else if(countC==max){
			nuc='N';
		}
		sum+=countC;
		
		if(max>cov-sum){
			return nuc;
		}

		if(countG>max){
			max=countG;
			nuc='G';
		}else if(countG==max){
			nuc='N';
		}
		sum+=countG;
		
		if(max>cov-sum){
			return nuc;
		}
		
		if(countN>=max){
			return 'N';
		}
		
		return nuc;
	}

	
	protected int getCGfull(){
		return F1;
	}
	
	protected int getCGplus(){
		return P1;
	}
	
	protected int getCGminus(){
		return M1;
	}
	
	protected int getCGnone(){
		return U1;
	}
	
	protected int getGCfull(){
		return f1;
	}
	
	protected int getGCplus(){
		return p1;
	}
	
	protected int getGCminus(){
		return m1;
	}
	
	protected int getGCnone(){
		return u1;
	}
	
	protected int getCGfullRev(){
		return F2;
	}
	
	protected int getCGplusRev(){
		return P2;
	}
	
	protected int getCGminusRev(){
		return M2;
	}
	
	protected int getCGnoneRev(){
		return U2;
	}
	
	protected int getGCfullRev(){
		return f2;
	}
	
	protected int getGCplusRev(){
		return p2;
	}
	
	protected int getGCminusRev(){
		return m2;
	}
	
	protected int getGCnoneRev(){
		return u2;
	}
	
	
	
	
	protected char getConsensus(){
		return consensus;
	}
	
	protected int getMeth() {
		return meth;
	}

	protected int getUnmeth() {
		return unmeth;
	}

	protected int getCov() {
		return cov;
	}

	protected int getPos() {
		return pos;
	}

	protected String getChr() {
		return chr;
	}
	
	
}

class hpSeqPileupCounterPair{
	private final static String sep="\t"; 
	
	private enum Context{
		GC,CG,other,unknown;
	}
	
	private hpSeqPileupCounter p1;
	private hpSeqPileupCounter p2;
	private Context context;
	
	
	
	protected hpSeqPileupCounterPair(String chr, int pos, String r1, String r2){
		p2= new hpSeqPileupCounter();
		this.propagate(chr, pos, r1, r2);
		
	}
	
	protected hpSeqPileupCounterPair(String s){
		this(s.split("\t")[0],Integer.parseInt(s.split("\t")[1]),s.split("\t")[2], s.split("\t")[3]);
	}
	
	protected hpSeqPileupCounterPair(){
		this("",-616,"","");
	}
	
	protected void propagate(String chr, int pos, String r1, String r2){
		p1=p2;
		p2= new hpSeqPileupCounter(chr, pos, r1, r2);
		
		
		if(p1.getPos()+1==p2.getPos() && p1.getChr().equals(p2.getChr())){
			if(p1.getConsensus()=='G' && p2.getConsensus()=='C'){
				context=Context.GC;
			}else if(p1.getConsensus()=='C' && p2.getConsensus()=='G'){
				context=Context.CG;
			}else{
				context=Context.other;
			}
		}else{
			context=Context.unknown;
		}
	}
	
	protected boolean couldBeMethylated(){
		return p1.getConsensus()=='C' || p1.getConsensus()=='G';
	}
	
	public String toString(){
		StringBuilder res= new StringBuilder(p1.getChr());
		res.append(sep).append(p1.getPos()).append(sep).append(p1.getConsensus());
		res.append(sep).append(p1.getCov());
		
		switch (p1.getConsensus()) {
		case 'C':
		case 'G':
			res.append(sep).append(p1.getMeth()).append(sep).append(p1.getUnmeth());
			break;
		default:
			res.append(sep).append("NA").append(sep).append("NA");
			break;
		}
		
		

		switch (context) {
		case other:
			//				add context
			res.append(sep).append(Context.other);
			//				add symmetric methylation
			//				Full
			res.append(sep).append("NA");
			//				Plus
			res.append(sep).append("NA");
			//				Minus
			res.append(sep).append("NA");
			//				None
			res.append(sep).append("NA");
			break;
		case CG:
			//				add context
			res.append(sep).append(Context.CG);
			//				add symmetric methylation
			//				Full
			res.append(sep).append(p1.getCGfull()+p2.getCGfullRev());
			//				Plus
			res.append(sep).append(p1.getCGplus()+p2.getCGminusRev());
			//				Minus
			res.append(sep).append(p1.getCGminus()+p2.getCGplusRev());
			//				None
			res.append(sep).append(p1.getCGnone()+p2.getCGnoneRev());
			break;
		case GC:
			//				add context
			res.append(sep).append(Context.GC);
			//				add symmetric methylation
			//				Full
			res.append(sep).append(p1.getGCfull()+p2.getGCfullRev());
			//				Plus
			res.append(sep).append(p1.getGCplus()+p2.getGCminusRev());
			//				Minus
			res.append(sep).append(p1.getGCminus()+p2.getGCplusRev());
			//				None
			res.append(sep).append(p1.getGCnone()+p2.getGCnoneRev());
			break;	
		case unknown:
		default:
			//				add context
			res.append(sep).append(Context.unknown);
			//				add symmetric methylation
			//				Full
			res.append(sep).append("NA");
			//				Plus
			res.append(sep).append("NA");
			//				Minus
			res.append(sep).append("NA");
			//				None
			res.append(sep).append("NA");
			break;
		}


		
		return res.toString();
	}
}


class biQtabLine{
	
	private String 	ID,
					methPattern,
					refSeqName,
					sampleName="default";
	
	final private String header= "ID\tAlignment score\tSequence identity\tMethylation pattern\tMean methylation\tMissing sites\tConversion rate\tReference sequence\tSample name";
	
	private Integer	score,
					missingSites;
	
	private Double	seqIdentity, // could be done together with reference sequence 
					meanMeth,
					conversionRate;
	protected biQtabLine(SAMRecord read, SortedSet<Integer> cytosines,String sampleName){
	
		this();
		this.ID= read.getReadName();
		this.score= read.getMappingQuality();
		this.refSeqName= read.getReferenceName();
		this.sampleName= sampleName;
		
//		initialize seqIdentity, methPattern
		
		StringBuilder methP= new StringBuilder(512);
		
		
		final boolean fwd=!read.getReadNegativeStrandFlag();
		
		final boolean ct=(!read.getReadPairedFlag() || read.getFirstOfPairFlag());
		
		
		final char 	meth=fwd?'C':'G',
					unmeth=fwd?'T':'A';
		
		
		
		
		
		final String seq=read.getReadString();
		int readPos= 0;
		int	refPos= read.getReferencePositionAtReadPosition(readPos+1) - 1;
		missingSites=refPos<0?-1:0;
		
		for (Integer pos : cytosines){
			if(pos < refPos || pos > read.getAlignmentEnd() - 1){
				methP.append('x');
			}else{
				while(readPos<seq.length() && refPos < pos){
					++readPos;
					refPos= read.getReferencePositionAtReadPosition(readPos+1) - 1;
				}
				if(refPos>pos){
					++missingSites;
					methP.append('x');
				}else{
					final char nuc= seq.charAt(readPos);
					
					if(nuc == unmeth){
						methP.append('0');
					}else if( nuc == meth){
						methP.append('1');
					}else{
						++missingSites;
						methP.append('x');
					}
				}

			}
		}
		
		
		
		
		
		
		
		assert methP.length()!=cytosines.size();
		
		this.methPattern= methP.toString();
	}
	
	protected biQtabLine(){
//		set default values
		ID="readName";
		score=-1;
		seqIdentity=-1.0;
		methPattern="1x0";
		meanMeth=-1.0;
		missingSites=-1;
		conversionRate=-1.0;
		refSeqName="reference";
		sampleName="defaultSample";
	}
	
	public String toString(){
		return String.format("%s\t%s\t%1.4f\t%s\t%1.4f\t%s\t%1.4f\t%s\t%s", ID, 
				score, seqIdentity, methPattern, meanMeth, missingSites,
				conversionRate, refSeqName, sampleName);
	}


	protected String getRefSeqName() {
		return refSeqName;
	}




	protected void setRefSeqName(String refSeqName) {
		this.refSeqName = refSeqName;
	}




	protected String getSampleName() {
		return sampleName;
	}




	protected void setSampleName(String sampleName) {
		this.sampleName = sampleName;
	}




	protected String getHeader() {
		return header;
	}
	
	protected String getMethPattern() {
		return methPattern;
	}
}

class MyStruct{
	protected int pos;
	protected char c;
	
	public MyStruct(){
		this.pos = 0;
		this.c = ' ';
	}
	
	public int getPosition(){
		return pos;
	}
	
	public void setPosition(int pos){
		this.pos=pos;
	}
	
	public char getChar(){
		return c;
	}
	
	public void setChar(char c){
		this.c=c;
	}
}
































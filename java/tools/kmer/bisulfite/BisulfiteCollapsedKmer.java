package tools.kmer.bisulfite;

import java.util.ArrayList;
import java.util.BitSet;

import tools.kmer.KmerSet_binary_utils;


public class BisulfiteCollapsedKmer {

	private BitSet kmer;
	private char methylated,unmethylated;
	private Integer[] counts;
	private int referred,occurences;
	
	protected BisulfiteCollapsedKmer(BitSet kmer, String orig, char methylated, char unmethylated, boolean isTarget){
		this.kmer=kmer;
		this.methylated=methylated;
		this.unmethylated=unmethylated;
		ArrayList<Integer>countsTmp=new ArrayList<Integer>();
		if(isTarget){
			referred=1;
			occurences=0;
			for(char nuc : orig.toCharArray()){
				if(nuc==methylated || nuc==unmethylated){
					countsTmp.add(0);
				}
			}
		}else{
			referred=0;
			occurences=1;
			for(char nuc : orig.toCharArray()){
				if(nuc==methylated){
					countsTmp.add(1);
				}else if(nuc==unmethylated){
					countsTmp.add(0);
				}
			}	
		}
		counts=countsTmp.toArray(counts);
//		System.out.println(orig+"\t"+counts.size()+"\t"+methylated);
	}
	
//	protected BisulfiteKmer(String orig, char methylated, char unmethylated,KmerSet_binary_utils kt_utils, boolean isTarget){
//		this(kt_utils.stringToBitSet(orig.replace(methylated, unmethylated)),orig,methylated,unmethylated,isTarget);
//	}
	
	protected boolean addKmer(String orig)throws Exception{
		for(int i=0,c=0;i<orig.length();++i){
			final char nuc=orig.charAt(i);
			if(nuc==methylated){
				if(c>=counts.length){
					KmerSet_binary_utils kt_utils= new KmerSet_binary_utils(orig.length());
					throw new Exception(orig+"\t"+kt_utils.bitSetToString(kmer)+"\t"+counts.length+"\t"+c+"\t"+kmer.hashCode()+"\t"+kmer.toString());
				}
				counts[c]++;
				++c;
			}else if(nuc==unmethylated){
				++c;
			}
		}
		++occurences;
		return true;
	}
	
	protected void increaseReferred(){
		++referred;
	}
	
	public int hashCode(){
		return kmer.hashCode();
	}

	protected int getReferred() {
		return referred;
	}

	protected int getOccurences() {
		return occurences;
	}

	protected Integer getCounts(int pos) {
		return counts[pos];
	}
	
	
}

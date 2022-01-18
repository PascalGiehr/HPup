package tools.kmer;

import java.io.Serializable;
import java.util.BitSet;
import java.util.Collection;

import tools.sequences.sequenceUtils;

public class KmerSet_binary_utils implements Serializable{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private int kmerSize;
	
	public KmerSet_binary_utils(int kmerSize) {
		this.kmerSize=kmerSize;
	}
	
	protected int getKmerSize(){
		return kmerSize;
	}
	
	public void kmerize(String seq, Collection<BitSet> out){
		BitSet bs;
		if(this.getKmerSize()<=seq.length()){
			bs=this.stringToBitSet(seq.substring(0, this.getKmerSize()));
			out.add(this.kmerToUse(bs));
			for(int i=this.getKmerSize();i<seq.length();++i){
				bs= this.shiftBitSet(bs, seq.charAt(i));
				out.add(this.kmerToUse(bs));
			}
		}
	}
	
	public void kmerizeNoN(String seq, Collection<BitSet> out){
		String[] l= seq.split("N+");
		for(String s: l){
			this.kmerize(s, out);
//			if(this.getKmerSize()<=s.length()){
//				bs=this.stringToBitSet(s.substring(0, this.getKmerSize()));
//				out.add(this.kmerToUse(bs));
//				for(int i=this.getKmerSize();i<s.length();++i){
//					bs= this.shiftBitSet(bs, s.charAt(i));
//					out.add(this.kmerToUse(bs));
//				}
//			}
		}
	}
	
	public void kmerizeNoNForward(String seq, Collection<BitSet> out){
		String[] l= seq.split("N+");
		BitSet bs;
		for(String s: l){
			if(this.getKmerSize()<=s.length()){
				bs=this.stringToBitSet(s.substring(0, this.getKmerSize()));
				out.add(bs);
				for(int i=this.getKmerSize();i<s.length();++i){
					bs= this.shiftBitSet(bs, s.charAt(i));
					out.add(bs);
				}
			}
		}
	}
	
	public BitSet shiftBitSet(BitSet oldKmer, char nuc){
		//what goes wrong in the generateion of kmers
		BitSet kmer= oldKmer.get(3, kmerSize*3);
//		for(int i=0,j=3;j<kmerSize*3;i++,j++){
//			if(oldKmer.get(j)){
//				kmer.set(i);
//			}
//		}
		setPositionsInBitSet(kmer, (kmerSize-1)*3, nuc);
		kmer.set(kmerSize*3);
		return kmer;
	}
	
	public BitSet stringToBitSet(String kmer_seq){
		BitSet kmer= new BitSet(kmerSize*3+1);
		for(int i=0,j=0;i<kmerSize;i++){
			j=setPositionsInBitSet(kmer, j, kmer_seq.charAt(i));
		}
		kmer.set(kmerSize*3);
		return kmer;
	}
	
	protected BitSet kmerToUse(BitSet kmer){
//		System.out.println(bitSetToString(kmer));
		BitSet kmerRevComp= new BitSet(kmerSize*3+1);
		int j=(kmerSize-1)*3,i=0,result;
		char c,r;
		for(;j>=0;){
			c=getPositionInBitSet(kmer, i);
			r=sequenceUtils.complement(getPositionInBitSet(kmer, j));
//			System.out.println(c+"\t"+r);
			if((result=c-r)!=0){
				if(result<0){
					return (BitSet) kmer.clone();
				}else{
					i=setPositionsInBitSet(kmerRevComp, i, r);
					--j;/**
					if(!kmer.get(--j)){
						kmerRevComp.set(i++);
					}
					if(!kmer.get(--j)){
						kmerRevComp.set(i++);
					}
					if(!kmer.get(--j)){
						kmerRevComp.set(i++);
					}*/
					break;
				}
			}else{
				i=setPositionsInBitSet(kmerRevComp, i, r);
				j-=3;
				/**
				if(!kmer.get(--j)){
					kmerRevComp.set(i++);
				}
				if(!kmer.get(--j)){
					kmerRevComp.set(i++);
				}
				if(!kmer.get(--j)){
					kmerRevComp.set(i++);
				}*/
			}
			/**
			if(kmer.get(i)){
				if(kmer.get(j)){
					//rev comp smaller
					++i;
					--j;
					break;
				}else{
					kmerRevComp.set(i);
					//continue
				}
			}else{
				if(kmer.get(j)){
					//continue
					
				}else{
					// original smaller
					return kmer;
				}
			}*/
		}
		for(;j>=0;++i,--j){
			if(!kmer.get(j)){
				kmerRevComp.set(i);
			}
		}
		kmerRevComp.set(kmerSize*3);
//		System.out.println(bitSetToString(kmerRevComp));
		return kmerRevComp;
	}
	
	protected BitSet kmerToUse(String kmer_string){
		BitSet kmer= new BitSet(kmerSize*3+1);
		kmer.set(kmerSize*3);
		BitSet kmerRevComp= new BitSet(kmerSize*3+1);
		kmerRevComp.set(kmerSize*3);
		int result,j=kmerSize-1,ic=0;
		char c,r;
		for(int i=0;j>=0;++i,--j){
			c=kmer_string.charAt(i);
			r=sequenceUtils.complement(kmer_string.charAt(j));
			if((result=c-r)!=0){
				if(result<0){
					//construct the full orignal bitkmer and return it
					ic=setPositionsInBitSet(kmer, ic, c);
					++i;
					for(;i<kmerSize;++i){
						ic=setPositionsInBitSet(kmer, ic, kmer_string.charAt(i));
					}
					return kmer;
				}else{
					//construct the full rev.comp bitkmer and return it
					ic=setPositionsInBitSet(kmerRevComp, ic, r);
					--j;
					break;
				}
			}else{
				setPositionsInBitSet(kmer, ic, c);
				ic=setPositionsInBitSet(kmerRevComp, ic, r);
			}
		}
		for(;j>=0;--j){
			ic=setPositionsInBitSet(kmerRevComp, ic, sequenceUtils.complement(kmer_string.charAt(j)));
		}
		return kmerRevComp;
	}

	private int setPositionsInBitSet(BitSet kmer,int pos,char nuc){
		//assumes that the bits are not set
		switch (nuc) {
		case 'A':
		case 'a':
			//000
			pos+=3;
			break;
		case 'C':
		case 'c':
			//001
			pos+=2;
			kmer.set(pos++);
			break;
		case 'G':
		case 'g':
			//011
			pos++;
			kmer.set(pos++);
			kmer.set(pos++);
			break;
		case 'N':
		case 'n':
			//100
			kmer.set(pos++);
			pos++;
			pos++;
			break;
		case 'T':
		case 't':
			//111
			kmer.set(pos++);
			kmer.set(pos++);
			kmer.set(pos++);
			break;
		default:
			//110
			System.err.println("unknown char (shift) "+nuc);
			kmer.set(pos++);
			kmer.set(pos++);
			pos++;
			break;
		}
		return pos;
	}
	
	public String bitSetToString(BitSet kmer){
		StringBuffer s=new StringBuffer();
		for(int i=0;i<kmerSize*3;i+=3){
			s.append(getPositionInBitSet(kmer, i));
		}
		return s.toString();
	}
	
	private char getPositionInBitSet(BitSet kmer,int position){
		if(kmer.get(position)){
			if(kmer.get(position+1)){
				if(kmer.get(position+2)){
					//111
					return 'T';
				}else{
					//110
					return 'N';
				}
			}else{
				if(kmer.get(position+2)){
					//101
					System.err.println("# KmerSet_binary_util unknown coding: 101");
					return 'X';
				}else{
					//100
					return 'N';
				}
			}
		}else{
			if(kmer.get(position+1)){
				if(kmer.get(position+2)){
					//011
					return 'G';
				}else{
					//010
					System.err.println("# KmerSet_binary_util unknown coding: 010");
					return 'X';
				}
			}else{
				if(kmer.get(position+2)){
					//001
					return 'C';
				}else{
					//000
					return 'A';
				}
			}
		}
	}
}

package tools.kmer;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.util.zip.GZIPInputStream;

import tools.sequences.sequenceUtils;

public class kmerFunctions {
	
	public static KmerSet_binary readKmers(String kmerCountFile, int cutoff)throws Exception{
		BufferedReader in= new BufferedReader(new FileReader(kmerCountFile));
		return readKmers(in, cutoff);
	}
	
	public static KmerSet_binary readKmersZ(String kmerCountFile, int cutoff)throws Exception{
		BufferedReader in= new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(kmerCountFile))));
		return readKmers(in, cutoff);
	}
	
	public static KmerSet_binary readKmers(BufferedReader in, int cutoff)throws Exception{
		
		String[] l=in.readLine().split("\t");
		KmerSet_binary kmers= new KmerSet_binary(l[1].length());
		if(Integer.parseInt(l[0])>=cutoff){
			kmers.addSeq(l[1]);
		}
		int line=1;
		for(String s=in.readLine();s!=null;s=in.readLine(),line++){
			if(line%10000==0)
				System.err.print("\t"+line+"   "+kmers.size()+"\r");
			l=s.split("\t");
			if(Integer.parseInt(l[0])>=cutoff){
				kmers.addSeq(l[1]);
			}
		}
		System.err.println("\t"+line+"   "+kmers.size());
		
		return kmers;
	}

	public static KmerSet_binary readKmers(String kmerFile)throws Exception{
		return (KmerSet_binary)(new ObjectInputStream(new FileInputStream(kmerFile))).readObject();
	}

	
	public static ArrayList<String> kmerToUse(StringBuffer seq, int kmerSize){
		ArrayList<String> kmers= new ArrayList<String>();
		if(seq.length()>=kmerSize){
			StringBuffer kmer=new StringBuffer(seq.substring(0, kmerSize));
			kmers.add(kmerToUse(kmer));
			for(int i=kmerSize;i<seq.length();++i){
				kmer.deleteCharAt(0);
				kmer.append(seq.charAt(i));
				kmers.add(kmerToUse(kmer));
			}
		}
		return kmers;
	}
	
	public static ArrayList<String> kmers(StringBuffer seq, int kmerSize){
		ArrayList<String> kmers= new ArrayList<String>();
		if(seq.length()>=kmerSize){
			StringBuffer kmer= new StringBuffer(seq.substring(0, kmerSize));
			kmers.add(kmer.toString());
			for(int i=kmerSize;i<seq.length();++i){
				kmer.deleteCharAt(0);
				kmer.append(seq.charAt(i));
				kmers.add(kmer.toString());
			}
		}
		return kmers;
	}
	
	public static ArrayList<String> kmers(String seq, int kmerSize){
		return kmers(new StringBuffer(seq),kmerSize);
	}
	
	public static ArrayList<String> kmerToUse(String seq, int kmerSize){
		return kmerToUse(new StringBuffer(seq), kmerSize);
	}

	public static String kmerToUse(String kmer){
		return kmerToUse(new StringBuffer(kmer));
	}
	
	public static String kmerToUse(StringBuffer kmer){
		StringBuffer kmerRevComp= new StringBuffer();
		int j=kmer.length()-1,result;
		char c,r;
		for(int i=0;i<kmer.length();i++,j--){
			c=kmer.charAt(i);
			r=sequenceUtils.complement(kmer.charAt(j));
			if((result=c-r)!=0){
				if(result<0){
					return kmer.toString();
				}else{
					kmerRevComp.append(r);
					j--;
					break;
				}
			}else{
				kmerRevComp.append(r);
			}
		}
		for(;j>=0;j--){
			kmerRevComp.append(sequenceUtils.complement(kmer.charAt(j)));
		}
		return kmerRevComp.toString();
	}
	
	
}

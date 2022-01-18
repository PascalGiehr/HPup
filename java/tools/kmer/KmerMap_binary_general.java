package tools.kmer;

import java.io.Serializable;
import java.util.BitSet;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;

public class KmerMap_binary_general<T> extends HashMap<BitSet, T> implements Serializable{

	private static final long serialVersionUID = 1L;
	private int kmerSize;
	private KmerSet_binary_utils kt_utils;
//	private HashMap<BitSet, T> Buffer;
//	private final int bufferSize=20000000;
	
	public KmerMap_binary_general(int kmerSize){
		super(600000000);
		this.kmerSize=kmerSize;
		kt_utils= new KmerSet_binary_utils(kmerSize);
//		Buffer= new HashMap<BitSet, T>(bufferSize*2);
	}
	
	
	
//	public int addSequenceBuffered(String seq){
//		BitSet kmer=kt_utils.stringToBitSet(seq.substring(0, kmerSize));
//		
////		System.err.println(kt_utils.bitSetToString(kmer));
////		System.err.println(kt_utils.bitSetToString(kt_utils.kmerToUse(kmer))+" O");
//		
//		if(Buffer.containsKey(kmer)){
//			Buffer.put(kmer, Buffer.get(kmer)+1);
//		}else{
//			Buffer.put(kmer, 1);
//		}
//		
//		for(int i=kmerSize;i<seq.length();i++){
//			kmer=kt_utils.shiftBitSet(kmer, seq.charAt(i));
////			System.err.println(kt_utils.bitSetToString(kmer));
////			System.err.println(kt_utils.bitSetToString(kt_utils.kmerToUse(kmer))+" O");
//			
//			if(Buffer.containsKey(kmer)){
//				Buffer.put(kmer, Buffer.get(kmer)+1);
//			}else{
//				Buffer.put(kmer, 1);
//			}
//		}
//		if(Buffer.size()>bufferSize){
//			System.err.print("# flushing ...");
//			this.flushBuffer();
//			System.err.println(" Done");
//		}
//		return seq.length()-kmerSize+1;
//	}
	
	
//	public void flushBuffer(){
//		this.addBuffer(Buffer);
//		Buffer=new HashMap<BitSet, T>(bufferSize*2);
//	}
	
	protected void addBuffer(HashMap<BitSet, T> m2){
		Iterator<Entry<BitSet, T>> it=m2.entrySet().iterator();
		while(it.hasNext()){
			Entry<BitSet, T> e= it.next();
			//questionable whether the kmerToUse function should be here.
			//For the buffer it is definitely the best place, but it might
			//alter another kmer set..
			this.addKmer(kt_utils.kmerToUse(e.getKey()), e.getValue());
		}
	}
	
//	public boolean addKmer(BitSet kmer){
//		return addKmer(kmer, 1); 
//	}
	
	public boolean addKmer(String kmer,T count){
		if(kmer.length()!=kmerSize){
			System.err.println("Tried to add a kmer: "+kmer+", which isn't "+ kmerSize+ "bp");
			return false;
		}else{
			this.addKmer(kt_utils.stringToBitSet(kmer), count);
			return true;
		}
	}
	
	public T get(String kmer_string){
		if(kmer_string.length()!=kmerSize){
			System.err.println("Tried to find a kmer with wrong size: "+kmer_string);
			return null;
		}
		BitSet kmer=kt_utils.kmerToUse(kmer_string);
		if(this.containsKey(kmer)){
			return this.get(kmer);
		}else{
			return null;
		}
	}
	
	public T get(StringBuffer kmer_stringBuffer){
		return this.get(kmer_stringBuffer.toString());
	}
	
	public boolean addKmer(BitSet kmer,T count){
//		if(this.containsKey(kmer)){
//			this.put(kmer, this.get(kmer)+count);
//		}else{
			this.put(kmer, count);
//		}
		return true;
	}
	
	public void print(){
		/**System.err.print("Sorting ...");
		List<Entry<BitSet, Integer>> list =new LinkedList<Entry<BitSet,Integer>>(this.entrySet());
		Collections.sort(list,new Comparator<Entry<BitSet,Integer>>() {

			@Override
			public int compare(Entry<BitSet, Integer> arg0,
					Entry<BitSet, Integer> arg1) {
				final BitSet bs0=arg0.getKey(),bs1=arg1.getKey();
				for(int i=0;i<bs0.length();++i){
					if(bs0.get(i)){
						if(!bs1.get(i)){
							return 1;
						}
					}else{
						if(bs1.get(i)){
							return -1;
						}
					}
				}
				return 0;
			}
		});
		System.err.println(" Done");*/
		Iterator<Entry<BitSet, T>> it = this.entrySet().iterator();
		while(it.hasNext()){
			Entry<BitSet, T> e= it.next();
			System.out.println(e.getValue()+"\t"+kt_utils.bitSetToString(e.getKey()));
		}
	}
}

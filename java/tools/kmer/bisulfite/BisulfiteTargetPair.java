package tools.kmer.bisulfite;

import java.util.BitSet;

public class BisulfiteTargetPair {

	protected final int kmer;
	protected final BitSet path;
	
	public BisulfiteTargetPair(int kmer, BitSet path){
		this.kmer= kmer;
		this.path= path;
	}
}

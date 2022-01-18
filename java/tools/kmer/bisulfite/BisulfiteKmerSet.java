package tools.kmer.bisulfite;

import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;

import cern.colt.list.IntArrayList;
import cern.colt.list.ObjectArrayList;

//import com.google.common.hash.BloomFilter;
//import com.google.common.hash.Funnel;
//import com.google.common.hash.PrimitiveSink;

import tools.kmer.KmerSet_binary_utils;
import utils.perfectHash.EasyPerfectMinimalHashing;

public class BisulfiteKmerSet {
	//private HashMap<BitSet,BisulfiteKmerSparseMatrix> kmers;
	private final BisulfiteCollapsedGroup[] kmers;
	private final HashSet<BitSet> keys;
//	private final BloomFilter<BitSet> keys;
	final private KmerSet_binary_utils kt_utils;
	private final char unmethylated,methylated;
	final EasyPerfectMinimalHashing<BitSet> epmh;
	ObjectArrayList equivalences;
	
	public BisulfiteKmerSet(int kmerSize, EasyPerfectMinimalHashing<BitSet> epmh, final char unmethylated, final char methylated){
		this.epmh=epmh;
		kmers= new BisulfiteCollapsedGroup[epmh.size()];
		keys=new HashSet<BitSet>(epmh.size()*5/4,1);
		this.unmethylated=unmethylated;
		this.methylated=methylated;
		this.equivalences=new ObjectArrayList();
//		keys= BloomFilter.create(new Funnel<BitSet>() {
//			private static final long serialVersionUID = 1L;
//
//			@Override
//			public void funnel(BitSet bs, PrimitiveSink into){
//				for(int i=0;i<bs.length();++i){
//					into.putBoolean(bs.get(i));
//				}
//			}
//		}, epmh.size(), 0.001);
		kt_utils=new KmerSet_binary_utils(kmerSize);
	}
	
	public void trimToSize(){
		for (BisulfiteCollapsedGroup kmer : kmers) {
			kmer.trimToSize();
		}
	}
	
	public void transform(){
		for (BisulfiteCollapsedGroup kmer : kmers){
			kmer.transform();
		}
	}
	
	public void generateEquvivalences(BisulfiteTarget[] targets){
		equivalences= new ObjectArrayList();
		HashMap<IntArrayList, Integer> keyStore= new HashMap<IntArrayList, Integer>();
		for(int i=0;i<kmers.length;++i){
			BisulfiteCollapsedGroup kmer=kmers[i];
			for(int j=0;j<kmer.size();++j){
				IntArrayList targetsCollapsed= kmer.getTarget(j);
				int currentEquivalence;
				ObjectArrayList kmerList;
				if(keyStore.containsKey(targetsCollapsed)){
					currentEquivalence=keyStore.get(targetsCollapsed);
					kmerList=(ObjectArrayList)equivalences.get(currentEquivalence);
				}else{
					currentEquivalence=equivalences.size();
					keyStore.put(targetsCollapsed, currentEquivalence);
					kmerList=new ObjectArrayList();
					equivalences.add(kmerList);
					for(int k=0;k<targetsCollapsed.size();k+=2){
						targets[targetsCollapsed.get(k)].addEquivalence(currentEquivalence);
					}
				}
				kmerList.add(new int[]{i,j});
			}
		}
	}
	
	public BisulfiteTargetPair addKmerTarget(String orig, BisulfiteTargetKey target)throws Exception{
		return this.addKmerTarget(orig, kt_utils.stringToBitSet(orig.replace(methylated, unmethylated)), target);
	}
	
	public BisulfiteTargetPair addKmerTarget(String orig, BitSet kmer, BisulfiteTargetKey target)throws Exception{
		return this.addKmerTarget(new StringBuffer(orig), kmer, target);
	}

	public BisulfiteTargetPair addKmerTarget(StringBuffer orig, BitSet kmer, BisulfiteTargetKey target)throws Exception{
		BisulfiteCollapsedGroup bksm;
		final int hash=epmh.hash(kmer);
		if(kmers[hash]==null){
			bksm=new BisulfiteCollapsedGroup();
			kmers[hash]=bksm;
		}else{
			bksm=kmers[hash];
		}
		keys.add(kmer);
//		keys.put(kmer);
		return new BisulfiteTargetPair(hash, bksm.addTarget(orig, methylated, unmethylated, target));
	}
	
	public BisulfiteCollapsedGroup getKmer(int hash){
		return kmers[hash];
	}

	public boolean addAvailableKmer(StringBuffer orig)throws Exception{
		return this.addAvailableKmer(orig, kt_utils.stringToBitSet(orig.toString().replace(methylated, unmethylated)));
	}

	public boolean addAvailableKmer(StringBuffer orig, BitSet kmer){
		if(keys.contains(kmer)){
//		if(keys.mightContain(kmer)){
			//			try{
			kmers[epmh.hash(kmer)].addKmer(orig, methylated, unmethylated, 1);
			//			}catch(Exception e){
			//				throw new Exception(e.getMessage()+"\n"+orig+"\t"+orig.toString().replace(methylated, unmethylated)+"\t"+kt_utils.bitSetToString(kmer)+"\t"+kmer.hashCode()+"\t"+kmer.toString());				
			//			}
			return true;
		}else{
			return false;
		}
	}
	
	public int size(){
		return kmers.length;
	}
}

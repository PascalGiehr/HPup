package tools.kmer.bisulfite;

import java.util.BitSet;
import java.util.HashSet;

import cern.colt.list.IntArrayList;
import cern.colt.list.ObjectArrayList;

//import com.google.common.hash.BloomFilter;
//import com.google.common.hash.Funnel;
//import com.google.common.hash.PrimitiveSink;

import tools.kmer.KmerSet_binary_utils;
import utils.perfectHash.EasyPerfectMinimalHashing;

public class BisulfiteKmerSet2 {
	//private HashMap<BitSet,BisulfiteKmerSparseMatrix> kmers;
	private final BisulfiteCollapsedGroup2[] kmers;
	private final HashSet<BitSet> keys;
//	private final BloomFilter<BitSet> keys;
	final private KmerSet_binary_utils kt_utils;
	private final char unmethylated,methylated;
	final EasyPerfectMinimalHashing<BitSet> epmh;
	public ObjectArrayList equivalences;
	int[] targetGroups;
	
	public BisulfiteKmerSet2(int kmerSize, EasyPerfectMinimalHashing<BitSet> epmh, final char unmethylated, final char methylated, final int nrOfTargets){
		this.epmh=epmh;
		kmers= new BisulfiteCollapsedGroup2[epmh.size()];
		keys=new HashSet<BitSet>(epmh.size()*5/4,1);
		this.unmethylated=unmethylated;
		this.methylated=methylated;
		this.equivalences=new ObjectArrayList();
		this.targetGroups= new int[nrOfTargets];
		for(int i=0;i<targetGroups.length;++i){
			targetGroups[i]=-1;
		}
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
		for (BisulfiteCollapsedGroup2 kmer : kmers) {
			kmer.trimToSize();
		}
	}
	
	
	public void generateEquvivalences(){
		equivalences= new ObjectArrayList();
		IntArrayList[] tmpTarget= new IntArrayList[targetGroups.length];
		for(int i=0;i<tmpTarget.length;++i){
			int nextId=i;
			while(targetGroups[nextId]>0){
				nextId=targetGroups[nextId];
			}
			if(tmpTarget[nextId]==null){
				tmpTarget[nextId]= new IntArrayList();
				equivalences.add(tmpTarget[nextId]);
			}
			tmpTarget[nextId].add(i);
		}
		equivalences.trimToSize();
		for(int i=0; i<equivalences.size();++i){
			((IntArrayList) equivalences.get(i)).trimToSize();
		}
		System.err.println("Done");
	}
	
	public BisulfiteTargetPair2 addKmerTarget(String orig, BisulfiteTargetKey target)throws Exception{
		return this.addKmerTarget(orig, kt_utils.stringToBitSet(orig.replace(methylated, unmethylated)), target);
	}
	
	public BisulfiteTargetPair2 addKmerTarget(String orig, BitSet kmer, BisulfiteTargetKey target)throws Exception{
		return this.addKmerTarget(new StringBuffer(orig), kmer, target);
	}

	public BisulfiteTargetPair2 addKmerTarget(StringBuffer orig, BitSet kmer, BisulfiteTargetKey target)throws Exception{
		BisulfiteCollapsedGroup2 bksm;
		final int hash=epmh.hash(kmer);
		if(kmers[hash]==null){
			bksm=new BisulfiteCollapsedGroup2();
			kmers[hash]=bksm;
		}else{
			bksm=kmers[hash];
		}
		keys.add(kmer);
//		keys.put(kmer);
		return new BisulfiteTargetPair2(hash, bksm.addTarget(orig, methylated, unmethylated, target));
	}
	
	public BisulfiteCollapsedGroup2 getKmer(int hash){
		return kmers[hash];
	}

	public boolean addAvailableKmer(StringBuffer orig)throws Exception{
		return this.addAvailableKmer(orig, kt_utils.stringToBitSet(orig.toString().replace(methylated, unmethylated)));
	}

	public boolean addAvailableKmer(StringBuffer orig, BitSet kmer){
		if(keys.contains(kmer)){
//		if(keys.mightContain(kmer)){
			//			try{
			IntArrayList targets=kmers[epmh.hash(kmer)].addKmer(orig, methylated, unmethylated, 1);
			if(targets.size()>0){
				//possibly a new equivalence
				int repGroup=targets.get(0);
				for(int i=0;i<targets.size();++i){
					int curId=targets.get(i),
							curGroup=targetGroups[curId];
					if(curGroup<0){
						curGroup=curId;
					}else{
						//walk down the path
						int nextId=curGroup;
						while(targetGroups[nextId]>0){
							nextId=targetGroups[nextId];
						}
						curGroup=nextId;
					}
					if(repGroup<curGroup){
						targetGroups[curGroup]=repGroup;
					}else if(curGroup<repGroup){
						targetGroups[repGroup]=curGroup;
						repGroup=curGroup;
					}
				}
				//update ids
				for(int i=0;i<targets.size();++i){
					int curId=targets.get(i);
					if(curId!=repGroup){
						targetGroups[curId]=repGroup;
					}
				}
				return true;
			}else{
				return false;
			}
		}else{
			return false;
		}
	}
	
	public int size(){
		return kmers.length;
	}
}

package tools.kmer.bisulfite;

import java.util.BitSet;

public class BisulfiteCollapsedKmerSet {

	private BitSet include,exclude;
	private boolean isEmpty;
	private int count;
	
	public BisulfiteCollapsedKmerSet(){
		isEmpty=true;
		count=0;
	}
	
	public BisulfiteCollapsedKmerSet(BitSet include){
		this();
		this.include= (BitSet) include.clone();
		isEmpty=false;
		exclude=new BitSet();
	}
	
	private void initialize(){
		include=new BitSet();
		exclude=new BitSet();
		isEmpty=false;
	}
	
	public void union(BitSet toAdd){
		if(isEmpty){
			this.initialize();
		}
		include.or(toAdd);
	}
	
	public void intersect(BitSet toSlice){
		if(isEmpty){
			this.initialize();
		}
		include.and(toSlice);
	}
	
	public void subtract(BitSet toSubract){
		if(isEmpty){
			this.initialize();
		}
		exclude.or(toSubract);
		checkEmpty();
	}
	
	private void checkEmpty(){
		int nextPos=include.nextSetBit(0);
		while(nextPos!=-1){
			if(exclude.get(nextPos)){
				nextPos=include.nextSetBit(nextPos+1);
			}else{
				break;
			}
		}
		if(nextPos==-1){
			isEmpty=true;
			include=null;
			exclude=null;
		}
	}
	
	public boolean contains(BitSet set){
		if(isEmpty){
			return false;
		}
		boolean isIncluded=true,isExcluded=true;
		int nextPos=set.nextSetBit(0);
		while(nextPos!=-1 && isIncluded){
			if(isIncluded){
				if(!include.get(nextPos)){
					isIncluded=false;
				}
			}
			if(isExcluded){
				if(!exclude.get(nextPos)){
					isExcluded=false;
				}
			}
			nextPos=set.nextSetBit(nextPos+1);
		}
		return isIncluded && !isExcluded;
	}
	
	public int getCount(){
		return count;
	}
	
	public void addCount(){
		++count;
	}
	
	public void addCount(int add){
		count+=add;
	}
	
	public boolean isEmpty(){
		return isEmpty;
	}
}

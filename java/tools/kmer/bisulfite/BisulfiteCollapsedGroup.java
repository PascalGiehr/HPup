package tools.kmer.bisulfite;

import java.util.BitSet;

import cern.colt.list.IntArrayList;
import cern.colt.list.ObjectArrayList;

public class BisulfiteCollapsedGroup {
	
	private ObjectArrayList paths= new ObjectArrayList(),
			targets= new ObjectArrayList();
	private IntArrayList counts;
	private boolean transformed;
	
	protected BisulfiteCollapsedGroup(){
		transformed=false;
	}
	
	public void trimToSize(){
		paths.trimToSize();
		targets.trimToSize();
		for(int i=0;i<targets.size();++i){
			((ObjectArrayList)targets.get(i)).trimToSize();
		}
		if(transformed){
			counts.trimToSize();
		}
	}
	
	protected void transform(){
		ObjectArrayList collapsedKmerSets= new ObjectArrayList(),
				newTargets= new ObjectArrayList();
		
		
		for(long i=1;i<=Math.pow(2, paths.size());++i){
			BitSet combination=getCombination(i);
			BisulfiteCollapsedKmerSet cks= new BisulfiteCollapsedKmerSet();
			//add the included paths
			int nextIndex=combination.nextSetBit(0);
			while (nextIndex != -1){
				cks.intersect((BitSet) paths.get(nextIndex));
				nextIndex= combination.nextSetBit(nextIndex+1);
			}
			//subtract excluded paths
			nextIndex=combination.nextClearBit(0);
			while( nextIndex != -1 && !cks.isEmpty()){
				cks.subtract((BitSet) paths.get(nextIndex));
				nextIndex= combination.nextClearBit(nextIndex+1);
			}
			if(!cks.isEmpty()){
				//this is a non-empty group ==> save
				collapsedKmerSets.add(cks);
				ObjectArrayList targetList= new ObjectArrayList();
				nextIndex= combination.nextSetBit(0);
				while(nextIndex != -1){
					ObjectArrayList toAdd= (ObjectArrayList) targets.get(nextIndex);
					targetList.addAllOfFromTo(toAdd,0,toAdd.size()-1);
					nextIndex= combination.nextSetBit(nextIndex+1);
				}
				targetList.sort();
				newTargets.add(targetList);
			}
		}
		paths= collapsedKmerSets;
		targets= newTargets;
		counts= new IntArrayList();
		for(int i=0;i<paths.size();++i){
			counts.add(0);
		}
		transformed=true;
	}
	
	protected void addKmer(StringBuffer orig, final char methylated, final char unmethylated,int n){
		if(!transformed){
			System.err.println("trying to add a kmer to an untransformed BisulfiteCollapsedGroup, will transform");
			transform();
		}
		BitSet path=getPath(orig, methylated, unmethylated);
		for(int i=0;i<paths.size();++i){
			BisulfiteCollapsedKmerSet cks=((BisulfiteCollapsedKmerSet)paths.get(i));
			if(cks.contains(path)){
				counts.setQuick(i, counts.getQuick(i)+1);
			}
		}
	}
	
	protected BitSet addTarget(StringBuffer orig,char methylated, char unmethylated, BisulfiteTargetKey target) throws Exception{
		if(transformed){
			throw new Exception("adding a target to an already transformed BisulfiteCollapsedGroup");
		}
		BitSet path= getPath(orig, methylated, unmethylated);
		int pathPos=paths.indexOf(path, true);
		if(pathPos==-1){
			paths.add(path);
			ObjectArrayList targetList= new ObjectArrayList();
			targetList.add(target);
			targets.add(targetList);
		}else{
			((ObjectArrayList)targets.get(pathPos)).add(target);
		}
		return path;
	}
	
	private BitSet getPath(StringBuffer orig, char methylated, char unmethylated){
		BitSet path= new BitSet();
		int pos=0;
		char nuc;
		for(int i=0;i<orig.length();++i){
			nuc=orig.charAt(i);
			if(nuc==unmethylated){
				path.set(pos);
				++pos;
				path.set(pos);
				++pos;
			}else if(nuc==methylated){
				++pos;
				path.set(pos);
				++pos;
			}
		}
		return path;
	}
	/**
	 * adopted from http://stackoverflow.com/questions/2473597/bitset-to-and-from-integer-long
	 * @param value
	 * @return
	 */
	private BitSet getCombination(long value) {
		BitSet bits= new BitSet();
		int index= 0;
		while( value != 0L) {
			if(value % 2 !=0) {
				bits.set(index);
			}
			++index;
			value= value >>> 1;
		}
		return bits;
	}
	
	public int size(){
		return paths.size();
	}
	
	public IntArrayList getTarget(int i){
		final ObjectArrayList targetList=(ObjectArrayList) targets.get(i);
		IntArrayList target=new IntArrayList();
		for(int j=0;j<targetList.size();++j){
			final BisulfiteTargetKey btk= (BisulfiteTargetKey) targetList.get(j);
			target.add(btk.target);
			target.add(btk.pos);
		}
		return target;
	}
}
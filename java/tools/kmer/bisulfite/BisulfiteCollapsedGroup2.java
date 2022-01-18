package tools.kmer.bisulfite;

import java.util.BitSet;

import cern.colt.list.IntArrayList;
import cern.colt.list.ObjectArrayList;

public class BisulfiteCollapsedGroup2 {

	private ObjectArrayList paths= new ObjectArrayList(),
			read_paths= new ObjectArrayList(),
			target_kmers= new ObjectArrayList(),
			targets= new ObjectArrayList();
	private IntArrayList kmer_counts= new IntArrayList(),
			positions= new IntArrayList();
	
	protected BisulfiteCollapsedGroup2(){
		
	}
	
	public void trimToSize(){
		paths.trimToSize();
		read_paths.trimToSize();
		target_kmers.trimToSize();
		kmer_counts.trimToSize();
		targets.trimToSize();
		for(int i=0;i<targets.size();++i){
			((IntArrayList) targets.get(i)).trimToSize();
			((IntArrayList) target_kmers.get(i)).trimToSize();
		}
	}
	
	public int[] getCounts(int targetNr,int position){
		int[] counts= new int[]{0,0};
		if(targetNr>targets.size()){
			return counts;
		}
		int kmerPos=positions.indexOf(position);
		if(kmerPos==-1){
			return counts;
		}
		kmerPos=2*kmerPos;
		//irrelevant which letters are used...
		BitSet methylated=getPath(new StringBuffer("C"), 'C', 'T', false),
				unmethylated=getPath(new StringBuffer("T"),'C','T',false);
		
		IntArrayList target_kmer= (IntArrayList) target_kmers.get(targetNr);
		for(int i=0;i<target_kmer.size();++i){
			int kmerNr=target_kmer.get(i);
			BitSet path=(BitSet) read_paths.get(kmerNr);
			if(methylated.get(0)==path.get(kmerPos) && methylated.get(1)==path.get(kmerPos+1)){
				counts[0]+=kmer_counts.get(kmerNr);
			}else if(unmethylated.get(0)==path.get(kmerPos) && unmethylated.get(1)==path.get(kmerPos+1)){
				counts[1]+=kmer_counts.get(kmerNr);
			}else{
				System.err.println("measuring at non valid position");
			}
		}
		
		return counts;
	}
	
	protected IntArrayList addKmer(StringBuffer orig, final char methylated, final char unmethylated,int n){
		BitSet path=getPath(orig, methylated, unmethylated,false);
		//is this path already registred?
		int pathPos= read_paths.indexOf(path, true);
		if(pathPos==-1){
			//new path
			IntArrayList pathTargets= new IntArrayList();
			//check if it belongs to a target
			int readID=read_paths.size();
			for(int i=0;i<paths.size();++i){
				BitSet target_path= (BitSet) paths.get(i);
				int nextPos= path.nextSetBit(0);
				boolean included=true;
				while(nextPos!=-1 && included){
					if(!target_path.get(nextPos)){
						included=false;
					}
					nextPos=path.nextSetBit(nextPos+1);
				}
				if(included){
					((IntArrayList) target_kmers.get(i)).add(readID);
					pathTargets.addAllOf((IntArrayList) targets.get(i));
				}
			}
			if(pathTargets.size()>0){
				read_paths.add(path);
				kmer_counts.add(n);
				pathTargets.sort();
				return pathTargets;
			}else{
				return new IntArrayList();
			}
		}else{
			kmer_counts.setQuick(pathPos, kmer_counts.getQuick(pathPos)+n);
			return new IntArrayList();
		}
	}
	
	
	
	protected int addTarget(StringBuffer orig,char methylated, char unmethylated, BisulfiteTargetKey target) throws Exception{
		if(paths.size()==0){
			for(int i=0;i<orig.length();++i){
				final char nuc= orig.charAt(i);
				if(nuc==methylated || nuc==unmethylated){
					positions.add(i);
				}
			}
		}
		BitSet path= getPath(orig, methylated, unmethylated,true);
		int pathPos=paths.indexOf(path, true);
		if(pathPos==-1){
			paths.add(path);
			IntArrayList targetList= new IntArrayList();
			targetList.add(target.getTarget());
			targets.add(targetList);
			target_kmers.add(new IntArrayList());
			return paths.size()-1;
		}else{
			((IntArrayList)targets.get(pathPos)).add(target.getTarget());
			return pathPos;
		}
	}
	
	private BitSet getPath(StringBuffer orig, char methylated, char unmethylated,boolean isTarget){
		BitSet path= new BitSet();
		int pos=0;
		char nuc;
		for(int i=0;i<orig.length();++i){
			nuc=orig.charAt(i);
			if(nuc==methylated){
				path.set(pos);
				++pos;
				if(isTarget){
					path.set(pos);
				}
				++pos;
			}else if(nuc==unmethylated){
				++pos;
				path.set(pos);
				++pos;
			}
		}
		return path;
	}
	
	public int size(){
		return paths.size();
	}
	
	public String toString(){
		return "Size: "+size();
	}
}

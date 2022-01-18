package tools.kmer.bisulfite;

import java.util.ArrayList;
import java.util.BitSet;

import cern.colt.list.IntArrayList;
import cern.colt.list.ObjectArrayList;

public class BisulfiteCliqueKmer {

	ObjectArrayList targets; //contains IntArrayList:s (would be better with hashsets?)
	ObjectArrayList paths; //contains BitSet:s
	
	
	protected BisulfiteCliqueKmer(){
		targets= new ObjectArrayList();
		paths= new ObjectArrayList();
	}
	
	protected void trimToSize(){
		targets.trimToSize();
		paths.trimToSize();
		for(int i=0;i<targets.size();++i){
			((IntArrayList) targets.get(i)).trimToSize();
		}
	}
	
	protected void addTarget(StringBuffer orig, final char methylated, final char unmethylated, int targetNr){
		BitSet path= getPath(orig, methylated, unmethylated);
		
		int pathPos=paths.indexOf(path, true);
		
		if(pathPos==-1){
			paths.add(path);
			IntArrayList kmer= new IntArrayList();
			kmer.add(targetNr);
			targets.add(kmer);
		}else{
			IntArrayList kmer=(IntArrayList) targets.get(pathPos);
			if(!kmer.contains(targetNr)){
				kmer.add(targetNr);
			}
		}
	}
	
	private ArrayList<BitSet> 
	
	private BitSet getPath(StringBuffer orig, final char methylated, final char unmethylated){
		BitSet path= new BitSet();
		int pos=0;
		char nuc;
		for(int i=0;i<orig.length();++i){
			nuc=orig.charAt(i);
			if(nuc==unmethylated){
				path.set(pos);
				++pos;
			}else if(nuc==methylated){
				++pos;
			}
		}
		return path;
	}
}

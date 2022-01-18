package tools.kmer.bisulfite;

import java.util.HashMap;

public class BisulfiteKmer {

	private HashMap<Integer, Integer> pos2count;
	private int tot, occurences;
	
	
	protected BisulfiteKmer(){
		tot= 0;
		occurences= 0;
		pos2count= new HashMap<Integer, Integer>();
	}
	
	protected BisulfiteKmer(int[] positions){
		this();
		for(int pos : positions){
			pos2count.put(pos, 0);
		}
	}
	
	protected void addNewPos(Integer pos){
		if(!pos2count.containsKey(pos)){
			pos2count.put(pos, 0);
		}
	}
	
	protected void increaseTot(){
		++tot;
	}
	
	protected void increaseOccurences(){
		++occurences;
	}
	
	protected void increasePositions(String orig,char methylated,int n){
		this.increaseTot();
		for(Integer pos : pos2count.keySet()){
			if(orig.charAt(pos)==methylated){
				pos2count.put(pos, pos2count.get(pos)+n);
			}
		}
	}
	
	protected void increasePositions(String orig,char methylated){
		this.increasePositions(orig, methylated, 1);
	}
	
	protected int getTot(){
		return tot;
	}
	
	protected int getOccurences(){
		return occurences;
	}
	
	protected Integer[] getPos(){
		return pos2count.keySet().toArray(new Integer[]{});
	}
	
	protected Integer getMeth(int pos){
		return pos2count.get(pos);
	}
}

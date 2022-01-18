package tools.kmer.bisulfite;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashSet;

import cern.colt.list.ObjectArrayList;
import cern.colt.map.OpenIntDoubleHashMap;

public class BisulfiteKmerSparseMatrix {

	final int defaultValue=0,
			occurencesPos=-1,
			referredPos=-2,
			rawcountPos=-3;
	private HashSet<Integer> positions;
	
	private ObjectArrayList matrix;
	private ObjectArrayList paths;
	
//	private HashBasedTable<BitSet, Integer, Double> matrix;
	
//	private HashMap<BitSet, HashMap<Integer, Double>> matrix;
	
	protected BisulfiteKmerSparseMatrix(){
		positions=new HashSet<Integer>();
//		matrix= new HashMap<BitSet, HashMap<Integer,Double>>();
		matrix= new ObjectArrayList();
		paths= new ObjectArrayList();
//		matrix= HashBasedTable.create();
	}
	
	public void trimToSize(){
		paths.trimToSize();
		matrix.trimToSize();
		for(int i=0;i<matrix.size();++i){
			((OpenIntDoubleHashMap) matrix.get(i)).trimToSize();
		}
	}
	
	protected BitSet addTarget(StringBuffer orig,char methylated, char unmethylated,ArrayList<Integer> positionsToMeasure){
		BitSet path= getPath(orig, methylated, unmethylated);
		
//		if(matrix.containsRow(path)){
//			Map<Integer, Double> kmer= matrix.row(path);
//			kmer.put(referredPos,kmer.get(referredPos)+1);
//			for(Integer pos : positionsToMeasure){
//				positions.add(pos);
//				if(!kmer.containsKey(pos)){
//					kmer.put(pos, 0.0);
//				}
//			}
//		}else{
//			matrix.put(path, referredPos, 1.0);
//			matrix.put(path, occurencesPos, 0.0);
//			matrix.put(path, rawcountPos, 0.0);
//			for(Integer pos : positionsToMeasure){
//				positions.add(pos);
//				matrix.put(path, pos, 0.0);
//			}
//		}
		
		OpenIntDoubleHashMap kmer;
		int pathPos=paths.indexOf(path, true);
		//		kmer=(OpenIntDoubleHashMap) matrix.get(path.hashCode());
		if(pathPos==-1){
			kmer= new OpenIntDoubleHashMap();
			matrix.add(kmer);
			paths.add(path);
			kmer.put(referredPos, 1.0);
			kmer.put(occurencesPos, 0.0);
			kmer.put(rawcountPos, 0.0);
			for(Integer pos : positionsToMeasure){
				positions.add(pos);
				kmer.put(pos, 0.0);
			}
		}else{
			kmer= (OpenIntDoubleHashMap) matrix.get(pathPos);
			kmer.put(referredPos,kmer.get(referredPos)+1);
			for(Integer pos : positionsToMeasure){
				positions.add(pos);
				if(!kmer.containsKey(pos)){
					kmer.put(pos, 0.0);
				}
			}
		}
		return path;
	}
	
//	protected void addKmerBasic(StringBuffer orig, char methylated, char unmethylated, double n){
//		BitSet path= getPath(orig, methylated, unmethylated),tmp; //which positions are methylated/unconverted
//		OpenIntDoubleHashMap kmer;
//		for(int pathPos=0; pathPos<paths.size();++pathPos){
//			tmp=(BitSet) path.clone();
//			tmp.xor((BitSet) paths.get(pathPos)); //where are the differences to this key
//			tmp.and(path); //are any of the differences an unconverted position in new sequence
//			if(tmp.isEmpty()){ //if not this is a possible path
//				//for now simply add to each position
//				kmer=(OpenIntDoubleHashMap) matrix.get(pathPos);
//				kmer.put(occurencesPos, kmer.get(occurencesPos)+n);
//				for(Integer pos : positions){
//					if(orig.charAt(pos)==methylated){
//						if(kmer.containsKey(pos)){
//							kmer.put(pos, kmer.get(pos)+n);
//						}else{
//							kmer.put(pos,n);
//						}
//					}
//				}
//			}
//		}
//	}
	protected void addKmer(StringBuffer orig, final char methylated, final char unmethylated, final double n){
		//Do check if this is possible to add, this needs to be done if a bloom filter is to be used
		
		
		this.addKmerSimpleProb(orig, methylated, unmethylated, n);
	}
	
//	protected void addKmerBetaV1(StringBuffer orig, final char methylated, final char unmethylated, final double n){
//		BitSet path= getPath(orig, methylated, unmethylated),tmp; //which positions are methylated/unconverted
//		double sum=0;
//		Double singleValue,curMeth,curOccurences;
//		OpenIntDoubleHashMap kmer;
//		ArrayList<OpenIntDoubleHashMap> kmerToUpdate= new ArrayList<OpenIntDoubleHashMap>();
//		ArrayList<Double> weights= new ArrayList<Double>();
//		for(int pathPos=0;pathPos<paths.size();++pathPos){
////		for(BitSet key : matrix.keySet()){//calculate weights and save 
//			tmp=(BitSet) path.clone();
//			tmp.xor((BitSet) paths.get(pathPos)); //where are the differences to this key
//			tmp.and(path); //are any of the differences an unconverted position in new sequence
//			if(tmp.isEmpty()){ //if not this is a possible path
//				kmer=(OpenIntDoubleHashMap) matrix.get(pathPos);
//				singleValue=1.0;
//				curOccurences=kmer.get(occurencesPos);
//				if(curOccurences>0){
//					for(Integer pos : positions){
//						curMeth=kmer.get(pos);
//						if(curMeth==null){
//							kmer.put(pos, 0.0);
//							curMeth=0.0;
//						}
////						curMeth=kmer.get(pos);
//						//assumes that the whole chunk goes to this kmer... this should overestimate the probability
//						if(orig.charAt(pos)==methylated){
//							singleValue*=Probability.betaComplemented( curMeth+1, curOccurences-curMeth+1,(curMeth+n)/(curOccurences+n));
//						}else{
//							singleValue*=Probability.beta( curMeth+1, curOccurences-curMeth+1,(curMeth)/(curOccurences+n));
//						}
//					}
//				}
//				if(singleValue==0){
////					singleValue= Math.pow(0.5, positions.size());
//					System.err.println("Generated zero probability");
//				}
//				sum+=singleValue;
//				weights.add(singleValue);
//				kmerToUpdate.add(kmer);
//			}
//		}
//		
//		//normalize weights and update counts
//		for(int i=0;i<weights.size();++i){
//			singleValue=sum==0.0?n/weights.size():weights.get(i)*n/sum;
//			kmer=kmerToUpdate.get(i);
//			kmer.put(occurencesPos, kmer.get(occurencesPos)+singleValue);
//			for(Integer pos : positions){
//				if(orig.charAt(pos)==methylated){
//					if(kmer.containsKey(pos)){
//						kmer.put(pos, kmer.get(pos)+singleValue);
//					}else{
//						kmer.put(pos, singleValue);
//					}
//				}
//			}
//		}
//	}
	
	protected void addKmerSimpleProb(StringBuffer orig, char methylated, char unmethylated, double n){
		BitSet path= getPath(orig, methylated, unmethylated); //which positions are methylated/unconverted
		double sum=0;
		Double singleValue,curMeth,curOccurences;
		OpenIntDoubleHashMap kmer;
		ArrayList<OpenIntDoubleHashMap> kmerToUpdate= new ArrayList<OpenIntDoubleHashMap>();
//		ArrayList<Map<Integer, Double>> kmerToUpdate= new ArrayList<Map<Integer,Double>>();
//		Map<Integer, Double> kmer;
		ArrayList<Double> weights= new ArrayList<Double>();
		for(int pathPos=0;pathPos<paths.size();++pathPos){
//		for(BitSet key : matrix.keySet()){//calculate weights and save 
//		for(BitSet key : matrix.rowKeySet()){
			BitSet tmp=(BitSet) path.clone();
			tmp.xor((BitSet) paths.get(pathPos)); //where are the differences to this key
			tmp.and(path); //are any of the differences an unconverted position in new sequence
			if(tmp.isEmpty()){ //if not this is a possible path
				kmer= (OpenIntDoubleHashMap) matrix.get(pathPos);
				singleValue=1.0; // this was set to zero before... ERROR
				curOccurences=kmer.get(occurencesPos);
				if(curOccurences==0){
					//if no observation, assume 50% chance
					singleValue= Math.pow(0.5, positions.size());
				}else{
					for(Integer pos : positions){
						curMeth=kmer.get(pos);
						if(curMeth==null){
							kmer.put(pos, 0.0);
							curMeth=0.0;
						}
						if(orig.charAt(pos)==methylated){
							//these probabilities could be based on all kmers in a target
							singleValue*= curMeth/curOccurences;
						}else{
							singleValue*= (curOccurences-curMeth)/curOccurences;
						}
					}
				}
				sum+=singleValue;
				weights.add(singleValue);
				kmerToUpdate.add(kmer);
			}
		}
		//normalize weights and update counts
		for(int i=0;i<weights.size();++i){
			singleValue=sum==0.0?n/weights.size():weights.get(i)*n/sum;
			kmer=kmerToUpdate.get(i);
			kmer.put(occurencesPos, kmer.get(occurencesPos)+singleValue);
			for(Integer pos : positions){
				if(orig.charAt(pos)==methylated){
					final Double d=kmer.get(pos);
					if(d==null){
						kmer.put(pos, singleValue);
					}else{
						kmer.put(pos, d+singleValue);
					}
				}
			}
		}
	}
	
	protected double getOccurences(BitSet path)throws Exception{
		return this.getCount(path, occurencesPos);
	}
	
	protected double getReferred(BitSet path)throws Exception{
		return this.getCount(path, referredPos);
	}
	
	protected double getCount(BitSet path, int pos) throws Exception{
		int pathPos= paths.indexOf(path, true);
		if(pathPos==-1){
			throw new Exception("Path not available");
		}
		OpenIntDoubleHashMap kmer =(OpenIntDoubleHashMap) matrix.get(pathPos);
		
		Double d= kmer.get(pos);
		
		if(d==null){
			return 0.0;
		}
		return d;
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
			}else if(nuc==methylated){
				++pos;
			}
		}
		return path;
	}
	
}

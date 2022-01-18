package tools.kmer.bisulfite;

import java.io.BufferedWriter;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.HashSet;


public class BisulfiteTarget {

	private final String sep="\t";
	
//	private ArrayList<BisulfiteTargetPair> kmers;//, kmersCTR;
	private BisulfiteTargetPair2[] kmers;
	private HashSet<Integer> whichEquivalences;
	private Integer[] positions;
	private String name;
	private int kmerSize,nextPos;
//	private HashMap<Integer, HashMap<Integer, Integer>> map= new HashMap<Integer, HashMap<Integer,Integer>>();
//	private final char methylated='C',
//			unmethylated='T';
//			mapCTR= new HashMap<Integer, HashMap<Integer,Integer>>();
	
	public BisulfiteTarget(String name, String sequence, int kmerSize, ArrayList<Integer> positions){
		this.name=name;
		this.kmerSize=kmerSize;
		if(sequence.length()-kmerSize+1>0){
			kmers= new BisulfiteTargetPair2[sequence.length()-kmerSize+1];
		}else{
			kmers=new BisulfiteTargetPair2[0];
		}
		this.positions=positions.toArray(new Integer[]{});
		whichEquivalences= new HashSet<Integer>();
	}
	
	public int[][] getCountsSimple(BisulfiteKmerSet2 kmerDb){
		int[][] sum=new int[][] {new int[positions.length],new int[positions.length]};
				
//				sumMeth= new int[positions.length],
//				sumUnMeth= new int[positions.length];
		for(int kmerNr=0;kmerNr<kmers.length;++kmerNr){
			final int[] positionsToMeasure= positionsToMeasure(kmerNr);
			BisulfiteCollapsedGroup2 kmerGroup= kmerDb.getKmer(kmers[kmerNr].kmer);
			for(int i=0;i<positions.length;++i){
				if(positionsToMeasure[i]!=-1){
					int[] counts=kmerGroup.getCounts(kmers[kmerNr].path,positionsToMeasure[i]);
					sum[0][i]+=counts[0];
					sum[1][i]+=counts[1];
				}
			}
		}
		return sum;
	}
	
	public void addEquivalence(Integer eq){
		whichEquivalences.add(eq);
	}
	
	public void trimToSize(BisulfiteKmerSet2 kmerDb){
		for(BisulfiteTargetPair2 btp : kmers){
			kmerDb.getKmer(btp.kmer).trimToSize();
		}
	}
	
	public void add(BisulfiteTargetPair2 kmer){
		kmers[nextPos]=kmer;
		++nextPos;
	}
	
	private int[] positionsToMeasure(int kmerNr){
		//TODO: improve
//		ArrayList<Integer> positionsToMeasure= new ArrayList<Integer>();
		int[] positionsToMeasure= new int[positions.length];
		for(int i=0;i<positions.length;++i){
			final Integer posMod=positions[i]-kmerNr;
			if(posMod>-1 && posMod<kmerSize){
				positionsToMeasure[i]=posMod;
//				positionsToMeasure.add(posMod);
			}else{
				positionsToMeasure[i]=-1;
//				positionsToMeasure.add(-1);
			}
		}
		return positionsToMeasure;
	}
	
	
//	public ArrayList<Integer> getNextpositionsToMeasure(){
//		return this.positionsToMeasure(nextPos);
//	}

	

	
	public void printToStdOut(String dir, BisulfiteKmerSet2 kmerDb)throws Exception{
		//TODO: rewrite... refine ... make better
		BufferedWriter out= new BufferedWriter(new OutputStreamWriter(System.out));
		int[][] counts=getCountsSimple(kmerDb);
		for(int i=0;i<positions.length;++i){
			out.write(name+sep+dir+sep+positions[i]+sep+counts[0][i]+sep+counts[1][i]+"\n");
		}
		
//		BisulfiteTargetPair2 btp;
//		for(Integer pos : positions){
//			int kmerNr= pos-kmerSize+1; //initialize to the first kmer that overlap the position
//			kmerNr= kmerNr>-1?kmerNr:0; //catch out of bounds stuff
//			int end= pos+1;
//			end= nextPos<end?nextPos:end;
//			for(;kmerNr<end;++kmerNr){
//				btp=kmers[kmerNr];
//				BisulfiteCollapsedGroup2 kmer=kmerDb.getKmer(btp.kmer);
////				out.write(name+sep+dir+sep+pos+sep+kmerNr+sep+kmer.getReferred(btp.path)+sep+kmer.getOccurences(btp.path)+sep+kmer.getCount(btp.path, pos-kmerNr)+"\n");
//			}
//		}
		out.flush();
	}
}


























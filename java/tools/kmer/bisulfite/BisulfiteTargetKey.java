package tools.kmer.bisulfite;

public class BisulfiteTargetKey implements Comparable{

	int target,pos;

	public BisulfiteTargetKey(int target,int pos){
		this.target=target;
		this.pos=pos;
	}
	
	protected int getTarget() {
		return target;
	}

	protected int getPos() {
		return pos;
	}

	@Override
	public int compareTo(Object o) {
		if(this == o){
			return 0;
		}
		
		final BisulfiteTargetKey btk = (BisulfiteTargetKey) o;
		
		if(this.pos == btk.pos){
			return this.pos - btk.pos;
		}else{
			return this.target - btk.target;
		}
	}

	@Override
	public boolean equals(Object o){
		if (this == o){
			return true;
		}
		if (o instanceof BisulfiteTargetKey){
			final BisulfiteTargetKey btk = (BisulfiteTargetKey) o;
			return this.target == btk.target && this.pos == btk.pos;
		}
		return false;
	}
	
	@Override
	public int hashCode(){
		int hash= 1;
		hash = hash * 17 + target;
		hash = hash * 31 + pos;
		
		return hash;
	}
	
}

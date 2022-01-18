package tools.fastq;

import java.io.BufferedWriter;
import java.io.Serializable;
import java.util.ArrayList;

import tools.fasta.FastaSeq;
import tools.rocheQual.RocheQualSeq;

public class FastqSeq implements Serializable{

	
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private String header;
	private StringBuffer seq,quality;

	private FastqSeq(String header){
		this.header = header;
	}
	
		public FastqSeq(String header, String seq, String quality) {
		this(header);
		this.seq = new StringBuffer(seq);
		this.quality= new StringBuffer(quality);
	}
	
	public FastqSeq(String header, StringBuffer seq, StringBuffer quality){
		this(header);
		this.seq = new StringBuffer(seq);
		this.quality= new StringBuffer(quality);
	}
	
	public FastqSeq(FastqSeq fqs){
		this(fqs.header,fqs.seq,fqs.quality);
	}
	
	public FastqSeq(FastaSeq fs){
		this("@"+fs.getHeader().substring(1),fs.getSeq(),fs.getGoodSangerQualSeq());
	}
	
	private String newblerHeader(String library,String direction){
		String qname=this.getQname();
		return this.getHeader()+" dir="+direction+" library="+library+" template="+qname.substring(0, qname.length()-2);
	}
	
	public void toNewblerHeader(String library,String direction){
		this.setHeader(this.newblerHeader(library, direction));
	}
	
	public String getHeader() {
		return header;
	}
	public String getQname() {
		return header.substring(1).split(" ")[0];
	}
	
	public void addToQname(String extra){
		String[] l=header.split(" ");
		StringBuffer newHeader=new StringBuffer(l[0]);
		newHeader.append(extra);
		for(int i=1;i<l.length;++i){
			newHeader.append(" ");
			newHeader.append(l[i]);
		}
		this.setHeader(newHeader.toString());
	}

	public void setHeader(String header) {
		this.header = header;
	}

	public String getSeq() {
		return seq.toString();
	}
	
	public String getSeq(int start,int stop){
		return seq.substring(start, stop);
	}
	
	public int length(){
		return seq.length();
	}

	public void setSeq(String seq) {
		this.seq = new StringBuffer(seq);
	}
	
	public void setSeq(StringBuffer seq){
		this.seq = new StringBuffer(seq);
	}
	
	public String getQuality() {
		return quality.toString();
	}
	
	public String getQuality(int start,int stop){
		return quality.substring(start, stop);
	}

	public void setQuality(String quality) {
		this.quality = new StringBuffer(quality);
	}
	
	public void setQuality(StringBuffer quality) {
		this.quality = new StringBuffer(quality);
	}
	
	public void trimThis(int start, int end){
		setQuality(this.quality.substring(start, end));
		setSeq(this.seq.substring(start, end));
		
	}

	public void trimThis(int length){
		trimThis(0,length);
	}
	
	public void correct(ArrayList<FastqCorrection> corrections)throws Exception{
		String s="",q="";
		int lastpos= 0;
		for (FastqCorrection correction : corrections) {
			if(this.getSeq().charAt(correction.getPos())!=correction.getFrom()){
				throw new Exception("Trying to apply multiple changes and cannot correct \n"+this.toString()+"\nwith correction:\n"+correction.toString());
			}
			s+=getSeq().substring(lastpos, correction.getPos())+correction.getTo();
			q+=getQuality().substring(lastpos, correction.getPos())+((char)correction.getQual());
			lastpos=correction.getPos()+1;
		}
		this.setSeq(s+getSeq().substring(lastpos));
		this.setQuality(q+getQuality().substring(lastpos));
		
	}
	
	public void correct(FastqCorrection correction)throws Exception{
		if(this.getSeq().charAt(correction.getPos())!=correction.getFrom()){
			throw new Exception("cannot correct \n"+this.toString()+"\nwith correction:\n"+correction.toString());
		}
		this.setSeq(this.getSeq(0, correction.getPos())+correction.getTo()+this.getSeq().substring(correction.getPos()+1));
		this.setQuality(this.getQuality().substring(0, correction.getPos())+((char)correction.getQual())+this.getQuality().substring(correction.getPos()+1));
	}
	
	public void write(BufferedWriter out)throws Exception{
		out.write(header);
		out.write('\n');
		out.write(seq.toString());
		out.write("\n+\n");
		out.write(quality.toString());
		out.write('\n');
	}

	public String toString(){
		return header+"\n"+seq+"\n+\n"+quality;
	}
	
	
	
	public FastaSeq toFastaSeq(){
		return new FastaSeq(">"+this.getHeader().substring(1),this.seq);
	}
	
	public FastqSeq getSubFastqSeq(int start, int stop)throws Exception{
		if(start<0||stop>seq.length()||stop>quality.length()){
			throw new Exception("Limits for "+header.substring(1)+", with a length of "+seq.length()+", is out of bounds ("+start+", "+stop+")");
		}
		return new FastqSeq(this.getHeader(), this.getSeq(start, stop), this.getQuality(start, stop));
	}
	
	public RocheQualSeq toPhredQualSeq(){
		return this.toQualSeq(64);
	}
	
	public FastqSeq changeQualBase(int origBase,int newBase,int limit){
		return new FastqSeq(this.getHeader(),this.seq,this.changeQualBase(origBase, newBase, true,limit));
	}
	
	private StringBuffer changeQualBase(int origBase, int newBase,boolean shift,int limit){
		StringBuffer s=new StringBuffer();
		int n;
		final int diff =origBase-newBase;
		final String curQuality=this.getQuality();
		for(int i=0;i<curQuality.length();i++){
			n=curQuality.charAt(i)-diff;
			if(n<newBase){
				System.err.println("Negative quality values when changing base value in: "+ this.getQname());
				s.append((char) newBase);
			}else if(n>limit){
				System.err.println("High quality values when changing base value in: "+ this.getQname());
				s.append((char) limit);
			}else{
				s.append((char)n);
			}
		}
		return s;
	}
	
	private StringBuffer toSangerQual(int base){
		StringBuffer s=new StringBuffer();
		int n;
		for (String number : this.getQuality().split(" +")) {
			if(number.length()>0){
				n=Integer.parseInt(number);
				if(n>=0){
					s.append((char) (base+n));
				}else{
					System.err.println("negative quality values in: "+this.getQname());
					s.append((char) (base));
				}
			}
		}
		return s;
	}
	
	public void convertNumberToPhredQual(){
		setQuality(toSangerQual(64));
	}
	
	private RocheQualSeq toQualSeq(int base){
		StringBuffer seq=new StringBuffer();
		if(this.length()>0){
			seq.append(quality.charAt(0)-base);
			for(int i=1;i<this.length();i++){
				seq.append(" "+(quality.charAt(i)-base));
			}
		}
		return new RocheQualSeq(">"+header.substring(1), seq.toString());		
	}
	
	private StringBuffer reverseQual(){
		StringBuffer revQual=new StringBuffer(this.quality);
		revQual.reverse();
		return revQual;
	}
	
	private StringBuffer reverseComplementSeq(){
		StringBuffer revCompSeq= new StringBuffer();
		for(int i=0;i<seq.length();i++){
			revCompSeq.insert(0,this.complement(seq.charAt(i)));
		}
		return revCompSeq;
	}
	
	public FastqSeq reverseComplement(){
		return new FastqSeq(this.header, this.reverseComplementSeq(), this.reverseQual());
	}
	
	public void reverseComplementThis(){
		this.setSeq(this.reverseComplementSeq());
		this.setQuality(this.reverseQual());
	}
	
	public void reverseThis(){
		seq=seq.reverse();
		quality=quality.reverse();
	}
	
	private char complement(char c){
		switch (c) {
		case 'A':
			return 'T';
		case 'T':
			return 'A';
		case 'G':
			return 'C';
		case 'C':
			return 'G';
		case 'X':
			return 'X';
		case 'N':
			return 'N';
		default:
			System.err.println("unknown nucleotide: "+c+", will return "+c);
			return c;
		}
	}
	
	public String toString(String newHeader)throws Exception{
		return toString(0,this.length(),newHeader);
	}
	
	public String toString(int start,int stop)throws Exception{
		return this.toString(start, stop, header.substring(1));
	}
	
	public String toString(int start,int stop,String newHeader)throws Exception{
		if(start<0||stop>seq.length()){
			throw new Exception("Limits for "+header.substring(1)+", with a length of "+seq.length()+", is out of bounds ("+start+", "+stop+")");
		}
		
		return "@"+newHeader+"\n"+seq.substring(start, stop)+"\n+\n"+quality.substring(start,stop);
	}
	
	public int indexHammingDist(char[] index) throws Exception{
		char[] readIndex= header.substring(header.lastIndexOf(':')+1).toCharArray();
		
		if(index.length!=readIndex.length){
			throw new Exception("The parsed index ("+header.substring(header.lastIndexOf(':')+1)+") have a different length than the given index ("+(index.length)+")");
		}
		
		int dist=0;
		for(int i=0;i<index.length;++i){
			dist+=(index[i]!=readIndex[i]?1:0);
		}
		
		return dist;
	}
}

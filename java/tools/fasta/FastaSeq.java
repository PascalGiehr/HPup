package tools.fasta;

import java.io.BufferedWriter;

public class FastaSeq {

	private String header,qname;
	private StringBuffer seq;

	private FastaSeq(String header){
		this.header = header;
		this.updateQname();
	}
	
	public FastaSeq(String header, String seq) {
		this(header);
		this.seq = new StringBuffer(seq);
	}
	
	public FastaSeq(String header, StringBuffer seq){
		this(header);
		this.seq= new StringBuffer(seq);
	}
	
	public String getGoodSangerQualSeq(){
		return getSangerQualSeq(30);
	}
	
	private String getSangerQualSeq(int quality){
		return getQualSeq(quality,33);
	}
	
	private String getQualSeq(int quality, int base){
		StringBuffer qual=new StringBuffer();
		for(int i=0;i<this.length();i++){
			qual.append((char) (base+quality));
		}
		return qual.toString();
	}

	public String getHeader() {
		return header;
	}

	public void setHeader(String header) {
		this.header = header;
		this.updateQname();
	}

	public String getSeq() {
		return seq.toString();
	}
	
	public String getSeq(int start,int end){
		return seq.substring(start, end);
	}
	
	public int length(){
		return seq.length();
	}

	public void setSeq(String seq) {
		this.seq = new StringBuffer(seq);
	}
	
	public String getQname(){
		return qname;
	}
	
	public String toString(){
		return header+"\n"+seq;
	}
	
	public void write(BufferedWriter out)throws Exception{
		out.write(header);
		out.write('\n');
		out.write(seq.toString());
		out.write('\n');
		
	}
	
	public String toString(int start,int stop)throws Exception{
		if(start<0||stop>seq.length()){
			throw new Exception("Limits for "+qname+", with a length of "+seq.length()+", is out of bounds ("+start+", "+stop+")");
		}
		
		return header+"\n"+this.getSeq(start, stop);
	}
	
	private void updateQname(){
		qname=header.substring(1).split(" ")[0];
	}
	
}

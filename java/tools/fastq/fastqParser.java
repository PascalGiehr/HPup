package tools.fastq;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;

public class fastqParser {

	private BufferedReader reader;
	private boolean hasmore= true;
	private String s=null;
	private String currHit= null;
	private String nextHit=null;
	private StringBuffer currSeq=new StringBuffer();
	private StringBuffer nextSeq;
	private StringBuffer currQual;
	private StringBuffer nextQual;
	private String prefix=null;
	
	
	public fastqParser(String file) throws Exception{
		this(new BufferedReader(new FileReader(file)),"");
	}
	
	public fastqParser(String file, boolean gzipped)throws Exception{
		if(gzipped){
			this.reader=new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file))));
		}else{
			this.reader=new BufferedReader(new FileReader(file));
		}
		this.prefix="";
		for(s= reader.readLine();s!=null;s=reader.readLine())
			if(s.length()>0)
				if(s.charAt(0)=='@')
					break;
		readNext();
	}
	
	/**
	 * 
	 * @param reader
	 * @throws IOException
	 * @throws BioException
	 */
	public fastqParser(BufferedReader reader,String prefix) throws IOException, Exception{
		this.reader= reader;
		this.prefix= prefix;
		for(s= reader.readLine();s!=null;s=reader.readLine())
			if(s.length()>0)
				if(s.charAt(0)=='@')
					break;
		readNext();
	}
	
	public void close()throws Exception{
		reader.close();
	}
	
	public FastqSeq next() throws IOException, Exception{
		String header=nextHit();
		return new FastqSeq(header,getSeq(),getQual());
	}
	public String nextHit() throws IOException, Exception{
		if (!hasmore) {
			throw new Exception("No more hits 1");
		}
		currHit= nextHit;
		currSeq= new StringBuffer(nextSeq);
		currQual=new StringBuffer(nextQual);
		readNext();
		return currHit;
	}
	public String getSeq()throws IOException,Exception{
		/*if (currSeq.length()<1) {
			throw new Exception("No more hits 2");
		}*/
		return currSeq.toString();
	}
	public String getQual(){
		return currQual.toString();
	}
//	public String getHit()throws IOException,Exception{
//		/*if (!hasmore) {
//			throw new Exception("No more hits 3");
//		}*/
//		return currHit;
//	}
	public boolean hasNext() {
		return hasmore;
	}
	private void readNext() throws IOException, Exception{
		//assumes that both the sequence and the corresponding quality values only occupy one line each.
//		nextSeq=new StringBuffer();
//		nextQual=new StringBuffer();
		if (s==null) {
			nextSeq=new StringBuffer();
			nextQual=new StringBuffer();
			hasmore= false;
		} else {
			nextHit=s;
			nextSeq=new StringBuffer(reader.readLine());
			nextQual=new StringBuffer();
//			s=reader.readLine();
			for(s= reader.readLine();s!=null;s=reader.readLine())
				if(s.length()>0){
					if(s.charAt(0)=='+')
						break;
					else{
						nextSeq.append(s);
					}
				}
			if(!nextHit.substring(1).equals(s.substring(1))&&!(s.equals("+"))){
				System.err.println("Strange file:\n"+nextHit+" is not equal to "+s);
			}
			int line=0;
			for(s= reader.readLine();s!=null;s=reader.readLine(),line++){
				if(s.length()>0){
					if(s.startsWith("@")&&nextQual.length()>=nextSeq.length())
						if(s.startsWith("@"+prefix)){
							break;
						}else{
//							System.err.println("quality line:\n"+s);
							nextQual.append(s);
						}
					else{
						nextQual.append(s);
					}
				}
//				if(line>1){
//					System.err.println(s);
//				}
			}
//			nextQual=reader.readLine();
//			s=reader.readLine();
		}
	}
}

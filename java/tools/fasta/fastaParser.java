package tools.fasta;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

public class fastaParser {
	private BufferedReader reader;
	private boolean hasmore= true;
	String s=null;
	private String currHit= null;
	private String nextHit=null;
	private StringBuffer currSeq= new StringBuffer();
	private StringBuffer nextSeq=null;
	
	
	public fastaParser(String fileName) throws IOException, Exception{
		this(new BufferedReader(new FileReader(fileName)));
	}
	/**
	 * 
	 * @param reader
	 * @throws IOException
	 * @throws BioException
	 */
	public fastaParser(BufferedReader reader) throws IOException, Exception{
		this.reader= reader;
		for(s= reader.readLine();s!=null;s=reader.readLine())
			if(s.length()>0)
				if(s.charAt(0)=='>')
					break;
		readNext();
	}
	
	public void close() throws IOException{
		reader.close();
	}
	
	public FastaSeq next() throws IOException, Exception{
		String header=nextHit();
		return new FastaSeq(header,getSeq());
	}
	public String nextHit() throws IOException, Exception{
		if (!hasmore) {
			throw new Exception("No more hits 1");
		}
		currHit= nextHit;
		currSeq=new StringBuffer(nextSeq);
		readNext();
		return currHit;
	}
//	public String getSeqOrig()throws IOException,Exception{
//		/*if (currSeq.length()<1) {
//			throw new Exception("No more hits 2");
//		}*/
//		return currSeq;
//	}
	public String getSeq()throws IOException,Exception{
		/*if (currSeq.length()<1) {
			throw new Exception("No more hits 2");
		}*/
		return currSeq.toString();
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
		nextSeq=new StringBuffer();
		if (s==null) {
			hasmore= false;
		} else {
			nextHit=s;
			for(s= reader.readLine();s!=null;s=reader.readLine())
				if(s.length()>0){
					if(s.charAt(0)=='>')
						break;
					else{
						nextSeq.append(s);
					}
				}
		}
	}
}

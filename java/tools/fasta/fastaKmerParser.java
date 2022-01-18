package tools.fasta;

import java.io.BufferedReader;
import java.io.IOException;

public class fastaKmerParser {
	private BufferedReader reader;
	private boolean hasmore= true;
	private String s=null;
	private int curPos; 
	private int length;
	private String currKmer="";
	private StringBuffer nextKmer;
	private String curQname="",nextQname=null;
//	private boolean restart;
	
	/**
	 * 
	 * @param reader
	 * @throws IOException
	 * @throws BioException
	 */
	public fastaKmerParser(BufferedReader reader,int length) throws IOException, Exception{
		this.reader= reader;
		this.length=length;
		curPos=0;
		for(s= reader.readLine();s!=null;s=reader.readLine())
			if(s.length()>0)
				if(s.charAt(0)=='>')
					break;
		
		currKmer="";
		readNext();
		
	}
	
	public String getQname(){
		return curQname;
	}
	
	public String next() throws IOException, Exception{
		if (!hasmore) {
			throw new Exception("No more hits 1");
		}
		curQname=nextQname;
		currKmer=nextKmer.toString();
		readNext();
		return currKmer;
	}
//	public String getSeqOrig()throws IOException,Exception{
//		/*if (currSeq.length()<1) {
//			throw new Exception("No more hits 2");
//		}*/
//		return currSeq;
//	}
	public String getKmer()throws IOException,Exception{
		/*if (currSeq.length()<1) {
			throw new Exception("No more hits 2");
		}*/
		return currKmer;
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
	
	public void close()throws Exception{
		reader.close();
	}
	
	
	private void readNext() throws IOException, Exception{
		if (s==null) {
			hasmore= false;
		} else {
			if(s.charAt(0)=='>'){
				nextQname=s.split(" ")[0].substring(1);
				//restart
				s=reader.readLine();
				while(s.length()<length){
					String tmp=reader.readLine();
					if(tmp==null||(tmp.length()>0&&tmp.charAt(0)=='>')){
						s=tmp;
						break;
					}
					s+=tmp;
				}
				if(s.charAt(0)=='>'){
					readNext();
				}else{
					nextKmer=new StringBuffer(s.substring(0, length));
					curPos=length;
				}
			}else if(s.length()==curPos){
				//read new line
				s=reader.readLine();
				curPos=0;
				readNext();
			}else{
				nextKmer.deleteCharAt(0);
				nextKmer.append(s.charAt(curPos));
//				nextKmer=nextKmer.substring(1)+s.charAt(curPos);
				curPos++;
			}
		}
	}
}

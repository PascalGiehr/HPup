package tools.fasta;

import java.nio.ByteBuffer;
import java.nio.CharBuffer;
import java.nio.channels.SeekableByteChannel;
import java.nio.charset.Charset;
import java.nio.charset.CharsetDecoder;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.HashMap;
import java.util.Set;

public class RandomAccessFastaParser {

	SeekableByteChannel sbc;
	ByteBuffer bb= ByteBuffer.allocate(2048);
	CharBuffer cb;
	HashMap<String, Long> qname2position;
	long maxSeek;
	FastaSeq curSeq,nextSeq;
	boolean hasNext;
	
	final static String newline="\n",
			charmap="ISO-8859-1",
			figmentedQname="This is not a qname";
	
	final static CharsetDecoder decoder= (Charset.forName(charmap)).newDecoder();
	
	public RandomAccessFastaParser(String file) throws Exception{
		sbc= Files.newByteChannel(Paths.get(file), StandardOpenOption.READ);
		qname2position= new HashMap<String, Long>();
		hasNext= false;
		maxSeek=0;
		curSeq= new FastaSeq(">NotASequnce",new StringBuffer());
	}
	
	public Set<String> getAllNames()throws Exception{
		if(seek(figmentedQname)){
			throw new Exception("ERROR: Don't use \""+figmentedQname+"\" as a the first part of your header lines");
		}
		return qname2position.keySet();
	}
	
	
	public FastaSeq get(final String chr)throws Exception{
		if(curSeq.getQname().equals(chr)){
			return curSeq;
		}
		readNext(chr);
		return curSeq;
	}
	
	private void readNext(final String chr)throws Exception{
		if(!qname2position.containsKey(chr)){
			if(!seek(chr)){
				throw new Exception("ERROR: the fasta file doesn't contain the sequence:\n"+chr);
			}
		}
		sbc.position(qname2position.get(chr));
		StringBuffer seq= new StringBuffer();
		String header=null,s, remains="";
		String[] l;
		boolean completeSeq=false;
		int readBytes=0;
		while(!completeSeq && (readBytes=sbc.read(bb))>0){
			bb.flip();
//			s= new String(bb.array());
			s= decoder.decode(bb).toString();
			bb.clear();
			int i=0;
			if(header==null){
				s=remains+s;
			}
			l=s.split(newline);
			if(header==null){
				for(;i<l.length-1 && header==null;++i){
					if(l[i].startsWith(">"+chr)){
						header=l[i];
					}
				}
				if(header==null && i==l.length-1){
					//check last line
					if(s.endsWith(newline) && l[i].startsWith(">"+chr)){
						header=l[i];
					}else{
						remains=l[i];
					}
					++i;
				}
			}
			for(;i<l.length && !completeSeq;++i){
				if(!(completeSeq=l[i].startsWith(">"))){
					seq.append(l[i]);
				}
			}
		}
		if(sbc.position()-bb.capacity()>maxSeek){
			maxSeek=sbc.position()-bb.capacity();
		}
		curSeq= new FastaSeq(header, seq);
	}
	
	private boolean seek(final String chr)throws Exception{
		boolean chrFound=false;
		long curPosition;
//		if(sbc.position()<maxSeek){
//			sbc.position(maxSeek);
//		}
		sbc.position(maxSeek);
		int readBytes=0;
		while(!chrFound && (readBytes=sbc.read(bb))>0){
			bb.flip();
			String s= decoder.decode(bb).toString();
			curPosition=sbc.position()-bb.capacity();
			if(curPosition>maxSeek){
				maxSeek=curPosition;
			}
			bb.clear();
			
			String[] l=s.split(newline);
			
			//extend until the last row doesn't contain a broken qname
			if(l[l.length-1].startsWith(">")){
				StringBuffer sb= new StringBuffer(s);
				boolean extend=true;
				while(extend){
					if(extend=(!s.endsWith(newline))){
						final String lastLine=l[l.length-1];
						if(extend=(lastLine.startsWith(">") && !lastLine.contains(" "))){
							if(extend=((readBytes=sbc.read(bb))>0)){
								bb.flip();
								s= decoder.decode(bb).toString();
								long seeked=sbc.position()-bb.capacity();
								if(seeked>maxSeek){
									maxSeek=seeked;
								}
								bb.clear();
								sb.append(s);
								l=s.split(newline);
								if(l.length==1){
									l[0]=lastLine+l[0];
								}
							}
						}
					}
				}
				s=sb.toString();
				l=s.split(newline);
			}
			for(int i=0;i<l.length;++i){
				if(l[i].startsWith(">")){
					//extract qname
					String qname=l[i].substring(1).split(" ")[0];
					if(qname2position.containsKey(qname)){
						if(qname2position.get(qname)!=curPosition){
							throw new Exception("ERROR: the fasta file contains non-unique names (everything preceding the first space in a header line is considered the name):\n"+qname);
						}
					}
					qname2position.put(qname, curPosition);
					if(chr.equals(qname)){
						chrFound=true;
					}
				}
			}
		}
		return chrFound;
	}
}

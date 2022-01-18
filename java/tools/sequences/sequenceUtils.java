package tools.sequences;

public class sequenceUtils {

	public static char complement(char c){
		switch (c) {
		case 't':
			return 'a';
		case 'T':
			return 'A';
		case 'a':
			return 't';
		case 'A':
			return 'T';
		case 'g':
			return 'c';
		case 'G':
			return 'C';
		case 'c':
			return 'g';
		case 'C':
			return 'G';
		case 'N':
			return 'N';
		case 'n':
			return 'n';
			default:
//				System.err.println("tools.sequences.sequenceUtils.complement: Unknown char: "+c+", will return "+c);
				return c;
//			throw new Exception("Unknown character: "+c);
		}
	}
	
	public static String Complement(String in){
		StringBuffer seqOut=new StringBuffer();
		for (char c : in.toCharArray()) {
			seqOut.append(complement(c));
		}
		return seqOut.toString();
	}
	
	public static String reverseComplement(String in){
		StringBuffer seqOut=new StringBuffer();
		for (char c : in.toCharArray()) {
			seqOut.insert(0,complement(c));
		}
		return seqOut.toString();
	}
}

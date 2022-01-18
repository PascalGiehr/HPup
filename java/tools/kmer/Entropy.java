package tools.kmer;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class Entropy {
	protected List<Character> alphabet = null;

	
	public Entropy(List<Character> alphabet)
	{
		this.alphabet=alphabet;		
	}
	
	public Entropy(){
		alphabet=new ArrayList<Character>();
		alphabet.add('A');
		alphabet.add('C');
		alphabet.add('G');
		alphabet.add('T');
	}
	
	protected double entropy(List<Double>prob)
	{
		double sum=0;
		for(int i=0; i < prob.size();++i)
		{
			if(prob.get(i)!=0)
			sum+=prob.get(i)*(Math.log(prob.get(i))/Math.log(2));
		}
		sum*=-1;
		return sum;
	}
	
	protected List<Double> calcProb(String string)
	{
		List<Double>prob = new ArrayList<Double>();
		Map<Character, Integer> distr = new HashMap<Character, Integer>();
		for(int i=0; i <string.length();++i)
		{
			if(!distr.containsKey(string.charAt(i)))
			{
				distr.put(string.charAt(i),1);
			}
			else
			{
				distr.put(string.charAt(i),distr.get(string.charAt(i))+1);
			}
		}
		for(int i=0; i < this.alphabet.size();++i)
		{
			if(distr.get(this.alphabet.get(i))==null)
			{
				prob.add(0.0);
			}
			else
			prob.add((double)distr.get(this.alphabet.get(i))/(double)string.length());
		}
		return prob;
	}
	
	public double calcMaxEntropy()
	{
		List<Double> prob = new ArrayList<Double>();
		double alphabetLength = this.alphabet.size();
		for(int i=0; i<alphabetLength;++i)
		{
			prob.add(1.0/alphabetLength);
		}
		return entropy(prob);
	}
	
	
	public double calcEntropy(String string)
	{
		return entropy(calcProb(string));
	}
}

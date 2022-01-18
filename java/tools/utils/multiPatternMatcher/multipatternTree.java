package tools.utils.multiPatternMatcher;

import java.util.ArrayList;


//this structure is completely expanded and contains 4^length nodes

public class multipatternTree {

	private multipatternTreeNode root;
	
	public multipatternTree(){
		root= new multipatternTreeNode();
		root.targets=null;
		root.a=new multipatternTreeNode();
		root.c=new multipatternTreeNode();
		root.g=new multipatternTreeNode();
		root.t=new multipatternTreeNode();
	}
	
	public void consume(String seq, String tag){
		root.consume(seq.toUpperCase(),tag);
	}
	
	public ArrayList<String> get(String seq){
		return root.get(seq);
	}
	
	public static class multipatternTreeNode{
		private ArrayList<String> targets;
		private multipatternTreeNode c,t,a,g;
		
		protected ArrayList<String> get(String seq){
			if(seq.length()==0){
				if(targets==null){
					return new ArrayList<String>();
				}else{
					return targets;
				}
			}
			switch (seq.charAt(0)) {
			case 'A':
				if(a==null)
					return new ArrayList<String>();
				return a.get(seq.substring(1));
			case 'C':
				if(c==null)
					return new ArrayList<String>();
				return c.get(seq.substring(1));
			case 'G':
				if(g==null)
					return new ArrayList<String>();
				return g.get(seq.substring(1));
			case 'T':
				if(t==null)
					return new ArrayList<String>();
				return t.get(seq.substring(1));

			default:
				return new ArrayList<String>();
			}
		}
		
		protected void consume(String seq, String tag){
			if(seq.length()==0){
				if(targets==null){
					targets= new ArrayList<String>();
				}
				targets.add(tag);
			}else{
				if(a==null){
					a= new multipatternTreeNode();
					c= new multipatternTreeNode();
					g= new multipatternTreeNode();
					t= new multipatternTreeNode();
				}
				switch (seq.charAt(0)) {
				case 'A':
					a.consume(seq.substring(1), tag);
					break;
				case 'C':
					c.consume(seq.substring(1), tag);
					break;
				case 'G':
					g.consume(seq.substring(1), tag);
					break;
				case 'T':
					t.consume(seq.substring(1), tag);
					break;
				case 'R':
					a.consume(seq.substring(1), tag);
					g.consume(seq.substring(1), tag);
					break;
				case 'Y':
					c.consume(seq.substring(1), tag);
					t.consume(seq.substring(1), tag);
					break;
				case 'S':
					g.consume(seq.substring(1), tag);
					c.consume(seq.substring(1), tag);
					break;
				case 'W':
					a.consume(seq.substring(1), tag);
					t.consume(seq.substring(1), tag);
					break;
				case 'K':
					g.consume(seq.substring(1), tag);
					t.consume(seq.substring(1), tag);
					break;
				case 'M':
					a.consume(seq.substring(1), tag);
					c.consume(seq.substring(1), tag);
					break;
				case 'B':
					c.consume(seq.substring(1), tag);
					g.consume(seq.substring(1), tag);
					t.consume(seq.substring(1), tag);
					break;
				case 'D':
					a.consume(seq.substring(1), tag);
					g.consume(seq.substring(1), tag);
					t.consume(seq.substring(1), tag);
					break;
				case 'H':
					a.consume(seq.substring(1), tag);
					c.consume(seq.substring(1), tag);
					t.consume(seq.substring(1), tag);
					break;
				case 'V':
					a.consume(seq.substring(1), tag);
					c.consume(seq.substring(1), tag);
					g.consume(seq.substring(1), tag);
					break;
				case 'N':
					a.consume(seq.substring(1), tag);
					c.consume(seq.substring(1), tag);
					g.consume(seq.substring(1), tag);
					t.consume(seq.substring(1), tag);
					break;
				default:
					break;
				}
			}
		}
	}
}

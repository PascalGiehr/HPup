package tools.utils;

import java.util.ArrayList;

public class ObjectMatrix {
	private int height;
	private int width;
	private ArrayList matrix;

	/**
	 * creates an empty matrix with the defined height and width
	 * indecies starts with zero, height and width defined in the same way as for java arrays
	 * An matrix with height 3 has 3 rows with indecies 0,1,2
	 * @param height
	 * @param width
	 */
	public ObjectMatrix(int height,int width){
		this.height=height;
		this.width=width;
		matrix=new ArrayList();
		for (int i = 0; i < height*width; i++) {
			matrix.add(i,null);
		}
	}
	/**
	 * indicies starts with zero!
	 * @param i= rownumber
	 * @param j= columnnumber
	 * @return object at place i;j
	 * @throws Exception if indicies out of range
	 */
	public Object get(int i,int j) throws Exception{
		if(i>=height)
			throw new Exception("rowindex out of range");
		if(j>=width)
			throw new Exception("columnindex out of range");
		return matrix.get(i*width+j);
	}
	/**
	 * indicies starts with zero!
	 * @param i= row number
	 * @param j= column number
	 * @throws Exception if indicies out of range
	 */
	public void set(int i,int j,Object c)throws Exception{
		if(i>=height)
			throw new Exception("rowindex out of range");
		if(j>=width)
			throw new Exception("columnindex out of range");
		matrix.set(i*width+j,c);
	}
	public int getHeight() {
		return height;
	}
	public int getWidth() {
		return width;
	}
}

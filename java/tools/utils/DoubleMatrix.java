package tools.utils;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Iterator;

public class DoubleMatrix implements Serializable{
	private int height;
	private int width;
	private double[] matrix;

	/**
	 * creates an empty matrix with the defined height and width
	 * indecies starts with zero, height and width defined in the same way as for java arrays
	 * An matrix with height 3 has 3 rows with indecies 0,1,2
	 * @param height
	 * @param width
	 */
	public DoubleMatrix(int height,int width){
		this.height=height;
		this.width=width;
		matrix=new double[height*width];
	}
	public DoubleMatrix(DoubleMatrix b){
		this(b.getHeight(),b.getWidth());
		for (int i = 0; i < matrix.length; i++) {
			matrix[i]=b.get(i);
		}
	}
	/**
	 * indicies starts with zero!
	 * @param i= rownumber
	 * @param j= columnnumber
	 * @return object at place i;j
	 * @throws Exception if indicies out of range
	 */
	public double get(int i,int j) throws Exception{
		if(i>=height)
			throw new Exception("rowindex out of range");
		if(j>=width)
			throw new Exception("columnindex out of range");
		return matrix[i*width+j];
	}
	public void set(DoubleMatrix b)throws Exception{
		if(!(this.width==b.getWidth()&&this.height==b.getHeight()))
			throw new Exception("Matrix dimensions disagree");
		for (int i = 0; i < matrix.length; i++) {
			matrix[i]=b.get(i);
		}
	}
	/**
	 * indicies starts with zero!
	 * @param i= row number
	 * @param j= column number
	 * @throws Exception if indicies out of range
	 */
	public void set(int i,int j,double c)throws Exception{
		if(i>=height)
			throw new Exception("rowindex out of range");
		if(j>=width)
			throw new Exception("columnindex out of range");
		matrix[i*width+j]=c;
	}
	public int getHeight() {
		return height;
	}
	public int getWidth() {
		return width;
	}
	public DoubleMatrix getRow(int row)throws Exception{
		if(row>=this.height)
			throw new Exception("Out of Bounds");
		DoubleMatrix ans=new DoubleMatrix(1,this.width);
		int tmp=row*this.width;
		for(int i=0;i<this.width;i++){
			ans.set(i,matrix[tmp+i]);
		}
		return ans;
	}
	public double euclidianLengthFromRowTo(int row,DoubleMatrix b)throws Exception{
		if((this.width!=b.getHeight()||this.width!=b.getWidth())&&this.width!=b.getHeight()*b.getWidth()){
			throw new Exception("Matrix dimensions disagree");
		}
		double length=0;
		int tmp=row*width;
		for(int i=0;i<this.width;i++){
			length+=(b.get(i)-matrix[tmp+i])*(b.get(i)-matrix[tmp+i]);
		}
		return Math.sqrt(length);
	}
	public double squaredEuclidianLengthFromRowTo(int row,DoubleMatrix b)throws Exception{
		if((this.width!=b.getHeight()||this.width!=b.getWidth())&&this.width!=b.getHeight()*b.getWidth()){
			throw new Exception("Matrix dimensions disagree");
		}
		double length=0;
		int tmp=row*width;
		for(int i=0;i<this.width;i++){
			length+=(b.get(i)-matrix[tmp+i])*(b.get(i)-matrix[tmp+i]);
		}
		return length;
	}
	public double euclidianLength()throws Exception{
		if(!(height==1||width==1)){
			throw new Exception("Not a vector");
		}
		double length=0;
		for (int i = 0; i < matrix.length; i++) {
			length+=matrix[i]*matrix[i];
		}
		return Math.sqrt(length);
	}
	/**
	 * subtracts b from this matrix and returns a new matrix
	 * @param b
	 * @return
	 * @throws Exception
	 */
	public DoubleMatrix subtract(DoubleMatrix b)throws Exception{
		if(this.height!=b.getHeight()||this.width!=b.getWidth())
			throw new Exception("Matrix dimensions disagree");
		DoubleMatrix ans=new DoubleMatrix(this.height,this.width);
		for (int i = 0; i < matrix.length; i++) {
			ans.set(i,matrix[i]-b.get(i));
		}
		return ans;
	}
	/**
	 * subtracts b from this matrix
	 * @param b
	 * @throws Exception
	 */
	public void subtractFromThis(DoubleMatrix b)throws Exception{
		if(this.height!=b.getHeight()||this.width!=b.getWidth())
			throw new Exception("Matrix dimensions disagree");
		for (int i = 0; i < matrix.length; i++) {
			matrix[i]-=b.get(i);
		}
	}
	/**
	 * adds b from this matrix and returns a new matrix
	 * @param b
	 * @return
	 * @throws Exception
	 */
	public DoubleMatrix add(DoubleMatrix b)throws Exception{
		if(this.height!=b.getHeight()||this.width!=b.getWidth())
			throw new Exception("Matrix dimensions disagree");
		DoubleMatrix ans=new DoubleMatrix(this.height,this.width);
		for (int i = 0; i < matrix.length; i++) {
			ans.set(i,matrix[i]+b.get(i));
		}
		return ans;
	}
	/**
	 * add b from this matrix
	 * @param b
	 * @throws Exception
	 */
	public void addToThis(DoubleMatrix b)throws Exception{
		if(this.height!=b.getHeight()||this.width!=b.getWidth())
			throw new Exception("Matrix dimensions disagree");
		for (int i = 0; i < matrix.length; i++) {
			matrix[i]+=b.get(i);
		}
	}
	public void addRowToThis(DoubleMatrix b,int row)throws Exception{
		if(!(this.width==b.getWidth()||this.height==b.getWidth()))
			throw new Exception("Matrix dimensions disagree");
		int tmp=row*b.getWidth();
		for (int i = 0; i < matrix.length; i++) {
			matrix[i]+=b.get(tmp+i);
		}
	}
	public void add(int i,int j,double c)throws Exception{
		if(i>=height)
			throw new Exception("rowindex out of range");
		if(j>=width)
			throw new Exception("columnindex out of range");
		matrix[i*width+j]+=c;
	}
	public DoubleMatrix divideThis(double b){
		for (int i = 0; i < matrix.length; i++) {
			matrix[i]/=b;
		}
		return this;
	}
	public DoubleMatrix multiplyThis(double b){
		for (int i = 0; i < matrix.length; i++) {
			matrix[i]*=b;
		}
		return this;
	}
	public DoubleMatrix dotMultiplySelf(){
		for (int i = 0; i < matrix.length; i++) {
			matrix[i]=matrix[i]*matrix[i];
		}
		return this;
	}
	public void setToZero(){
		for (int i = 0; i < matrix.length; i++) {
			matrix[i]=0;
		}
	}
	public double[] max(){
		double[] max=new double[]{matrix[0],0};
		for (int i = 1; i < matrix.length; i++) {
			if(matrix[i]>max[0]){
				max[0]=matrix[i];
				max[1]=i;
			}
		}
		return max;
	}
	public double min(){
		double min=matrix[0];
		for (int i = 1; i < matrix.length; i++) {
			if(matrix[i]<min){
				min=matrix[i];
			}
		}
		return min;
	}
	public ArrayList<Integer> sort()throws Exception{
		ArrayList<Integer> ans=new ArrayList<Integer>();
		for (int i = 0; i < matrix.length; i++) {
			ans.add(new Integer(i));
		}
		return this.sort(ans);
	}
	private ArrayList<Integer> sort(ArrayList<Integer> a){
		if(a.size()<=1){
			return a;
		}else{
			double pivot=matrix[a.get(0).intValue()];
			ArrayList<Integer> pivotList=new ArrayList<Integer>();
			ArrayList<Integer> more=new ArrayList<Integer>();
			ArrayList<Integer> less=new ArrayList<Integer>();
			for (Integer i : a) {
				if(matrix[i.intValue()]==pivot){
					pivotList.add(i);
				}else if(matrix[i.intValue()]<pivot){
					less.add(i);
				}else{
					more.add(i);
				}
			}
			pivotList.addAll(more);
			pivotList.addAll(0,less);
			return pivotList;
		}
	}
	protected void set(int i,double c){
		matrix[i]=c;
	}
	public double get(int i){
		return matrix[i];
	}
	public String toString(){
		String ans="";
		for (int i = 0; i < height; i++) {
			ans+=matrix[i*width];
			for (int j = 1; j < width; j++) {
				ans+="\t"+matrix[i*width+j];
			}
			ans+="\n";
		}
		return ans;
	}
}

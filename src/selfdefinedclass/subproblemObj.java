package selfdefinedclass;

import java.util.List;
import java.util.*;
import java.lang.Integer;

public class subproblemObj {
	private List<Double> X;
	private List<Integer> Z;
	
	public subproblemObj(List<Double> X, List<Integer> Z){
		this.X = new ArrayList<Double>(X);
		this.Z = new ArrayList<Integer>(Z);
	}
	
	public List<Double> getX(){
		return X;
	}
	
	public List<Integer> getZ(){
		return Z;
	}
}

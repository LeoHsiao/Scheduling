package selfdefinedclass;
import java.util.*;
import java.lang.Integer;

public class BackupObj {
	private List<List<Double>> X;
	private List<List<List<Integer>>> Z;
	
	public BackupObj(List<List<Double>> X, List<List<List<Integer>>> Z){
		this.X.equals(X);
		this.Z.equals(Z);
	}
	
	public List<List<Double>> getX(){
		return X;
	}
	
	public List<List<List<Integer>>> getZ(){
		return Z;
	}
}

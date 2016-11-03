import java.util.*;
import java.io.*;
import java.lang.Integer;
import selfdefinedclass.Task;
import selfdefinedclass.BackupObj;
import selfdefinedclass.subproblemObj;

public class TaskScheduling {

	public static final int MAX_ITER = 350;
	public static final int MAX_REP  = 20;
	public static final double EPS  = 0.00001;
	public static final double distLimit = 50;
	
	public static void main(String[] argv) throws FileNotFoundException {
		int datan, K, T;
		double upperbound;
		/*Scanner scanner = new Scanner(System.in);
		K = scanner.nextInt();
		T = scanner.nextInt();
		
		List<Double> handling_k = new ArrayList<Double>();
		for(int k = 0; k < K; k++){
			handling_k.add(scanner.nextDouble());
		}
		
		List<Double> workload_t = new ArrayList<Double>();
		for(int t = 0; t < T; t++){
			workload_t.add(scanner.nextDouble());
		}
		
		List<List<Double>> rewards_tk = new ArrayList<List<Double>>();
		for(int t = 0; t < T; t++){
			rewards_tk.add(new ArrayList<Double>());
		}
		for(int k = 0; k < K; k++){
			for(int t = 0; t < T; t++){
				rewards_tk.get(t).add(scanner.nextDouble());
			}
		}
		
		List<Double> penalty_k = new ArrayList<Double>();
		for(int k = 0; k < K; k++){
			penalty_k.add(scanner.nextDouble());
		}
		List<List<Double>> travel_kk = new ArrayList<List<Double>>();
		for(int k = 0; k < K; k++){
			travel_kk.add(new ArrayList<Double>());
		}
		for(int k1 = 0; k1 < K; k1++){
			for(int k2 = k1; k2 < K; k2++){
				if(k1 == k2)
					travel_kk.get(k1).add(0.0);
				else{
					travel_kk.get(k1).add(scanner.nextDouble());
					travel_kk.get(k2).add(scanner.nextDouble());
				}
			}
		}
		List<Double> N_k = new ArrayList<Double>();
		for(int k = 0; k < K; k++){
			N_k.add(scanner.nextDouble());
		}
		double gamma = 0.0;
		gamma = scanner.nextDouble();*/
		//fileio
		String fileName = "./src/test2.txt";
		File input_file = new File(fileName);
		Scanner scanner = new Scanner(input_file);
			K = scanner.nextInt();
			
			T = scanner.nextInt();
			
			List<List<Double>> rewards_tk = new ArrayList<List<Double>>();
			for(int t = 0; t < T; t++){
				List<Double> temp = new ArrayList<Double>();
				for(int k = 0; k < K; k++){
					temp.add(scanner.nextDouble());
				}
				rewards_tk.add(temp);
			}
			
			List<Double> handling_k = new ArrayList<Double>();
			for(int k = 0; k < K; k++){
				handling_k.add(scanner.nextDouble());
			}
			
			List<Double> workload_t = new ArrayList<Double>();
			for(int t = 0; t < T; t++){
				workload_t.add(scanner.nextDouble());
			}
			
			List<Double> penalty_k = new ArrayList<Double>();
			for(int k = 0; k < K; k++){
				penalty_k.add(scanner.nextDouble());
			}
			
			List<Double> N_k = new ArrayList<Double>();
			for(int k = 0; k < K; k++){
				N_k.add(scanner.nextDouble());
			}
			
			List<List<Double>> travel_kk = new ArrayList<List<Double>>();
			for(int k1 = 0; k1 < K; k1++){
				List<Double> temp = new ArrayList<Double>();
				for(int k2 = 0; k2 < K; k2++){
					temp.add(scanner.nextDouble());
				}
				travel_kk.add(temp);
			}
			
			double gamma = 0.0;
			gamma = scanner.nextDouble();			
			
			List<Double> lambda_xk = new ArrayList<Double>();
			List<Double> lambda_yk = new ArrayList<Double>();
			
			for(int k = 0; k < K; k++){
				lambda_xk.add(0.0);
				lambda_yk.add(0.0);
			}
			
			List<Double> actualWorkload_t = new ArrayList<Double>();
			for(int t = 0; t < T; t++){
				actualWorkload_t.add(workload_t.get(t) * gamma);
			}
			
			double U = 0.0, U_old=0.0;
			double delta = 100.0;
			int count=0,step=0;
			
			List<List<Double>> sol_X;
			List<List<Integer>> sol_Y;
			List<List<List<Integer>>> sol_Z;
			
			List<Double> lambda1_old = new ArrayList<Double>();
			List<Double> lambda2_old = new ArrayList<Double>();
			//start subgradient
			while(true){
				if(count >20){
					delta = delta/2.0;
					count =0;
				}
				System.out.println("step:" +step +", delta:" +delta + "\nlambda1:" +lambda_xk+ "\nlambda2:" +lambda_yk);
				List<List<Double>> X_tk  = new ArrayList<List<Double>>();
				List<List<Integer>> Y_tk = new ArrayList<List<Integer>>();
				List<List<List<Integer>>> Z_tkk = new ArrayList<List<List<Integer>>>();
				//solve t subproblem
				for(int t = 0; t < T; t++){
					//using a greedy approach with highest CP value
					List<Task> tasks = new ArrayList<Task>();
					for(int k =0; k < K; k++){
						tasks.add(new Task(k, rewards_tk.get(t).get(k) - lambda_xk.get(k) - penalty_k.get(k) + lambda_yk.get(k)));
					}
					subproblemObj returnObj = subproblem(K, tasks, handling_k, travel_kk, gamma * workload_t.get(t));
					X_tk.add(returnObj.getX());
					
					List<Integer> Y_k = new ArrayList<Integer>();
					for(int k = 0; k < K; k++){
						Y_k.add((int)( Math.ceil(X_tk.get(t).get(k))));
					}
					Y_tk.add(Y_k);
					List<Integer> Xlist = returnObj.getZ();
					List<List<Integer>> Z_kk = new ArrayList<List<Integer>>();
					for(int k = 0; k < K; k++){
						Z_kk.add(new ArrayList<Integer>(Collections.nCopies(K, 0)));
					}
					//nothing chose
					if(Xlist.size() <= 1){
						Z_tkk.add(Z_kk);
						continue;
					}
					for(int i = 0; i < Xlist.size(); i++){
						if(i+1 < Xlist.size()){
							Z_kk.get(Xlist.get(i)).set(Xlist.get(i+1), 1);
						}
						else{
							Z_kk.get(Xlist.get(i)).set(Xlist.get(0), 1);
						}
					}
					Z_tkk.add(Z_kk);
				}
				System.out.println(X_tk);
				U = 0.0;
				// renew L
				for(int t = 0; t < T; t++){
					for(int k = 0; k < K; k++){
						U += rewards_tk.get(t).get(k)*X_tk.get(t).get(k);
						U += penalty_k.get(k)*Y_tk.get(t).get(k);
						
					}
				}
				for(int k = 0; k < K; k++){
					double sum_X=0.0;
					int    sum_Y=0;
					for(int t = 0; t < T; t++){
						sum_X += X_tk.get(t).get(k);
						sum_Y += Y_tk.get(t).get(k);
						
					}
					U += lambda_xk.get(k)*(1-sum_X);
					U += lambda_yk.get(k)*(N_k.get(k)-sum_Y);
				}
				if( U > U_old){
					delta = delta/2.0;
					//continue;
				}
				else{
					if(step > MAX_ITER){
						//solution
						sol_X = new ArrayList<List<Double>>(X_tk);
						sol_Y = new ArrayList<List<Integer>>(Y_tk);
						sol_Z = new ArrayList<List<List<Integer>>>(Z_tkk);
						upperbound = U;			
						break;
					}
					//calculate subgradient & lambda_new
					for(int k=0;k<K;k++){
						int    sum_Y_tk = 0;
						double sum_X_tk = 0;
						for(int t=0;t<T;t++){
							sum_Y_tk += Y_tk.get(t).get(k);
							sum_X_tk += X_tk.get(t).get(k);
						}
						lambda_xk.set(k,lambda_xk.get(k) + delta * (1 - sum_X_tk));
						lambda_yk.set(k,lambda_yk.get(k) + delta * (N_k.get(k) - sum_Y_tk));
					}
					if(step!=0){
						if((Math.abs(U_old - U)/U < EPS) && U > U_old ){
							sol_X = new ArrayList<List<Double>>(X_tk);
							sol_Y = new ArrayList<List<Integer>>(Y_tk);
							sol_Z = new ArrayList<List<List<Integer>>>(Z_tkk);
							upperbound = U;			
							break;
						}
						else{
							delta/= 2.0;
							count ++;
							if(count >20){
								System.out.println("count > 20");
								sol_X = new ArrayList<List<Double>>(X_tk);
								sol_Y = new ArrayList<List<Integer>>(Y_tk);
								sol_Z = new ArrayList<List<List<Integer>>>(Z_tkk);
								upperbound = U;
								break;
							}
						}
					}
				}
				U_old = U;
				step ++;
			}
			List<Double> capacity = new ArrayList<Double>();
			for(int t = 0; t < T; t++){
				capacity.add(gamma * workload_t.get(t));
			}
			
			
			BackupObj outObj = backup(K, T, rewards_tk, penalty_k, handling_k, travel_kk, capacity, sol_X, sol_Z);
			List<List<Double>> tempX = outObj.getX();
			List<List<List<Integer>>> tempZ = outObj.getZ();
			for(int t = 0; t < T; t++){
				System.out.println(tempX.get(t));
				for(int k = 0; k < K; k++){
					System.out.println(tempZ.get(t).get(k));
				}
				System.out.println();
			}
			
		
		scanner.close();
	}
	
	private static subproblemObj subproblem(int K, List<Task> tasks, List<Double> handling, List<List<Double>> traveling, double capacity){
		double max_val = -1.0;
		List<Double> max_X_k = new ArrayList<Double>();
		List<Integer> max_choose = new ArrayList<Integer>();
		for(int k = 0; k < K; k++){
			if(tasks.get(k).getValue() <= 0.0)
				continue;
			List<Double> X_k = new ArrayList<Double>(Collections.nCopies(K, 0.0));
			List<Integer> cur_choose = new ArrayList<Integer>();
			X_k.set(k,1.0);
			cur_choose.add(k);
			double leftCapacity = capacity - handling.get(k);
			if(capacity < 0) 
				continue;
			double value = tasks.get(k).getValue();
			int idx, n, cur_index;
			double max_cp, max_time=0;
			double cp;
			while(leftCapacity > 0){
				max_cp = 0;
				idx = -1;
				cur_index = -1;
				n = cur_choose.size();
				//find the point with highest CP value
				for(int i = 0; i < n; i++){
					int k1 = cur_choose.get(i);
					for(int kk = 0; kk < K; kk++){
						double traveltime = 0.0;
						if(i+1 != n){
							traveltime = traveling.get(k1).get(kk) + traveling.get(kk).get(cur_choose.get(i+1)) - traveling.get(k1).get(cur_choose.get(i+1));
						}
						else if(i+1 == n && n != 1){
							traveltime = traveling.get(k1).get(kk) + traveling.get(kk).get(cur_choose.get(0)) - traveling.get(k1).get(cur_choose.get(0));
						}
						else{
							traveltime = 2.0 * traveling.get(k1).get(kk);
						}
						if(cur_choose.contains(kk)) 
							continue;
						if( traveltime > leftCapacity) 
							continue;
						//cur_cp = (Reward - Penalty)/(Handling +travel)
						if(leftCapacity > traveltime){
							//can do whole task
							cp = tasks.get(kk).getValue()/traveltime;
							if(cp > max_cp){
								max_cp = cp;
								idx = kk;
								cur_index = i+1;
								max_time = traveltime;
							}
						}
						/*else{
							prop = (capacity - alpha * dist_k.get(kk)) / handling.get(k);
							cp = (prop * tasks.get(kk).getValue()) / capacity;
							if(cp > max_cp){
								max_cp = cp;
								idx = kk;
								max_time = capacity;
							}
						}*/
					}
				}
				if(idx == -1) 
					break;
				else if(max_cp <0) 
					break;
				else{
					//add that point to choosen points @k
					value += (max_cp * max_time);
					leftCapacity -= max_time;
					X_k.set(idx, 1.0);
					cur_choose.add(cur_index, idx);
				}
			}
			if(value > max_val){
				//change set of point starting with k @t
				max_val = value;
				max_X_k.clear();
				max_X_k.addAll(X_k);
				max_choose.clear();
				max_choose.addAll(cur_choose);
			}
		}
		if(max_val <= 0){
			max_X_k = new ArrayList<Double>(Collections.nCopies(K, 0.0));
			max_choose.add(-1);
		}
		System.out.println(max_X_k);
		System.out.println(max_choose);
		return (new subproblemObj(max_X_k, max_choose));
	}
	
	private static BackupObj backup(int K, int T, List<List<Double>> rewards, List<Double> penalty, List<Double> handling, List<List<Double>> traveling, List<Double> capacity, List<List<Double>> X, List<List<List<Integer>>> Z){
		List<Integer> S = new ArrayList<Integer>();
		for(int k = 0; k < K; k++){
			S.add(k);
		}
		for(int k = 0; k < K; k++){
			double sumk = 0.0;
			for(int t = 0; t < T; t++){
				sumk = sumk + X.get(t).get(k);
			}
			if(sumk == 1.0){
				S.remove(S.indexOf(k));
			}
			if(sumk > 1.0){
				List<Task> temp = new ArrayList<Task>();
				for(int t = 0; t < T; t++){
					if(X.get(t).get(k) == 1.0)
					temp.add(new Task(t, rewards.get(t).get(k)));
				}
				Task maxDay = Collections.max(temp);
				for(int t = 0; t < T; t++){
					if(t != maxDay.getId()){
						X.get(t).set(k, 0.0);
						int kto = -1, kfrom = -1;
						for(int k1 = 0; k1 < K; k1++){
							if(Z.get(t).get(k).get(k1) == 1){
								kto = k1;
								Z.get(t).get(k).set(k1, 0);
							}
							if(Z.get(t).get(k1).get(k) == 1){
								kfrom = k1;
								Z.get(t).get(k1).set(k, 0);
							}
						}
						Z.get(t).get(kfrom).set(kto, 1);
					}
				}
				S.remove(S.indexOf(k));
			}
		}
		List<Double> Used = new ArrayList<Double>();
		for(int t = 0; t < T; t++){
			Used.add(0.0);
		}
		for(int k = 0; k < K; k++){
			for(int t = 0; t < T; t++){
				if(X.get(t).get(k) == 1){
					double travel = traveling.get(k).get( Z.get(k).indexOf(Collections.max(Z.get(t).get(k))));
					Used.set(t, Used.get(t) + travel + rewards.get(t).get(k));
				}
			}
		}
		for(int t = 0; t < T; t++){
			if(Used.get(t) > capacity.get(t)){
				System.out.println("capacity error(Backup)");
				break;
			}
		}
		
		while(S.size() > 0){
			List<Double> costList = new ArrayList<Double>();
			List<Integer> tList = new ArrayList<Integer>();
			List<Integer> fromList = new ArrayList<Integer>();
			List<Integer> toList = new ArrayList<Integer>();
			for(int i = 0; i < S.size(); i++){
				List<Task> reward_k = new ArrayList<Task>();
				for(int t = 0; t < T; t++){
					reward_k.add(new Task(t, rewards.get(t).get(S.get(i))));
				}
				int t = Collections.max(reward_k).getId();
				int startk = -1;
				for(int k = 0; k < K; k++){
					if(X.get(t).get(k) == 1){
						startk = k;
						break;
					}
				}
				double mintime = Double.MAX_VALUE;
				int minfrom = -1, minto = -1;
				int k1 = startk;
				if(k1 != -1){
					for(int k = 0; k < K; k++){
						System.out.println("Z size "+Z.size());
						System.out.println("Z_t size "+Z.get(t).size());
						System.out.println("k1  "+k1);
						System.out.println("Z_t_k1 size "+Z.get(t).get(k1).size());
						if(Z.get(t).get(k1).get(k) == 1){
							if(Used.get(t) + traveling.get(k1).get(i) + traveling.get(i).get(k) - traveling.get(k1).get(k) < capacity.get(t)){
								if(traveling.get(k1).get(i) + traveling.get(i).get(k) - traveling.get(k1).get(k) < mintime){
									mintime = traveling.get(k1).get(i) + traveling.get(i).get(k) - traveling.get(k1).get(k);
									minfrom = startk;
									minto = k;
								}
							}
							k1 = k;
							if(k1 == startk)
								break;
						}
					}
				}
				if(minfrom > -1){
					double maxReward = Collections.max(reward_k).getValue();
					reward_k.remove(reward_k.indexOf(Collections.max(reward_k)));
					double cost = maxReward - Collections.max(reward_k).getValue();
					costList.add(cost);
					fromList.add(minfrom);
					toList.add(minto);
					tList.add(t);
				}
				else{
					costList.add(Double.MIN_VALUE);
					fromList.add(minfrom);
					toList.add(minto);
					tList.add(t);
				}
			}
			if(Collections.max(costList) > 0){
				int i = costList.indexOf(Collections.max(costList));
				int chosenk = S.get(i), fromk = fromList.get(i), tok = toList.get(i), t = tList.get(i);
				double xk = 0.0;
				for(int t1 = 0; t1 < T; t1++){
					xk = xk + X.get(t1).get(chosenk);
				}
				if(Used.get(t) + (handling.get(chosenk)*(1.0-xk)) + traveling.get(fromk).get(chosenk) + traveling.get(chosenk).get(tok) - traveling.get(fromk).get(tok) < capacity.get(t)){
					X.get(t).set(chosenk, (1.0-xk));
					Z.get(t).get(fromk).set(tok, 0);
					Z.get(t).get(fromk).set(chosenk, 1);
					Z.get(t).get(chosenk).set(tok, 1);
					S.remove(i);
				}
				else{
					double temp = (capacity.get(t) - (Used.get(t) + traveling.get(fromk).get(chosenk) + traveling.get(chosenk).get(tok) - traveling.get(fromk).get(tok))) / handling.get(chosenk);
					X.get(t).set(chosenk, temp);
					Z.get(t).get(fromk).set(tok, 0);
					Z.get(t).get(fromk).set(chosenk, 1);
					Z.get(t).get(chosenk).set(tok, 1);
				}
			}
		}
		BackupObj returnObj = new BackupObj(X, Z);
		return returnObj;
	}
	
}

			
	
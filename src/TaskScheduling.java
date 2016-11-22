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
		int K, T;
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
		String fileName = "./src/test3.txt";
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
			
			List<Double> lambda_xk = new ArrayList<Double>(Collections.nCopies(K, (Double)0.0));
			List<Double> lambda_yk = new ArrayList<Double>(Collections.nCopies(K, (Double)0.0));
			
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
				//System.out.println("step:" +step +", delta:" +delta + "\nlambda1:" +lambda_xk+ "\nlambda2:" +lambda_yk);
				List<List<Double>> X_tk  = new ArrayList<List<Double>>();
				List<List<Integer>> Y_tk = new ArrayList<List<Integer>>();
				List<List<List<Integer>>> Z_tkk = new ArrayList<List<List<Integer>>>();
				//solve t subproblem
				for(int t = 0; t < T; t++){
					//using a greedy approach with highest CP value
					subproblemObj returnObj = subproblem(K, rewards_tk.get(t), handling_k, travel_kk, gamma * workload_t.get(t),lambda_xk,lambda_yk,penalty_k);
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
								//System.out.println("count > 20");
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
				
					System.err.println(tempX.get(t));
				
			}
			double MAXObj = 0.0;
			for(int t = 0; t < T; t++){
				for(int k = 0; k < K; k++){
					if(tempX.get(t).get(k) > 0.0){
						MAXObj = MAXObj + tempX.get(t).get(k) * rewards_tk.get(t).get(k) - penalty_k.get(k);
					}
				}
			}
			System.err.println("reward: " + MAXObj);
			
			PrintWriter output = new PrintWriter(System.out);
			for(int k = 0; k < K; k++){
				for(int t = 0; t < T; t++){
					output.println(tempX.get(t).get(k));
				}
			}
			for(int t = 0; t < T; t++){
				int thisk = -1, nextk = -1;
				for(int k = 0; k < K; k++){
					for(int k1 = 0; k1 < K; k1++){
						if(tempZ.get(t).get(k).get(k1) == 1){
							thisk = k;
							nextk = k1;
							break;
						}
					}
					if(thisk != -1)
						break;
				}
				if(thisk != -1){
					output.println(thisk);
					while(nextk != thisk){
						output.println(nextk);
						nextk = tempZ.get(t).get(nextk).indexOf(Collections.max(tempZ.get(t).get(nextk)));
					}
				}
			}
		output.close();
		scanner.close();
	}
	
	private static subproblemObj subproblem(int K, List<Double> Reward_k, List<Double> HandlingT_k, List<List<Double>> traveling, double workload,List<Double> lambda_xk, List<Double> lambda_yk, List<Double> Penalty_k){
		double max_val = 0.0;
		List<Integer> max_choose = new ArrayList<Integer>();
		List<Double> max_X_k = new ArrayList<Double>();
		for(int k=0;k<K;k++){
			List<Double> X_k = new ArrayList<Double>(Collections.nCopies(K, 0.0));
			List<Integer> cur_choose = new ArrayList<Integer>();
			X_k.set(k,1.0);
			cur_choose.add(k);
			double capacity = workload - HandlingT_k.get(k);
			double value = Reward_k.get(k)-lambda_xk.get(k)-Penalty_k.get(k)+lambda_yk.get(k);
			if(capacity <0) continue;
			if(value <0) continue;
			int idx,n;
			double prop=0;
			double max_cp,max_time=0;
			double time,cp;
			while(capacity>0){
				max_cp = 0;
				idx = -1;
				n = cur_choose.size();
				//find the point with highest CP value
				for(int i=0;i<n;i++){
					List<Double> dist_k = traveling.get(cur_choose.get(i));
					Double alpha = 2.0; //MST to spanning tree ratio
					for(int kk=0;kk<K;kk++){
						if(cur_choose.contains(kk)) continue;
						if(alpha*dist_k.get(kk)>capacity) continue;
						//cur_cp = (Reward - Penalty)/(Handling +travel)
						if(capacity > (alpha*dist_k.get(kk)+HandlingT_k.get(k))){
							time = alpha*dist_k.get(kk)+HandlingT_k.get(k);
							//can do whole task
							cp = (Reward_k.get(kk)-lambda_xk.get(kk)-Penalty_k.get(kk)+lambda_yk.get(kk))/time;
							if(cp > max_cp){
								max_cp = cp;
								idx = kk;
								prop = 1;
								max_time = time;
							}
						}
						else{
							prop = (capacity - alpha*dist_k.get(kk)) / HandlingT_k.get(k);
							cp = (prop*(Reward_k.get(kk)-lambda_xk.get(kk))-(Penalty_k.get(kk)-lambda_yk.get(kk)))/capacity;
							if(cp > max_cp){
								max_cp = cp;
								idx = kk;
								max_time = capacity;
							}
						}
					}
				}
				if(idx == -1) break;
				else if(max_cp <0) break;
				else{
					//add that point to choosen points @k
					value += max_cp;
					capacity -= max_time;
					X_k.set(idx,prop);
					cur_choose.add(idx);
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
		return (new subproblemObj(max_X_k, max_choose));
	}
	
	private static BackupObj backup(int K, int T, List<List<Double>> rewards, List<Double> penalty, List<Double> handling, List<List<Double>> traveling, List<Double> capacity, List<List<Double>> X, List<List<List<Integer>>> Z){
		List<Integer> S = new ArrayList<Integer>();
		//把所有task放進一個集合Ｓ
		for(int k = 0; k < K; k++){
			S.add(k);
		}
		//計算所有task目前被做的次數
		for(int k = 0; k < K; k++){
			double sumk = 0.0;
			for(int t = 0; t < T; t++){
				sumk = sumk + X.get(t).get(k);
			}
			//System.out.println(sumk);
			//做的次數剛好等於1的話就把那個task從集合Ｓ中移除
			if(sumk == 1.0){
				S.remove(S.indexOf(k));
			}
			//做的次數大於1的就調整成1，並把該task從Ｓ中移除
			if(sumk > 1.0){
				//System.out.println(k);
				List<Task> temp = new ArrayList<Task>();
				//把所有會把那個task完整做完的day抓出來,並找到當中reward最大的
				for(int t = 0; t < T; t++){
					if(X.get(t).get(k) == 1.0)
					temp.add(new Task(t, rewards.get(t).get(k)));
				}
				Task maxDay = Collections.max(temp);
				//調整被移除的日期的路徑
				//後來發現這裡有一個bug，就是如果那個task做的總次數大1，但是沒有一天完整做完（＝1）會有error
				for(int t = 0; t < T; t++){
					if(X.get(t).get(k) > 0.0 && t != maxDay.getId()){
						//System.out.println(t);
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
		//計算剩的capacity
		List<Double> Used = new ArrayList<Double>();
		for(int t = 0; t < T; t++){
			Used.add(0.0);
		}
		for(int k = 0; k < K; k++){
			for(int t = 0; t < T; t++){
				if(X.get(t).get(k) == 1){
					double travel = traveling.get(k).get( Z.get(t).get(k).indexOf(Collections.max(Z.get(t).get(k))));
					Used.set(t, Used.get(t) + travel + rewards.get(t).get(k));
				}
			}
		}
		
		for(int t = 0; t < T; t++){
			if(Used.get(t) > capacity.get(t)){
				//System.out.println("capacity error(Backup)");
				break;
			}
		}
		//把不足1的調整到做完，並移除
		while(S.size() > 0){
			List<Double> costList = new ArrayList<Double>();
			List<Integer> tList = new ArrayList<Integer>();
			List<Integer> fromList = new ArrayList<Integer>();
			List<Integer> toList = new ArrayList<Integer>();
			//找到Ｓ中的各個task加進各天的reward除capacity最大的值，存進costlist（有考慮不能全部放進去的話就只算能放進去的部分的reward)
			//上次meeting好像有發現如果某天只有一個task會有問題，但實際狀況應該不會只有一個task，因為會有起始點之類的，不是很確定
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
				//如果沒有可以放進去的day就存-1，這裡min_value好像有問題
				else{
					costList.add(Double.MIN_VALUE);
					fromList.add(minfrom);
					toList.add(minto);
					tList.add(t);
				}
			}
			//找到其中reward/cost最大的放進去，如果該task做完就從Ｓ中移除
			if(Collections.max(costList) > 0){
				int i = costList.indexOf(Collections.max(costList));
				int chosenk = S.get(i), fromk = fromList.get(i), tok = toList.get(i), t = tList.get(i);
				double xk = 0.0;
				for(int t1 = 0; t1 < T; t1++){
					xk = xk + X.get(t1).get(chosenk);
				}
				if(Used.get(t) + (handling.get(chosenk)*(1.0-xk)) + traveling.get(fromk).get(chosenk) + traveling.get(chosenk).get(tok) - traveling.get(fromk).get(tok) < capacity.get(t)){
					if(X.get(t).get(chosenk) > 0.0){
						X.get(t).set(chosenk, X.get(t).get(chosenk) + (1.0-xk));
					}
					else{
						X.get(t).set(chosenk, (1.0-xk));
						Z.get(t).get(fromk).set(tok, 0);
						Z.get(t).get(fromk).set(chosenk, 1);
						Z.get(t).get(chosenk).set(tok, 1);
					}
					
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

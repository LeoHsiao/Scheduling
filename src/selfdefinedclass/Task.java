package selfdefinedclass;

public class Task implements Comparable<Task>{
	private int taskId;
	private Double value;
	
	public Task(int Id, double value){
		this.value = value;
		this.taskId = Id;
	}
	
	public int getId(){
		return taskId;
	}
	
	public double getValue() {
        return value;
    }
	
	@Override
    public int compareTo(Task emp) { 
        return this.value.compareTo(emp.getValue());
    }
}

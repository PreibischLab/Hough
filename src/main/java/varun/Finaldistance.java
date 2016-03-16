package varun;

public abstract class Finaldistance {
	
	double[] realpoints;
	double[] funcparamone;
	double funcparamtwo;
	double funcparamthree;
	
	public Finaldistance(double[] realpoints, double[] funcparamone, double funcparamtwo, double funcparamthree){
		
		this.realpoints = realpoints;
		this.funcparamone= funcparamone;
		this.funcparamtwo = funcparamtwo;
		this.funcparamthree = funcparamthree;
		
	}
	
public double Circlefunctiondist(double[] realpoints, double[] funcparamone, double funcparamtwo){
	
	double distance;

	distance = Math.abs(Math.sqrt(Math.pow((realpoints[1] - funcparamone[1]),2) 
			+ Math.pow((realpoints[0] - funcparamone[0]),2) ) - funcparamtwo);

	return distance;
	
	
}
	
}

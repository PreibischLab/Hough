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
	
abstract public double Circlefunction(double[] realpoints, double[] funcparamone, double funcparamtwo);

	
}

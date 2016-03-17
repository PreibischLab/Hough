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
	
public Finaldistance ( double [] realpoints, double funcparamtwo, double funcparamthree){
	
	this.realpoints = realpoints;
	this.funcparamtwo = funcparamtwo;
	this.funcparamthree = funcparamthree;
	
}
	
}

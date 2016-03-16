package varun;

public class Finalfunction extends Finaldistance {

	public Finalfunction(double[] realpoints, double[] funcparamone, double funcparamtwo, double funcparamthree) {
		
		super(realpoints, funcparamone, funcparamtwo, funcparamthree);
		
	}
	
	
	public double Circlefunction(double[] realpoints, double[] funcparamone, double funcparamtwo){
		
		double distance;

		distance = Math.abs(Math.sqrt(Math.pow((realpoints[1] - funcparamone[1]),2) 
				+ Math.pow((realpoints[0] - funcparamone[0]),2) ) - funcparamtwo);

		return distance;
		
		
	}
	

	
	
	

}

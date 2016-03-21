package varun;

public class Finalfunction extends Finaldistance {

	public Finalfunction(double[] realpoints, double[] funcparamone, double funcparamtwo, double funcparamthree) {

		super(realpoints, funcparamone, funcparamtwo, funcparamthree);

	}

	public Finalfunction(double[] realpoints, double funcparamtwo, double funcparamthree) {

		super(realpoints, funcparamtwo, funcparamthree);

	}

	public Finalfunction(double[] realpoints, double funcparamtwo, double funcparamthree, double funcparamfour) {

		super(realpoints, funcparamtwo, funcparamthree, funcparamfour);

	}
	
	public Finalfunction(double[] realpoints, double funcparamtwo, double funcparamthree, double funcparamfour,
			double funcparamfive) {

		super(realpoints, funcparamtwo, funcparamthree, funcparamfour, funcparamfive);

	}
	
	public Finalfunction(double[] realpoints, double funcparamtwo, double funcparamthree, double funcparamfour,
			double funcparamfive, double funcparamsix) {

		super(realpoints, funcparamtwo, funcparamthree, funcparamfour, funcparamfive, funcparamsix);

	}
	

	public double Circlefunctiondist() {

		// funcparamone = center fo the circle, funcparamtwo = radius;

		double distance;

		distance = Math.abs(Math
				.sqrt(Math.pow((realpoints[1] - funcparamone[1]), 2) + Math.pow((realpoints[0] - funcparamone[0]), 2))
				- funcparamtwo);

		return distance;

	}

	public double Linefunctiondist() {

		// funcparamtwo = slope, funcparamthree = constant

		double distance;

		double minX = (realpoints[0] - funcparamtwo * (funcparamthree - realpoints[1]))
				/ (1 + funcparamtwo * funcparamtwo);
		double minY = minX * funcparamtwo + funcparamthree;

		distance = Math.pow((minX - realpoints[0]), 2) + Math.pow((minY - realpoints[1]), 2);

		return Math.sqrt(distance);
	}

	public double Sinfunction() {
		double function;

		function = funcparamtwo * Math.sin(Math.toRadians(realpoints[0] * funcparamthree + funcparamfour));

		return function;
	}

	public double DerivSinfunction() {
		double function;

		function = funcparamtwo * funcparamthree * Math.cos(Math.toRadians(realpoints[0] * funcparamthree + funcparamfour));

		return function;
	}
	
	public double Quadfunction() {
		double function;
		
		function = funcparamtwo*realpoints[0]*realpoints[0]+funcparamthree*realpoints[0] + funcparamfour;
		
		return function;
	}

	public double DerivQuadfunction(){
		double function;
		
		function = 2*funcparamtwo*realpoints[0] + funcparamthree;
		
		return function;
	}
	
	public double Cubicfunction(){
		double function;
		
		function = funcparamtwo*Math.pow(realpoints[0], 3) + funcparamthree*Math.pow(realpoints[0], 2) 
		         + funcparamfour*realpoints[0] + funcparamfive;
		
		return function;
	}
	
	public double DerivCubicfunction(){
		double function;
		
		function = 3*funcparamtwo*Math.pow(realpoints[0], 2) + 2*funcparamthree*realpoints[0] + funcparamfour;
		
		return function;
	}
	
	public double Biquadfunction(){
		double function;
		
		function = funcparamtwo*Math.pow(realpoints[0], 4) + funcparamthree*Math.pow(realpoints[0], 3) 
		         + funcparamfour*Math.pow(realpoints[0], 2) + funcparamfive*realpoints[0] +funcparamsix ;
		
		return function;
	}
	
	public double DerivBiquadfunction(){
		double function;
		
		function = 4*funcparamtwo*Math.pow(realpoints[0], 3) + 3*funcparamthree*Math.pow(realpoints[0], 2) 
		         + 2*funcparamfour*realpoints[0] + funcparamfive;
		
		return function;
	}
	
	

}

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

		function = funcparamtwo * Math.sin(Math.toRadians(realpoints[0] + funcparamthree));

		return function;
	}
	
	public double DerivSinfunction() {
		double function;

		function = funcparamtwo  * Math.cos(Math.toRadians(realpoints[0] + funcparamthree));

		return function;
	}
	
    public double Gaussianfunction(){
    	double function;
    	
    	function = funcparamtwo*Math.exp(-realpoints[0]*realpoints[0]/funcparamthree);
    	
    	return function;
    }
    
    public double DerivGaussianfunction(){
    	double function;
    	
    	function = funcparamtwo*(-2*realpoints[0]/(Math.pow(funcparamthree, 2)))*Math.exp(-realpoints[0]*realpoints[0]/funcparamthree);
    	
    	return function;
    }
	
	public double HalfCirclefunction(){
		double function;
		
		function = Math.sqrt((funcparamtwo*funcparamtwo-Math.pow(realpoints[0]-funcparamone[0],2)))+funcparamone[1];
		
		return function;
	}
	
	public double DerivHalfCirclefunction(){
		double function;
		
		function = -(realpoints[0]-funcparamone[0])/Math.sqrt((funcparamtwo*funcparamtwo-
				Math.pow(realpoints[0]-funcparamone[0],2)));
		return function;
	}
	
	public double OtherHalfCirclefunction(){
		double function;
		
		function = -Math.sqrt(funcparamtwo*funcparamtwo-Math.pow(realpoints[0]-funcparamone[0],2))+funcparamone[1];
		
		return function;
	}
	
	public double DerivOtherHalfCirclefunction(){
		double function;
		
		function = (realpoints[0]-funcparamone[0])/Math.sqrt(funcparamtwo*funcparamtwo-Math.pow(realpoints[0]-funcparamone[0],2));
		return function;
	}
	

	
	
	public double Quadfunction() {
		double function;
		
		function = funcparamtwo*Math.pow(realpoints[0], 2)+funcparamthree*realpoints[0] + funcparamfour;
		
		return function;
	}

	public double DerivQuadfunction(){
		double function;
		
		function = 2*funcparamtwo*realpoints[0] + funcparamthree;
		
		return function;
	}
	public double DerivsecQuadfunction(){
		double function;
		
		function = 2*funcparamtwo;
		
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
	
	public double DerivsecCubicfunction(){
		double function;
		
		function = 6*funcparamtwo*realpoints[0] + 2*funcparamthree;
		
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
	
	
	
	public double Mysin(double x){
		double function;
		
		function = x - Math.pow(x, 3)/6 + Math.pow(x, 5)/120 -Math.pow(x, 7)/5040;  
		
		return function;
	}
	
	public double Mycos(double x){
		double function;
		
		function = 1 - Math.pow(x, 2)/2 + Math.pow(x, 4)/24 - Math.pow(x, 6)/720;
		
		return function;
	}
	/** Quadratic function, x=secondrealpos, A*x^2+B*x+C **/
/*	
	Finalfunction Quadfunction = new Finalfunction(secondrealpos, 0.01,0.02, 4);
	final double functionvalue = Quadfunction.Quadfunction();
	final double functionderiv = Quadfunction.DerivQuadfunction();
	
*/
	/** Cubic function, x=secondrealpos, A*x^3+B*x^2+C*x+D **/
/*	
	Finalfunction Cubicfunction = new Finalfunction(secondrealpos, 1.1,0.1,-4.1, 4);
	final double functionvalue = Cubicfunction.Cubicfunction();
	final double functionderiv = Cubicfunction.DerivCubicfunction();
	
	/** Biquad function, x=secondrealpos, A*x^4+B*x^3+C*x^2+D*x+E **/
/*	
	Finalfunction Biquadfunction = new Finalfunction(secondrealpos,2.1,-2.4,-2.4,1.8,4.0);
	final double functionvalue = Biquadfunction.Biquadfunction();
	final double functionderiv = Biquadfunction.DerivBiquadfunction();
*/	
	

}

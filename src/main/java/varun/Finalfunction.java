package varun;

public class Finalfunction extends Finaldistance {

	public Finalfunction(double[] realpoints, double[] funcparamone, double funcparamtwo, double funcparamthree) {
		
		super(realpoints, funcparamone, funcparamtwo, funcparamthree);
		
	}
	
public double Circlefunction(double[] realpoints, double[] funcparamone, double funcparamtwo){
	
	double function;
	
	function = Math.sqrt(funcparamtwo*funcparamtwo - Math.pow((realpoints[0] - funcparamone[0]),2)) +funcparamone[1] ;
			
	return function;
}
	
	
public double Quadfunction(double[] realpoints, double funcparamtwo, double funcparamthree){
	
	double function;
	
	function = funcparamtwo*realpoints[0]*realpoints[0] +funcparamthree;
	
	
	return function;
}
	
public double Quadfunctionderiv(double[] realpoints, double funcparamtwo){
	
	double function;
	
	function = 2*funcparamtwo*realpoints[0];
	
	
	return function;
}

public double Newton (Finalfunction function, Finalfunction derivfunction, double initialguess, double epsilon, int Maxiter){
	
	
	
	double[] guess = {initialguess,0};
	for (int j= 0; j < Maxiter; j++) { 
	double fval = function.Quadfunction(guess, function.funcparamtwo, function.funcparamthree);
	double fder = derivfunction.Quadfunctionderiv(guess, derivfunction.funcparamtwo);
	double dx = fval/fder;
	initialguess -= dx;
	guess[0] = initialguess;
	guess[1] = 0;
	
	if (Math.abs(dx) < epsilon)
		 return guess[0]; 
	
	}
	return guess[0];
}







}

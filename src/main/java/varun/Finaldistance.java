package varun;

import net.imglib2.Cursor;
import net.imglib2.Localizable;
import net.imglib2.type.numeric.RealType;

public abstract class Finaldistance {

	double[] realpoints;
	double[] funcparamone;
	double funcparamtwo;
	double funcparamthree;
	double funcparamfour;
	double funcparamfive;
	double funcparamsix;

	public Finaldistance(double[] realpoints, double[] funcparamone, double funcparamtwo, double funcparamthree) {

		this.realpoints = realpoints;
		this.funcparamone = funcparamone;
		this.funcparamtwo = funcparamtwo;
		this.funcparamthree = funcparamthree;

	}

	public Finaldistance(double[] realpoints, double funcparamtwo, double funcparamthree) {

		this.realpoints = realpoints;
		this.funcparamtwo = funcparamtwo;
		this.funcparamthree = funcparamthree;

	}

	public Finaldistance(double[] realpoints, double funcparamtwo, double funcparamthree, double funcparamfour) {

		this.realpoints = realpoints;
		this.funcparamtwo = funcparamtwo;
		this.funcparamthree = funcparamthree;
		this.funcparamfour = funcparamfour;

	}
	
	public Finaldistance(double[] realpoints, double funcparamtwo, double funcparamthree, double funcparamfour, double funcparamfive) {

		this.realpoints = realpoints;
		this.funcparamtwo = funcparamtwo;
		this.funcparamthree = funcparamthree;
		this.funcparamfour = funcparamfour;
		this.funcparamfive = funcparamfive;

	}
	
	public Finaldistance(double[] realpoints, double funcparamtwo, double funcparamthree, double funcparamfour, 
			double funcparamfive, double funcparamsix) {

		this.realpoints = realpoints;
		this.funcparamtwo = funcparamtwo;
		this.funcparamthree = funcparamthree;
		this.funcparamfour = funcparamfour;
		this.funcparamfive = funcparamfive;
		this.funcparamsix = funcparamsix;

	}
	public static double Generalfunctiondist(double[] secondpos, double[] firstpos) {
		double distance ;
	
		distance = Math.pow((secondpos[0] - firstpos[0]), 2) + Math.pow((secondpos[1] - firstpos[1]), 2);

		return Math.sqrt(distance);
	}	

	public static <T extends RealType<T>> double  disttocurve(double[] secondrealpos,double[] realpos,
			double functionvalue, double functionderiv){
            //x= secondrealpos[0], f(x) = function value or = y for implicit curves
		Finalfunction Normalline = new Finalfunction(realpos, -1.0 / functionderiv,
				secondrealpos[0] / functionderiv + functionvalue);
	final double	 distanceline = Normalline.Linefunctiondist();
	
	return distanceline;
	}
	
	public static <T extends RealType<T>> double  disttocurvetangent(double[] secondrealpos,double[] realpos,
			double functionvalue, double functionderiv){
		Finalfunction Tangentline = new Finalfunction(realpos,  functionderiv,
				 functionvalue- secondrealpos[0]*functionderiv);
	final double	 distanceline = Tangentline.Linefunctiondist();
	
	return distanceline;
	}
	
}

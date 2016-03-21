package varun;

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
		double distance = 0;

		distance = Math.pow((secondpos[0] - firstpos[0]), 2) + Math.pow((secondpos[1] - firstpos[1]), 2);

		return Math.sqrt(distance);
	}	

}

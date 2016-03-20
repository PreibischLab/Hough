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

	public static double Generalfunctiondist(double[] secondpos, double[] firstpos) {
		double distance = 0;

		distance = Math.pow((secondpos[0] - firstpos[0]), 2) + Math.pow((secondpos[1] - firstpos[1]), 2);

		return Math.sqrt(distance);
	}

}

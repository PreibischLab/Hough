package varun;



public class ButterflyFilter {

	public static double[][] createButterflyKernel(double dtheta, double drho, double length){
		int size = 3;
		
		double[][] butterflyKernel = new double[size][size];
		double l = drho/Math.sin(Math.toRadians(dtheta));
        butterflyKernel[0][0]  = 0;
        butterflyKernel[0][1] = -(2*l+length)/2;
        butterflyKernel[0][2] = 0;
        butterflyKernel[1][0]  = l;
        butterflyKernel[1][1] = length;
        butterflyKernel[1][2] = l;
        butterflyKernel[2][0]  = 0;
        butterflyKernel[2][1] = -(2*l+length)/2;
        butterflyKernel[2][2] = 0;
		
		return butterflyKernel;
	}
	
}
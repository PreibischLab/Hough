package peakFitter;

import mISC.Tree.Distance;
import net.imglib2.Point;

public final class GaussianLine implements MTFitFunction {

	@Override
	public double val(double[] x, double[] a, double[] b) {
		final int ndims = x.length;
		return a[2 * ndims + 1] * Etotal(x, a, b);
	}

	@Override
	public double grad(double[] x, double[] a, double[] b, int k) {
		final int ndims = x.length;
		final double slope = b[ndims];

		if (k < ndims) {

			return 2 * b[k] * (x[k] - a[k]) * a[2 * ndims + 1] * Estart(x, a, b) ;

		}

		

		else if (k >= ndims  && k <=  ndims + 1) {
			int dim = k - ndims ;
			return 2 * b[dim] * (x[dim] - a[k]) * a[2 * ndims + 1] * Eend(x, a, b);

		}
		
        if (k == 2 * ndims) 

        	return Ederivds(x, a, b);
		
        
		else if (k == 2 * ndims + 1)
			
			return Etotal(x, a, b);

		else
			return 0;

	}

	/*
	 * PRIVATE METHODS
	 */

	/*
	 * @ Define a line analytically as a sum of gaussians, the parameters to be
	 * determined are the start and the end points of the line
	 * 
	 */

	private static final double Estart(final double[] x, final double[] a, final double[] b) {

		double sum = 0;
		double di;
		for (int i = 0; i < x.length; i++) {
			di = x[i] - a[i];
			sum += b[i] * di * di;
		}

		return Math.exp(-sum);

	}

	private static final double Eend(final double[] x, final double[] a, final double[] b) {

		double sum = 0;
		double di;
		int ndims = x.length;
		for (int i = 0; i < x.length; i++) {

			di = x[i] - a[i + ndims];
			sum += b[i] * di * di;
		}

		return Math.exp(-sum);

	}

	private static final double Etotal(final double[] x, final double[] a, final double[] b) {

		return Estart(x, a, b) + Eend(x, a, b) + Esum(x, a, b);

	}

	private static final double Esum(final double[] x, final double[] a, final double[] b) {

		final int ndims = x.length;
		final double slope = b[ndims];
		double[] minVal = new double[ndims];
		double[] maxVal = new double[ndims];

		if (slope >= 0) {
			for (int i = 0; i < x.length; i++) {
				minVal[i] = a[i];
				maxVal[i] = a[ndims + i];
			}

		}

		if (slope < 0) {
			minVal[0] = a[0];
			minVal[1] = a[ndims + 1];
			maxVal[0] = a[ndims];
			maxVal[1] = a[1];

		}

		double sum = 0;
		double sumofgaussians = 0;
		double di;

		double[] steppos = new double[ndims];
		double dsstart = a[2 * ndims];
		double dxstart = dsstart/ Math.sqrt( 1 + slope * slope) ;
		double dystart = slope * dxstart;

		steppos[0] = minVal[0];
		steppos[1] = minVal[1];

		while (true) {

			steppos[0] += dxstart;
			steppos[1] += dystart;
			sum = 0;
			for (int i = 0; i < x.length; i++) {
				di = x[i] - steppos[i];
				sum += b[i] * di * di;
			}
			sumofgaussians += Math.exp(-sum);

			if (steppos[0] >= maxVal[0] || steppos[1] >= maxVal[1] )
				break;
		}

		return sumofgaussians;
	}

	

	private static final double Ederivds(final double[] x, final double[] a, final double[] b) {

		final int ndims = x.length;
		final double slope = b[ndims];
		double[] minVal = new double[ndims];
		double[] maxVal = new double[ndims];

		if (slope >= 0) {
			for (int i = 0; i < x.length; i++) {
				minVal[i] = a[i];
				maxVal[i] = a[ndims + i];
			}

		}

		if (slope < 0) {
			minVal[0] = a[0];
			minVal[1] = a[ndims + 1];
			maxVal[0] = a[ndims];
			maxVal[1] = a[1];

		}

	
		double sumofgaussians = 0;
		double di, de;

     		double[] steppos = new double[ndims];
     		double[] endpos = new double[ndims];
     		double dsstart = a[2 * ndims];
    		double dxstart = dsstart/ Math.sqrt( 1 + slope * slope) ;
    		double dystart = slope * dxstart;

    		steppos[0] = minVal[0];
    		steppos[1] = minVal[1];
    		
    		endpos[0] = maxVal[0];
    		endpos[1] = maxVal[1];
    		
    		while (true){
    			steppos[0] += dxstart;
    			steppos[1] += dystart;
    			endpos[0] -= dxstart;
    			endpos[1] -= dystart;
    			
    			double sum = 0;
    			double dsum = 0;
    			double esum = 0 ;
    			double endsum = 0;
    			for (int i = 0; i < x.length; i++) {
    				di = x[i] - steppos[i];
    				sum += b[i] * di * di;
    				de = x[i] - endpos[i]; 
    				endsum += b[i] * de * de;
    				
    				if (i == 0)
    					dsum += 2 * b[i] * di/ Math.sqrt( 1 + slope * slope)  ;
    				if (i == 1)
    					dsum += 2  * b[i] * di * slope/ Math.sqrt( 1 + slope * slope) ;
    				
    				if (i == 0)
    					esum += -2  * b[i] * de/ Math.sqrt( 1 + slope * slope)  ;
    				if (i == 1)
    					esum += -2 *  b[i] * de * slope/ Math.sqrt( 1 + slope * slope) ;
    				
    			}
    			sumofgaussians += 0.5*(dsum * Math.exp(-sum) + esum * Math.exp(-endsum)) ;

    			if (steppos[0] >= maxVal[0] || steppos[1] >= maxVal[1] || endpos[0] <= minVal[0] || endpos[1] <= minVal[1]   )
    				break;
    			
    		}
     		
			

		return  sumofgaussians ;
	}

	
	
}

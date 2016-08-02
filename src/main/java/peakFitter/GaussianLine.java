package peakFitter;

import net.imglib2.Point;

public final class GaussianLine implements MTFitFunction {

	@Override
	public double val(double[] x, double[] a, double[] b) {
		final int ndims = x.length;
		return a[2*ndims]*E(x,a,b);
	}

	@Override
	public double grad(double[] x, double[] a, double[] b, int k) {
		final int ndims = x.length;
		if (k < ndims  ){
	
		return 2 * b[k] * (x[k] - a[k])*a[2*ndims] * E(x, a, b);
		}
		
		else if (k >= ndims && k < 2*ndims ){
			int dim = k - ndims;
			return 2 * b[dim] * (x[dim] - a[k])*a[2*ndims] * E(x, a, b);
			
		}
		
		else if (k == 2*ndims )
			return E(x,a,b);
		
		else
			return 0;

	}
	
	
	
	/*
	 * PRIVATE METHODS
	 */

	private static final double E(final double[] x, final double[] a, final double[] b){
		
		final int ndims = x.length;
		final double slope = b[ndims];
		 double[] minVal = new double[ndims];
		 double[] maxVal = new double[ndims];
		
		
			for (int i = 0; i < x.length; i++) {
				minVal[i] = a[i];
				maxVal[i] = a[ndims + i];
			}
		
		
		
		final double intercept = b[ndims + 1];
		double sum = 0;
		double sumofgaussians = 0;
		double di;
		
		
		
		
		final double stepsize =  Math.sqrt(b[0] * b[0] + b[1] * b[1]);
		final double[] steppos = new double[ndims];
		int count = 0;
		if (slope >= 0){
		while (true) {
			
			steppos[0] = minVal[0] + count * stepsize / Math.sqrt(1 + slope * slope);
			steppos[1] = minVal[1] + count * stepsize * slope / Math.sqrt(1 + slope * slope);
			sum = 0;
			for (int i = 0; i < x.length; i++) {
				di = x[i] - steppos[i];
				sum += b[i] * di * di;
			}
			sumofgaussians +=Math.exp(-sum);
			
			count++;
			if (steppos[0] > maxVal[0] || steppos[1] > maxVal[1])
				break;
		}
		
		
		}
		
		int negcount = 0;
		if (slope < 0){
			while (true) {
				
				steppos[0] = minVal[0] + negcount * stepsize / Math.sqrt(1 + slope * slope);
				steppos[1] = maxVal[1] + negcount * stepsize * slope / Math.sqrt(1 + slope * slope);
				sum = 0;
				for (int i = 0; i < x.length; i++) {
					di = x[i] - steppos[i];
					sum += b[i] * di * di;
				}
				sumofgaussians +=Math.exp(-sum);
				

				negcount++;
				if (steppos[0] > maxVal[0] || steppos[1] < minVal[1])
					break;
			}
			
			
			
		}
		
		
		
		return sumofgaussians;
	}
	

}

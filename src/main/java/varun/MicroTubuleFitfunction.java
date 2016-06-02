package varun;



// For each point along the HT determined line, the fit function is modelled as a gaussian

public  class MicroTubuleFitfunction implements MTFitFunction {

	public final double val(final double[] x, final double[] a) {
		
		
		
			return a[0] * E(x, a);
		
		

	}

	@Override
	public double grad(double[] x, double[] a, int k) {
		
		final int ndims = x.length;
		if (k == 0) {
			// With respect to A
			return E(x, a);

		} else if (k <= ndims) {
			// With respect to xi
			int dim = k - 1;
			return 2 * a[dim+ndims] * (x[dim] - a[dim+1]) * a[0] * E(x, a);

		} else {
			// With respect to ai
			int dim = k - ndims - 1;
			double di = x[dim] - a[dim+1];
			return - di * di * a[0] * E(x, a);
		}
		
		
	}

	/*
	 * PRIVATE METHODS 
	 */

	private static final double E(final double[] x, final double[] a) {
		final int ndims = x.length;
		double sum = 0;
		double di;
		
			
		for (int i = 0; i < x.length; i++) {
			di = x[i] - a[i + 1];
			sum += a[i + ndims + 1] * di * di;
		}
		
		return Math.exp(sum);
	
		
	}
}


package peakFitter;

public class MTfitDLM implements MTFitFunction {
	

	
		 
	
			/*
			 * METHODS
			 */

			@Override
			public final double val(final double[] x, final double a, final double[] fixedparam) {
				return a * E(x, fixedparam);
			}

			/**
			 * Partial derivatives indices are ordered as follow:
			 * <pre>k = 0       - A
			 *k = 1..n    - x_i (with i = k-1)
			 *k = n+1..2n - b_i (with i = k-n-1)</pre> 
			 */
			@Override
			public final double grad(final double[] x, final double a, final double[]  fixedparam, final int k) {
				final int ndims = x.length;
				if (k == 0) {
					// With respect to A
					return E(x, fixedparam);

				} else if (k <= ndims) {
					// With respect to xi
					int dim = k - 1;
					return 2 * fixedparam[dim+ndims] * (x[dim] - fixedparam[dim+1]) * a * E(x, fixedparam);

				} else {
					// With respect to ai
					int dim = k - ndims - 1;
					double di = x[dim] - fixedparam[dim+1];
					return - di * di * a * E(x, fixedparam);
				}
			}

			
			/*
			 * PRIVATE METHODS
			 */

			private static final double E(final double[] x, final double[] fixedparam) {
				final int ndims = x.length;
				double sum = 0;
				double di;
				for (int i = 0; i < x.length; i++) {
					di = x[i] - fixedparam[i+1];
					sum += fixedparam[i+ndims+1] * di * di;
				}
				return Math.exp(-sum);
			}
			

		}


		
		
		




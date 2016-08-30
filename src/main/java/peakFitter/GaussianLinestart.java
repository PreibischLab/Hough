package peakFitter;

public final class GaussianLinestart implements MTFitFunction {
	
	


		@Override
		public double val(double[] x, double[] a, double[] b) {
			final int ndims = x.length;
			return a[2 * ndims + 1 ] * E(x, a, b);
		}

		@Override
		public double grad(double[] x, double[] a, double[] b, int k) {
			final int ndims = x.length;
			final double slope = b[ndims];
			
			
			
			if (k < ndims) {
	       return  2 * b[k]* (x[k] - a[k]) * a[2 * ndims + 1] * E(x, a, b);
				
			}

			if (k == ndims){
			 int dim = k -  ndims ;
			 
			 double dx = a[k]/ Math.sqrt(1 + slope *slope);
			 double dy = slope * dx;
			 
			return -2 * b[dim] * (x[dim] - (a[dim] - dx)) * a[2 * ndims + 1]
				    * E(x, a, b)/ Math.sqrt(1 + slope *slope)
				   - 2 * b[dim + 1] * (x[dim + 1] - (a[dim + 1] - dy)) * a[2 * ndims + 1]
					*slope* E(x, a, b)/ Math.sqrt(1 + slope *slope) ;
			}
			else if (k >= ndims + 1 && k <= 2 * ndims) {
				int dim = k - ndims - 1;
				return 2 * b[dim] * (x[dim] - a[k]) * a[2 * ndims + 1] * E(x, a, b);

			}
			
			
			
			

			
			else if (k == 2 * ndims + 1)
				return E(x, a, b);
		
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

		private static final double E(final double[] x, final double[] a, final double[] b) {

			final int ndims = x.length;
			final double slope = b[ndims];
			double[] minVal = new double[ndims];
			double[] maxVal = new double[ndims];

			if (slope >= 0) {
				for (int i = 0; i < x.length; i++) {
					minVal[i] = a[i];
					maxVal[i] = a[ndims + i + 1];
				}

			}

			if (slope < 0) {
				minVal[0] = a[0];
				minVal[1] = a[ndims + 2];
				maxVal[0] = a[ndims + 1];
				maxVal[1] = a[1];

			}

			double sum = 0;
			double sumofgaussians = 0;
			double di;
			

			
			double[] steppos = new double[ndims];
			double[] endpos = new double[ndims];
			
			
			 double dsstart = a[ndims];
			 double dxstart = dsstart / Math.sqrt(1 + slope * slope);
			 double dystart = slope * dxstart;
			 
		
			    steppos[0] = minVal[0] ;
				steppos[1] = minVal[1]  ;
				for (int i = 0; i < x.length; i++) {
					di = x[i] - steppos[i];
					sum += b[i] * di * di;
				}
				sumofgaussians += Math.exp(-sum);
				
				steppos[0] += dxstart ;
				steppos[1] += dystart ;
				sum = 0;
				for (int i = 0; i < x.length; i++) {
					di = x[i] - steppos[i];
					sum += b[i] * di * di;
				}
				sumofgaussians += Math.exp(-sum);
				
				endpos[0] = maxVal[0];
				endpos[1] = maxVal[1];
				sum = 0;
				for (int i = 0; i < x.length; i++) {
					di = x[i] - endpos[i];
					sum += b[i] * di * di;
				}
				sumofgaussians += Math.exp(-sum);
				
				while (true) {

					steppos[0] += dxstart ;
					steppos[1] += dystart ;
					sum = 0;
					for (int i = 0; i < x.length; i++) {
						di = x[i] - steppos[i];
						sum += b[i] * di * di;
					}
					sumofgaussians += Math.exp(-sum);

					if (steppos[0] >= maxVal[0]  || steppos[1] >= maxVal[1] )
						break;
				}
				
			


			return sumofgaussians;
		}
		
		
		
		

	}




package peakFitter;


	

	public class GaussianLinemid implements MTFitFunction {
		@Override
		public double val(double[] x, double[] a, double[] b) {
			final int ndims = x.length;
			return a[2 * ndims] * Etotal(x, a, b);
		}

		@Override
		public double grad(double[] x, double[] a, double[] b, int k) {
			final int ndims = x.length;

			if (k < ndims) {

				return 2 * b[k] * (x[k] - a[k]) *a[2 * ndims] *  (Estart(x, a, b) ) ;

			}

			

			else if (k >= ndims  && k <=  ndims + 1) {
				int dim = k - ndims ;
				return 2 * b[dim] * (x[dim] - a[k]) *a[2 * ndims] *  (Eend(x, a, b) );

			}
			
	        

			else if (k == 2 * ndims)
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
			double[] minVal = new double[ndims];
			double[] maxVal = new double[ndims];

			
				for (int i = 0; i < x.length; i++) {
					minVal[i] = a[i];
					maxVal[i] = a[ndims + i];
				}

			

			double sum = 0;
			double sumofgaussians = 0;
			double di;

			double[] steppos = new double[ndims];
			
			double dist = Distance(maxVal, minVal);
			
			double distmid = 0;
			double distend = 0;
		
			
			
			steppos[0] = minVal[0];
			steppos[1] = minVal[1];

			double[] endpos = new double[ndims];
			
			endpos[0] = maxVal[0];
			endpos[1] = maxVal[1];
			
			
			double[] midpos = {(endpos[0] + steppos[0]) * 0.5,(endpos[1] + steppos[1]) * 0.5 };
			double[] newmidpos ={(endpos[0] + steppos[0]) * 0.5,(endpos[1] + steppos[1]) * 0.5 };
		
			while (true){
				
				for (int d  = 0; d < ndims; ++d){
					midpos[d] = 0.5*(steppos[d] + midpos[d]);
					newmidpos[d] = 0.5*(endpos[d] + newmidpos[d]);
				}
				
		
				distend = Distance(newmidpos,endpos);
				distmid = Distance(midpos, steppos);
				
				if (Math.abs(distmid - distend) < 1.0e-10){
					
					break;
				}
			}
			
			final double ds = distmid/2;
			double slope = (maxVal[1] - minVal[1]) / (maxVal[0] - minVal[0]);
			double dxstart = ds / Math.sqrt(1 + slope * slope);
			double dystart = slope * dxstart;
			
			while (true) {

				steppos[0] += dxstart ;
				steppos[1] += dystart ;
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

		

		public static double Distance(final double[] cordone, final double[] cordtwo) {

			double distance = 0;
			final double ndims = cordone.length;

			for (int d = 0; d < ndims; ++d) {

				distance += Math.pow((cordone[d] - cordtwo[d]), 2);

			}
			return Math.sqrt(distance);
		}
		
	}



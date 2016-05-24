package varun;
import java.util.ArrayList;

import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import varun.OverlayLines.Simulatedline;
import varun.PerformWatershedding.Lineobjects;
public class LengthDetection {

	
	
	/**
	 * Try to build a best guess parameter set to serve as a starting point for optimization.
	 * <p>
	 * The principle of this guess relies on supposing the intensity distribution is 
	 * actually a multivariate normal distribution from which numbers are 
	 * drawn, with a frequency I. We try to retrieve parameters from the distribution
	 * using a maximum likelihood approach. For the normal distribution, parameter
	 * estimates can be derived analytically, using the covariance matrix.  
	 * 
	 * @param X  the coordinates of observations
	 * @param I  the pixel value of observations
	 * @return  a double array containing crude estimates for parameters in this order: 
	 * <pre> 0.			A
	 * 1 → ndims		x₀ᵢ
	 * ndims+1 → 2 × ndims	cᵢ = 1 / σᵢ² </pre>
	 */
	private final double[] makeBestGuess(
			final ArrayList<Simulatedline> linelist,
			final double[][] X, 
			final double[] I) {
		
		int ndims = I.length;
		
		
		final double[][] Xcor;
		for (int index = 0; index < linelist.size(); ++index){
			 double [] position = new double[ndims];
			 FloatType Intensity;
			 position = linelist.get(index).point;
			 Intensity = linelist.get(index).Value;
			
		}
		
		
		double[] start_param = new double[2*ndims+1];

		double[] X_sum = new double[ndims];
		for (int j = 0; j < ndims; j++) {
			X_sum[j] = 0;
			for (int i = 0; i < X.length; i++) {
				X_sum[j] += X[i][j] * I[i];
			}
		}

		double I_sum = 0;
		double max_I = Double.NEGATIVE_INFINITY;
		for (int i = 0; i < X.length; i++) {
			I_sum += I[i];
			if (I[i] > max_I) {
				max_I = I[i];
			}

		}

		start_param[0] = max_I;

		for (int j = 0; j < ndims; j++) {
			start_param[j+1] = X_sum[j] / I_sum;
		}

		for (int j = 0; j < ndims; j++) {
			double C = 0;
			double dx;
			for (int i = 0; i < X.length; i++) {
				dx = X[i][j] - start_param[j+1];
				C += I[i] * dx * dx;
			}
			C /= I_sum;
			start_param[ndims + j + 1] = 1 / C;
		}
		
		return start_param;		
	}
	
}

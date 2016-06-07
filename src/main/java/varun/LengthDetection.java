package varun;

import java.util.ArrayList;

import com.sun.tools.javac.util.Pair;

import net.imglib2.Cursor;
import net.imglib2.Interval;
import net.imglib2.Localizable;
import net.imglib2.Point;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.kdtree.HyperPlane;
import net.imglib2.algorithm.localization.LevenbergMarquardtSolver;
import net.imglib2.algorithm.neighborhood.DiamondShape.NeighborhoodsIterableInterval;
import net.imglib2.algorithm.neighborhood.Neighborhood;
import net.imglib2.algorithm.neighborhood.RectangleNeighborhood;
import net.imglib2.algorithm.neighborhood.RectangleShape;
import net.imglib2.algorithm.region.hypersphere.HyperSphere;
import net.imglib2.algorithm.region.hypersphere.HyperSphereCursor;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.roi.RectangleRegionOfInterest;
import net.imglib2.roi.RegionOfInterest;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Intervals;
import net.imglib2.view.Views;
import varun.OverlayLines.Simulatedline;
import varun.PerformWatershedding.Lineobjects;

public class LengthDetection {
	public static final class Labelparam {
		final int Label;
		final double[] point;
		final FloatType Value;
		final double slope;
		final double intercept;
		

		protected Labelparam(final int Label, final double[] point, final FloatType Value, final double slope,
				final double intercept) {
			this.Label = Label;
			this.point = point;
			this.Value = Value;
			this.slope = slope;
			this.intercept = intercept;

		}
	}

	

	private final RandomAccessibleInterval<FloatType> inputimg;
	private final RandomAccessibleInterval<IntType> intimg;
	private final int ndims;
	
	public LengthDetection(RandomAccessibleInterval<FloatType> inputimg,RandomAccessibleInterval<IntType> intimg){
	
		this.inputimg = inputimg;
        this.intimg = intimg;
		this.ndims = inputimg.numDimensions();
		
		
	}

	

	/**
	 * Input a guess list containing the list of centroids determined by the HT,
	 * and read from the PointSamplelist of local maximas the psf of the
	 * microscope and
	 * 
	 * @returns a doube array startparam[] containing guess for the Gaussian at
	 *          one centroid point
	 * 
	 *          start_param[0] = Amplitude of Gaussian start_param[1 → ndims] = (x_i, y_i) (by HT) Centroids of Gaussian
	 *          start_param[ndims+1 → 2 × ndims] = 1 / σ_i^2 Sigmas of the Gaussian
	 * 
	 */

	public double[] makeBestGuess(final Localizable point, final Float Value, final double[][] X, final double[] I) {
		double[] start_param = new double[2 * ndims + 1];

		
		for (int d = 0; d < ndims; ++d) {

			start_param[d + 1] = point.getDoublePosition(d);
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

	
	// Get final parameters for the Gaussian fit along the line, input the
	// centroid position and value where the Gaussian has to be fitted
	// in the image
	public double[] Getfinalparam(final Localizable point, final Float Value, final double[] typical_sigma) throws Exception {

		

		RandomAccess<IntType> ranac = intimg.randomAccess();
		
		ranac.setPosition(point);
		
		int label = ranac.get().get();
		
		
		 PointSampleList<FloatType> datalist = gatherData(point, label);
		
		
		
		final Cursor<FloatType> listcursor = datalist.localizingCursor();
		
		
		
		double[][] X 	= new double[(int) datalist.size()][ndims];
		double[] I 		= new double[(int) datalist.size()];
		int index = 0;
		while(listcursor.hasNext()){
			listcursor.fwd();
			
			for (int d = 0; d < ndims; d++) {
				X[index][d] = listcursor.getDoublePosition(d);
			}
			
			I[index] = listcursor.get().getRealDouble();
			
			index++;
		}

		final double[] start_param = makeBestGuess( point, Value,X,I);
		// Correct for too large sigmas: we drop estimate and replace it by user input
				for (int j = 0; j < ndims; j++) {
					if ( start_param [ j + ndims + 1 ] <  1 / ( typical_sigma[j] * typical_sigma[j] ) ) {
						start_param [ j + ndims + 1 ] = 1 / ( typical_sigma[j] * typical_sigma[j] );
					}
				}
		final double[] finalparam = start_param.clone();
		int maxiter = 400;
		double lambda = 0.0001;
		double termepsilon = 0.01;

	
			LevenbergMarquardtSolverLocal.solve(X, finalparam, I, new GaussianMultiDLM(), lambda, termepsilon,
					maxiter);
			// NaN protection: we prefer returning the crude estimate that NaN
			for (int j = 0; j < finalparam.length; j++) {
				if (Double.isNaN(finalparam[j]))
					finalparam[j] = start_param[j];
			}
			
			

		return finalparam;

	}

	
	private PointSampleList<FloatType> gatherData(final Localizable point, final int label){
		
		final PointSampleList<FloatType> datalist = new PointSampleList<FloatType>(ndims);
		
		Cursor<IntType> intcursor = Views.iterable(intimg).localizingCursor();
		RandomAccess<FloatType> ranac = inputimg.randomAccess();
		
		while(intcursor.hasNext()){
			intcursor.fwd();
			int i = intcursor.get().get();
			if (i == label) {
				final Point newpoint = new Point(intcursor);
				ranac.setPosition(newpoint);
				datalist.add(newpoint, ranac.get().copy());
			}
		}
		
		return datalist;
	}
	
	
	
	

}

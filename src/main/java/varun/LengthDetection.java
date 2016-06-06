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

	private class Observation {
		public double[] I;
		public double[][] X;
	}

	private final RandomAccessibleInterval<FloatType> inputimg;
	private final RandomAccessibleInterval<IntType> intimg;
	private final ArrayList<Labelparam> guessline;
	private final PointSampleList<FloatType> centroidlist;
	private final int ndims;
	
	public LengthDetection(RandomAccessibleInterval<FloatType> inputimg, RandomAccessibleInterval<IntType> intimg,
			PointSampleList<FloatType> centroidlist){
	
		
		
		this.inputimg = inputimg;
		this.intimg = intimg;
		this.guessline = null;
		this.centroidlist = centroidlist;
		this.ndims = inputimg.numDimensions();
		
		
	}

	public LengthDetection(RandomAccessibleInterval<FloatType> inputimg, RandomAccessibleInterval<IntType> intimg,
			ArrayList<Labelparam> guessline) {
		
		this.inputimg = inputimg;
		this.intimg = intimg;
		this.guessline = guessline;
		this.centroidlist = null;
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
	 *          start_param[0] = A start_param[1 → ndims] = x_i, y_i (by HT)
	 *          start_param[ndims+1 → 2 × ndims] = 1 / σ_i^2
	 * 
	 */

	public double[] makeBestGuess(final double[] psf, final Localizable point, final Float Value) {
		double[] start_param = new double[2 * ndims + 1];

		start_param[0] = Value;
		for (int d = 0; d < ndims; ++d) {

			start_param[d + 1] = point.getDoublePosition(d);
		}

		for (int d = 0; d < ndims; ++d) {

			start_param[ndims + d + 1] = 1 / psf[d];
		}

		return start_param;
	}

	// @Deprecated version with labels to make the best guess

	public double[] makeBestGuess(final double[] psf, int listindex, int label) {
		double[] start_param = new double[2 * ndims + 1];

		for (int d = 0; d < ndims; ++d) {
			if (guessline.get(listindex).Label == label)
				start_param[d + 1] = guessline.get(listindex).point[d];

		}

		if (guessline.get(listindex).Label == label)
			start_param[0] = guessline.get(listindex).Value.get();

		for (int d = 0; d < ndims; ++d) {

			start_param[ndims + d + 1] = 1 / psf[d];
		}

		return start_param;
	}

	// Get final parameters for the Gaussian fit along the line, input the
	// centroid position and value where the Gaussian has to be fitted
	// in the image
	public double[] Getfinalparam(final double[] psf, final Localizable point, final Float Value) throws Exception {

		double[] pad_size = new double[ndims];
		for (int i = 0; i < ndims; i++) {
			pad_size[i] =  (2 * psf[i] + 1);
		}

		//final Observation data = gatherObservationData(point, pad_size);
		final PointSampleList<FloatType> datalist = gatherData(point, pad_size);
		
		final double[] start_param = makeBestGuess(psf, point, Value);
		
		final Cursor<FloatType> listcursor = datalist.localizingCursor();
		
		
		int n_pixels = 1;
		for (int i = 0; i < pad_size.length; i++) {
			n_pixels *= (2 * pad_size[i] + 1);
		}
		
		double[][] X 	= new double[n_pixels][ndims];
		double[] I 		= new double[n_pixels];
		int index = 0;
		while(listcursor.hasNext()){
			listcursor.fwd();
			
			for (int d = 0; d < ndims; d++) {
				X[index][d] = listcursor.getDoublePosition(d);
			}
			
			I[index] = listcursor.get().getRealDouble();
			index++;
			
		}

	

		final double[] finalparam = start_param.clone();
		int maxiter = 400;
		double lambda = 1e-3;
		double termepsilon = 1e-3;

	
			LevenbergMarquardtSolverLocal.solve(X, finalparam, I, new MicroTubuleFitfunction(), lambda, termepsilon,
					maxiter);
		

		return finalparam;

	}

	// @Deprectaed version with labels
	public double[] Getfinalparam(final double[] psf, int listindex, int ndims, int label) {

		// Determine the size of the data to gather
		int[] pad_size = new int[ndims];
		long[] pointValue = new long[ndims];
		for (int i = 0; i < ndims; i++) {
			pad_size[i] = (int) Math.ceil(2 * psf[i]);
		}

		for (int i = 0; i < ndims; i++) {
			pointValue[i] = (long) guessline.get(listindex).point[i];
		}
		// Gather data around peak
		final Point newpoint = new Point(ndims);
		newpoint.setPosition(pointValue);

		final Observation data = gatherObservationData(newpoint, pad_size, guessline.get(listindex).Label);

		final double[] start_param = makeBestGuess(psf, listindex, label);

		final double[][] X = data.X;
		final double[] I = data.I;

		final double[] finalparam = start_param.clone();
		int maxiter = 300;
		double lambda = 1e-3;
		double termepsilon = 1e-1;

		try {
			LevenbergMarquardtSolverLocal.solve(X, finalparam, I, new MicroTubuleFitfunction(), lambda, termepsilon,
					maxiter);
		} catch (Exception e) {
			e.printStackTrace();
		}

		return finalparam;

	}

	// @Deprectaed version with labels
	private Observation gatherObservationData(final Localizable point, final int[] pad_size, int label) {

		int ndims = inputimg.numDimensions();
		int n_pixels = 1;
		for (int i = 0; i < pad_size.length; i++) {
			n_pixels *= (2 * pad_size[i] + 1);
		}

		double[][] tmp_X = new double[n_pixels][ndims];
		double[] tmp_I = new double[n_pixels];
		int index = 0;

		// Initialize position
		int[] offset = new int[ndims];
		int[] size = new int[ndims];
		for (int i = 0; i < offset.length; i++) {
			offset[i] = point.getIntPosition(i) - pad_size[i];
			size[i] = 2 * pad_size[i] + 1;
		}

		long radius = 0;
		for (int i = 0; i < pad_size.length; i++) {
			radius += pad_size[i] * pad_size[i];
		}

		// create a view on the source with this interval

		HyperSphere<FloatType> sphere = new HyperSphere<FloatType>(inputimg, point, (int) Math.sqrt(radius));

		HyperSphereCursor<FloatType> cursorsphere = sphere.localizingCursor();

		RandomAccess<IntType> intranac = intimg.randomAccess();

		boolean outoflabel = false;
		boolean outofbounds = false;
		long[] pos = new long[ndims];
		while (cursorsphere.hasNext()) {

			cursorsphere.fwd();
			cursorsphere.localize(pos);
			for (int d = 0; d < ndims; d++) {
				pos[d] += offset[d];
				if (pos[d] < 0 || pos[d] >= inputimg.dimension(d)) {
					outofbounds = true;
					break;
				}
			}
			if (outofbounds) {
				outofbounds = false;
				continue;
			}

			intranac.setPosition(pos);
			final int intlabel = intranac.get().get();

			if (intlabel != label) {
				outoflabel = true;
				break;
			}

			if (outoflabel) {
				outoflabel = false;
				continue;
			}

			for (int i = 0; i < ndims; i++) {
				tmp_X[index][i] = pos[i];
			}

			tmp_I[index] = cursorsphere.get().getRealDouble();
			index++;
		}

		// Now we possibly resize the arrays, in case we have been too close to
		// the
		// image border.
		double[][] X = null;
		double[] I = null;
		if (index == n_pixels) {
			// Ok, we have gone through the whole square
			X = tmp_X;
			I = tmp_I;
		} else {
			// Re-dimension the arrays
			X = new double[index][ndims];
			I = new double[index];
			System.arraycopy(tmp_X, 0, X, 0, index);
			System.arraycopy(tmp_I, 0, I, 0, index);
		}

		Observation obs = new Observation();
		obs.I = I;
		obs.X = X;
		return obs;
	}

	private PointSampleList<FloatType> gatherData(final Localizable point, final double[] pad_size){
		
		final PointSampleList<FloatType> datalist = new PointSampleList<FloatType>(ndims);
		
		// initialize position
		
		int[] size = new int[ndims];
		
		for (int i = 0; i < ndims; i++) {
		size[i] = (int) (2 * pad_size[i] + 1);
		}
		
		
		int span = 0;
		for (int i = 0; i < ndims; i++) {
			span = size[i]/ndims;
			}
		
		
		RandomAccess<FloatType> ranac = inputimg.randomAccess();
		
		ranac.setPosition(point);
		
		HyperSphere<FloatType> sphere = new HyperSphere<FloatType>(inputimg, point, span);

		HyperSphereCursor<FloatType> cursorsphere = sphere.localizingCursor();
		
		boolean outofbounds = false;
		double[] pos = new double[ndims];
		
		while(cursorsphere.hasNext()){
			cursorsphere.fwd();
			cursorsphere.localize(pos);
			
			for (int d = 0; d < ndims; d++) {
			
				if (pos[d] < 0 || pos[d] >= inputimg.dimension(d)) {
					outofbounds = true;
					break;
				}
			}
			if (outofbounds) {
				outofbounds = false;
				continue;
			}
			
			Point newpoint = new Point(cursorsphere);
			datalist.add(newpoint, cursorsphere.get());
			
		}
		
		return datalist;
	}
	
	
	/**
	 * Gather observation data around the point where the Gaussian has to be
	 * fitted by creating a hypersphere around the chosen point, the radius of
	 * the sphere is calculated by the input pad_size array, which should be the
	 * sigma[] of the psf
	 * 
	 * 
	 **/
	private Observation gatherObservationData(final Localizable point, final int[] pad_size) {

		int ndims = inputimg.numDimensions();
		int n_pixels = 1;
		for (int i = 0; i < pad_size.length; i++) {
			n_pixels *= (2 * pad_size[i] + 1);
		}

		double[][] tmp_X = new double[n_pixels][ndims];
		double[] tmp_I = new double[n_pixels];
		int index = 0;

		// Initialize position
		int[] offset = new int[ndims];
		int[] size = new int[ndims];
		for (int i = 0; i < offset.length; i++) {
			offset[i] = point.getIntPosition(i) - pad_size[i];
			size[i] = 2 * pad_size[i] + 1;
		}

		long radius = 0;
		for (int i = 0; i < pad_size.length; i++) {
			radius += pad_size[i] * pad_size[i];
		}

		// create a view on the source with this interval

		HyperSphere<FloatType> sphere = new HyperSphere<FloatType>(inputimg, point, (int) Math.sqrt(radius));

		HyperSphereCursor<FloatType> cursorsphere = sphere.localizingCursor();

		boolean outoflabel = false;
		boolean outofbounds = false;
		long[] pos = new long[ndims];
		while (cursorsphere.hasNext()) {

			cursorsphere.fwd();
			cursorsphere.localize(pos);
			for (int d = 0; d < ndims; d++) {
				pos[d] += offset[d];
				if (pos[d] < 0 || pos[d] >= inputimg.dimension(d)) {
					outofbounds = true;
					break;
				}
			}
			if (outofbounds) {
				outofbounds = false;
				continue;
			}

			if (outoflabel) {
				outoflabel = false;
				continue;
			}

			for (int i = 0; i < ndims; i++) {
				tmp_X[index][i] = pos[i];
			}

			tmp_I[index] = cursorsphere.get().getRealDouble();
			index++;
		}

		// Now we possibly resize the arrays, in case we have been too close to
		// the
		// image border.
		double[][] X = null;
		double[] I = null;
		if (index == n_pixels) {
			// Ok, we have gone through the whole square
			X = tmp_X;
			I = tmp_I;
		} else {
			// Re-dimension the arrays
			X = new double[index][ndims];
			I = new double[index];
			System.arraycopy(tmp_X, 0, X, 0, index);
			System.arraycopy(tmp_I, 0, I, 0, index);
		}

		Observation obs = new Observation();
		obs.I = I;
		obs.X = X;
		return obs;
	}

}

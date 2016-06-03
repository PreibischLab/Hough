package varun;

import java.util.ArrayList;

import com.sun.tools.javac.util.Pair;

import net.imglib2.Cursor;
import net.imglib2.Interval;
import net.imglib2.Localizable;
import net.imglib2.Point;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.kdtree.HyperPlane;
import net.imglib2.algorithm.localization.LevenbergMarquardtSolver;
import net.imglib2.algorithm.neighborhood.DiamondShape.NeighborhoodsIterableInterval;
import net.imglib2.algorithm.neighborhood.Neighborhood;
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

		protected Labelparam(final int Label, final double[] point, final FloatType Value, final double slope) {
			this.Label = Label;
			this.point = point;
			this.Value = Value;
			this.slope = slope;

		}
	}

	private  class Observation {
		public double[] I;
		public double[][] X;
	}

	private final RandomAccessibleInterval<FloatType> inputimg;
	private final RandomAccessibleInterval<IntType> intimg;
	private final ArrayList<Labelparam> guessline;
	private final int ndims;

	public LengthDetection(
			RandomAccessibleInterval<FloatType> inputimg, 
			RandomAccessibleInterval<IntType> intimg,
			ArrayList<Labelparam> guessline) {
		this.inputimg = inputimg;
		this.intimg = intimg;
		this.guessline = guessline;
		this.ndims = inputimg.numDimensions();
	}

	/**
	 * Input a guess list containing the list of centroids determined by the HT,
	 * the psf of the microscope and
	 * 
	 * @returns a doube array startparam[] containing guess for the Gaussian at
	 *          one centroid point
	 * 
	 *          start_param[0] = A start_param[1 → ndims] = x_i, y_i (by HT)
	 *          start_param[ndims+1 → 2 × ndims] = 1 / σ_i^2
	 * 
	 */

	public double[] makeBestGuess(final double[] psf, int listindex,
			int label) {
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

	public double[] Getfinalparam(final double[] psf, int listindex,
			int ndims, int label) {

		
		
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

		final Observation data = gatherObservationData(newpoint, pad_size,
				guessline.get(listindex).Label);

		final double[] start_param = makeBestGuess(psf, listindex, label);

		final double[][] X = data.X;
		final double[] I = data.I;

		final double[] finalparam = start_param.clone();
		int maxiter = 500;
		double lambda = 1e-2;
		double termepsilon = 1e-1;

		try {
			LevenbergMarquardtSolverLocal.solve(X, finalparam, I, new MicroTubuleFitfunction(), lambda, termepsilon,
					maxiter);
		} catch (Exception e) {
			e.printStackTrace();
		}

		return finalparam;

	}

	/**
	 * Gather observation data around the point where the Gaussian has to be
	 * fitted by creating a hypersphere around the chosen point, the radius of
	 * the sphere is calculated by the input pad_size array, which should be the
	 * sigma[] of the psf
	 * 
	 * 
	 **/

	public  Observation gatherObservationData( final Localizable point, final int[] pad_size, int label) {

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

		HyperSphere<FloatType> sphere = new HyperSphere<FloatType>(inputimg, point, radius);

		HyperSphereCursor<FloatType> cursorsphere = sphere.localizingCursor();

		RandomAccess<IntType> intranac = intimg.randomAccess();

		boolean outoflabel = false;
		boolean outofbounds = false;
		long[] pos = new long[ndims];
		while (cursorsphere.hasNext()) {

			cursorsphere.fwd();
			cursorsphere.localize(pos);

			intranac.setPosition(pos);
			final int intlabel = intranac.get().get();

			if (intlabel != label) {
				outoflabel = true;
				break;
			}
			// Offset position array and check if we are not out of bounds
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
		double[][] X = tmp_X;
		double[] I = tmp_I;

		Observation obs = new Observation();
		obs.I = I;
		obs.X = X;
		return obs;
	}

	public  double[] getExactparam( double[] start_param, double[] typical_sigma, int label) {
		int n = inputimg.numDimensions();

		Cursor<IntType> intCursor = Views.iterable(intimg).localizingCursor();
		RandomAccess<FloatType> ranac = inputimg.randomAccess();
		ArrayList<Simulatedline> listdata = new ArrayList<Simulatedline>();

		while (intCursor.hasNext()) {
			intCursor.fwd();

			int i = intCursor.get().get();
			if (i == label) {

				final double[] position = new double[n];
				ranac.setPosition(intCursor);
				ranac.localize(position);
				final Simulatedline line = new Simulatedline(label, position, ranac.get());
				listdata.add(line);

			}

		}

		double[][] X = new double[listdata.size()][n];
		double[] I = new double[listdata.size()];

		for (int index = 0; index < listdata.size(); ++index) {
			for (int d = 0; d < n; ++d) {
				if (listdata.get(index).Label == label)
					X[index][d] = listdata.get(index).point[d];

			}

		}
		for (int index = 0; index < listdata.size(); ++index) {
			if (listdata.get(index).Label == label)
				I[index] = listdata.get(index).Value.get();
		}

		int maxiter = 500;
		double lambda = 1e-2;
		double termepsilon = 1e-1;

		final double[] finalparam = start_param.clone();

		try {
			LevenbergMarquardtSolverLocal.solve(X, finalparam, I, new MicroTubuleFitfunction(), lambda, termepsilon,
					maxiter);
		} catch (Exception e) {
			e.printStackTrace();
		}

		return finalparam;
	}

}

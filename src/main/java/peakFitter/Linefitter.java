package peakFitter;

import java.util.ArrayList;
import java.util.Random;

import houghandWatershed.Finalfunction;
import houghandWatershed.PerformWatershedding;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.Point;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import peakFitter.GaussianMaskFit.Endfit;
import preProcessing.GetLocalmaxmin;
import sun.tools.tree.ThrowStatement;

public class Linefitter {

	private final RandomAccessibleInterval<FloatType> inputimg;
	private final RandomAccessibleInterval<IntType> intimg;
	private final int ndims;

	public Linefitter(RandomAccessibleInterval<FloatType> inputimg, RandomAccessibleInterval<IntType> intimg) {

		this.inputimg = inputimg;
		this.intimg = intimg;
		this.ndims = inputimg.numDimensions();
		assert inputimg.numDimensions() == intimg.numDimensions();

	}

	public static enum Noise {
		nonoise, noise;

	}

	protected Noise Noise;

	private final double[] MakeLineguess(double slope, double intercept, double[] psf, int label) throws Exception {

		final double[] realpos = new double[ndims];
		double sigmasq, sigma = 1.0;
		sigmasq = sigma * sigma;
		RandomAccessibleInterval<FloatType> imgout = new ArrayImgFactory<FloatType>().create(inputimg, new FloatType());

		final Cursor<FloatType> outcursor = Views.iterable(imgout).localizingCursor();
		long[] newposition = new long[ndims];
		RandomAccess<IntType> ranac = intimg.randomAccess();
		final RandomAccess<FloatType> ranacinput = inputimg.randomAccess();
		double[] minVal = { Double.MAX_VALUE, Double.MAX_VALUE };
		double[] maxVal = { -Double.MIN_VALUE, -Double.MIN_VALUE };
		final double maxintensity = GetLocalmaxmin.computeMaxIntensityinlabel(inputimg, intimg, label);
		while (outcursor.hasNext()) {

			outcursor.fwd();
			outcursor.localize(realpos);
			ranacinput.setPosition(outcursor);
			// To set the pixel intensity as the shortest distance to the curve
			double distance = 0;
			double intensity = 0;

			Finalfunction linefunction = new Finalfunction(realpos, slope, intercept);
			distance = linefunction.Linefunctiondist();

			intensity = (1 / (sigma * Math.sqrt(2 * Math.PI))) * Math.exp(-distance * distance / (2 * sigmasq));

			intensity *= ranacinput.get().get();

			ranac.setPosition(outcursor);
			int i = ranac.get().get();

			if (i == label) {

				outcursor.get().setReal(intensity);
				if (ranacinput.get().get() / maxintensity >= 0.5) {
					outcursor.localize(newposition);

					long pointonline = (long) (newposition[1] - slope * newposition[0] - intercept);

					// To get the min and max co-rodinates along the line so we
					// have starting points to
					// move on the line smoothly

					if (pointonline == 0) {

						for (int d = 0; d < ndims; ++d) {
							if (outcursor.getDoublePosition(d) <= minVal[d])
								minVal[d] = outcursor.getDoublePosition(d);

							if (outcursor.getDoublePosition(d) >= maxVal[d])
								maxVal[d] = outcursor.getDoublePosition(d);

						}

					}

				}
			}
		}

		final double[] MinandMax = new double[2 * ndims + 1];

		if (slope >= 0) {
			for (int d = 0; d < ndims; ++d) {

				MinandMax[d] = minVal[d];
				MinandMax[d + ndims] = maxVal[d];
			}

		}

		if (slope < 0) {

			MinandMax[0] = minVal[0];
			MinandMax[1] = maxVal[1];
			MinandMax[2] = maxVal[0];
			MinandMax[3] = minVal[1];

		}
		MinandMax[2 * ndims] = maxintensity;
		//MinandMax[2 * ndims + 1] = 0.8;
		System.out.println("Label: " + label + " " + "Hough Detection: " + " StartX: " + MinandMax[0] + " StartY: "
				+ MinandMax[1] + " EndX: " + MinandMax[2] + " EndY: " + MinandMax[3]);

		for (int d = 0; d < ndims; ++d) {

			if (MinandMax[d] == Double.MAX_VALUE || MinandMax[d + ndims] == -Double.MIN_VALUE)
				return null;

		}

		return MinandMax;

	}

	// Get line parameters for fitting line to a line in a label

	public double[] Getfinallineparam(final int label, final double slope, final double intercept, final double[] psf,
			double minlength) throws Exception {

		PointSampleList<FloatType> datalist = gatherfullData(label);
		final Cursor<FloatType> listcursor = datalist.localizingCursor();
		double[][] X = new double[(int) datalist.size()][ndims];
		double[] I = new double[(int) datalist.size()];
		int index = 0;
		while (listcursor.hasNext()) {
			listcursor.fwd();

			for (int d = 0; d < ndims; d++) {
				X[index][d] = listcursor.getDoublePosition(d);
			}

			I[index] = listcursor.get().getRealDouble();

			index++;
		}

		final double[] start_param = MakeLineguess(slope, intercept, psf, label);
		if (start_param == null)
			return null;

		else {
			final double[] fixed_param = new double[ndims];

			for (int d = 0; d < ndims; ++d) {

				fixed_param[d] = 1.0 / Math.pow(psf[d], 2);
			}
			
			final double[] finalparamstart = start_param.clone();
			// LM solver part
			int maxiter = 500;
			double lambda = 1e-9;
			double termepsilon = 1e-2;
			final double[] inistartpos = { start_param[0], start_param[1] };
			final double[] iniendpos = { start_param[2], start_param[3] };

			double inicutoffdistance = Distance(inistartpos, iniendpos);

			if (inicutoffdistance > minlength) {

				LevenbergMarquardtSolverLine.solve(X, finalparamstart, fixed_param, I, new GaussianLinesimple(), lambda,
						termepsilon, maxiter);
			
				
				final double[] startpos = { finalparamstart[0], finalparamstart[1] };
				final double[] endpos = { finalparamstart[2], finalparamstart[3] };
				// NaN protection: we prefer returning the crude estimate than
				// NaN
				for (int j = 0; j < finalparamstart.length; j++) {
					if (Double.isNaN(finalparamstart[j]))
						finalparamstart[j] = start_param[j];
				}

				


				int iterations = 500;
/*
				double newslope = (endpos[1] - startpos[1]) / (endpos[0] - startpos[0]);
				double dist = Distance(endpos, startpos);
				double ds = finalparamstart[5];
				
				double dxstart = ds;
				double dystart = newslope * dxstart;
				
			final double maxintensity = finalparamstart[4];
				System.out.println("ds: " + finalparamstart[5] );
			*/
				System.out.println(
						"LM solver : " + " StartX: " + startpos[0] + " StartY:  " + startpos[1] );
				System.out.println("LM solver : " + " EndX: " + endpos[0] + " EndY:  " + endpos[1]);
				System.out.println(" Length:  " + Distance(startpos, endpos));
				double[] startfit = new double[ndims];
				double[] endfit = new double[ndims];
/*
				final double radius = // 2 * Math.sqrt(psf[0] * psf[0] + psf[1] * psf[1]);
						2 * Math.min(psf[0], psf[1]);

				final int numberofgaussians = (int) (radius / ds);
				

				startfit = peakFitter.GaussianMaskFit.sumofgaussianMaskFit(inputimg, intimg, startpos, psf, iterations,
						maxintensity, dxstart, dystart, newslope, numberofgaussians , Endfit.Start, label);

				endfit = peakFitter.GaussianMaskFit.sumofgaussianMaskFit(inputimg, intimg, endpos, psf, iterations,
						maxintensity, dxstart, dystart, newslope, numberofgaussians , Endfit.End, label);
*/
	//			if (startfit == null || endfit == null) {
					startfit = startpos;
					endfit = endpos;
		//		}
			//	System.out.println("Number of gaussians for mask fit:" + (1 + numberofgaussians)  );

				double[] refindedparam = { startfit[0], startfit[1], endfit[0], endfit[1] };

				return refindedparam;

			}

			else

				return null;

		}
	}

	private PointSampleList<FloatType> gatherfullData(final int label) {
		final PointSampleList<FloatType> datalist = new PointSampleList<FloatType>(ndims);

		RandomAccess<IntType> intranac = intimg.randomAccess();

		Cursor<FloatType> localcursor = Views.iterable(inputimg).localizingCursor();

		boolean outofbounds = false;

		while (localcursor.hasNext()) {
			localcursor.fwd();

			for (int d = 0; d < ndims; d++) {

				if (localcursor.getDoublePosition(d) < 0 || localcursor.getDoublePosition(d) >= inputimg.dimension(d)) {
					outofbounds = true;
					break;
				}
			}
			if (outofbounds) {
				outofbounds = false;
				continue;
			}
			intranac.setPosition(localcursor);
			int i = intranac.get().get();
			if (i == label) {

				Point newpoint = new Point(localcursor);
				datalist.add(newpoint, localcursor.get().copy());

			}
		}

		return datalist;
	}

	public double Distance(final double[] cordone, final double[] cordtwo) {

		double distance = 0;

		for (int d = 0; d < ndims; ++d) {

			distance += Math.pow((cordone[d] - cordtwo[d]), 2);

		}
		return Math.sqrt(distance);
	}

}

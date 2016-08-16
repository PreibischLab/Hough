package peakFitter;

import java.util.ArrayList;

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

	private final double[] MakeLineguess(double slope, double intercept, double [] psf, int label) throws Exception {

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

		final double[] MinandMax = new double[2 * ndims + 2];

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

		MinandMax[2 * ndims] = 0.5*(psf[0] + psf[1]) ;
		MinandMax[2 * ndims + 1] = maxintensity;

		for (int d = 0; d < ndims; ++d) {

			if (MinandMax[d] == Double.MAX_VALUE || MinandMax[d + ndims] == Double.MIN_VALUE)
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
			final double[] fixed_param = new double[2 * ndims];

			for (int d = 0; d < ndims; ++d) {

				fixed_param[d] = 1.0 / Math.pow(psf[d], 2);
			}
			fixed_param[ndims] = slope;
			fixed_param[ndims + 1] = intercept;
			final double[] finalparam = start_param.clone();

			// LM solver part
			int maxiter = 750;
			double lambda = 1e-4;
			double termepsilon = 1e-3;

			LevenbergMarquardtSolverLine.solve(X, finalparam, fixed_param, I, new GaussianLine(), lambda, termepsilon,
					maxiter);

			// NaN protection: we prefer returning the crude estimate than NaN
			for (int j = 0; j < finalparam.length; j++) {
				if (Double.isNaN(finalparam[j]))
					finalparam[j] = start_param[j];
			}

			final double[] startpos = { finalparam[0], finalparam[1] };
			final double[] endpos = { finalparam[2], finalparam[3] };

			double maxintensity = finalparam[5];
			double cutoffdistance = Distance(startpos, endpos);

			if (cutoffdistance > minlength) {

				int iterations = 100;

				final double[] pointnearstart = peakFitter.GaussianMaskFit.singlegaussianMaskFit(inputimg, intimg,
						startpos, psf, iterations, maxintensity);

				final double[] pointnearend = peakFitter.GaussianMaskFit.singlegaussianMaskFit(inputimg, intimg, endpos,
						psf, iterations, maxintensity);

				double newslope = (pointnearend[1] - pointnearstart[1]) / (pointnearend[0] - pointnearstart[0]);

				double[] startfit = new double[ndims];
				double[] endfit = new double[ndims];

				final double dx = finalparam[4]/ Math.sqrt( 1 + newslope * newslope);
				final double dy = newslope * finalparam[4];

				final int maxgaussian = 10;

				ArrayList<double[]> startlist = new ArrayList<double[]>();
				ArrayList<double[]> endlist = new ArrayList<double[]>();

				for (int numberofgaussians = 0; numberofgaussians < maxgaussian; numberofgaussians++) {
					startfit = peakFitter.GaussianMaskFit.sumofgaussianMaskFit(inputimg, intimg, startpos, psf,
							iterations, maxintensity, dx, dy, newslope, numberofgaussians, Endfit.Start);

					final double[] startandgauss = { startfit[0], startfit[1], numberofgaussians, newslope };
					startlist.add(startandgauss);
				}
				for (int numberofgaussians = 0; numberofgaussians < maxgaussian; numberofgaussians++) {
					endfit = peakFitter.GaussianMaskFit.sumofgaussianMaskFit(inputimg, intimg, endpos, psf,
							iterations, maxintensity, dx, dy, newslope, numberofgaussians, Endfit.End);
					final double[] endandgauss = { endfit[0], endfit[1], numberofgaussians, newslope };
					endlist.add(endandgauss);

				}

				double mindistance = Double.MAX_VALUE;
				double maxdistance = -Double.MIN_VALUE;
				int finalnumber = 0;
				int finalnumbersecond = 0;
				double startdistance = 0;
				double seconddistance = 0;
				double[] finalstartfit = new double[ndims];
				double[] finalendfit = new double[ndims];
				
				for (int listindex = 0; listindex < startlist.size(); ++listindex) {

					double[] startpoint = { startlist.get(listindex)[0], startlist.get(listindex)[1] };

					
					startdistance = Math.abs(startpoint[1]) + Math.abs(startpoint[0]);

					if (startdistance <= mindistance) {

						mindistance = startdistance;
						finalstartfit = startpoint;
						finalnumber = (int) startlist.get(listindex)[2];

					}

					for (int secondindex = 0; secondindex < endlist.size(); ++secondindex) {
						double[] endpoint = { endlist.get(secondindex)[0], endlist.get(secondindex)[1] };

					
						seconddistance = Math.abs(endpoint[1]) + Math.abs(endpoint[0]);

					//	System.out.println(startpoint[0] + " " + startpoint[1] + " " + endpoint[0] + " " + endpoint[1]);
					

						if (seconddistance >= maxdistance) {

							maxdistance = seconddistance;

							finalendfit = endpoint;

							finalnumbersecond = (int) endlist.get(secondindex)[2];
						}

					}
				}
				System.out.println("Number of gaussians used for start:" + (finalnumber) + " "
						+ "Number of gaussians used for end:" + (finalnumbersecond));
				final double[] refindedparam = { finalstartfit[0], finalstartfit[1], finalendfit[0], finalendfit[1],
						finalparam[4], newslope * finalparam[4], finalparam[5] };

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

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

		final double[] MinandMax = new double[2 * ndims + 2];

		if (slope >= 0) {
			for (int d = 0; d < ndims; ++d) {

				MinandMax[d] = minVal[d];
				MinandMax[d + ndims + 1] = maxVal[d];
			}

		}

		if (slope < 0) {

			MinandMax[0] = minVal[0];
			MinandMax[1] = maxVal[1];
			MinandMax[3] = maxVal[0];
			MinandMax[4] = minVal[1];

		}

		MinandMax[ndims] = 1.5;

		MinandMax[2 * ndims + 1] = maxintensity;
		System.out.println(
				MinandMax[0] + " " + MinandMax[1] + " " + MinandMax[2] + " " + MinandMax[3] + " " + MinandMax[4]);

		for (int d = 0; d < ndims; ++d) {

			if (MinandMax[d] == Double.MAX_VALUE || MinandMax[d + ndims] == Double.MIN_VALUE)
				return null;

		}

		return MinandMax;

	}

	// Get line parameters for fitting line to a line in a label

	public double[] Getfinallineparam(final int label, final double slope, final double intercept, final double[] psf,
			double minlength, double SNR, final Noise noiseorno) throws Exception {

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
			int maxiter = 350;
			double lambda = 1e-3;
			double termepsilon = 1e-1;

			LevenbergMarquardtSolverLine.solve(X, finalparam, fixed_param, I, new GaussianLine(), lambda, termepsilon,
					maxiter);

			// NaN protection: we prefer returning the crude estimate than NaN
			for (int j = 0; j < finalparam.length; j++) {
				if (Double.isNaN(finalparam[j]))
					finalparam[j] = start_param[j];
			}

			final double[] startpos = { finalparam[0], finalparam[1] };
			final double[] endpos = { finalparam[3], finalparam[4] };
			System.out.println(startpos[0] + " " + startpos[1] + " " + endpos[0] + " " + endpos[1] + " "
					+ Distance(startpos, endpos));
			double maxintensity = finalparam[4];
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

				long[] longstartfit = new long[ndims];
				long[] longendfit = new long[ndims];

				final double dx = finalparam[2] / Math.sqrt(1 + slope * slope);
				final double dy = slope * dx / Math.sqrt(1 + slope * slope);

				final double ds = Math.sqrt(dx * dx + dy * dy);
				int maxgaussian;
				switch (noiseorno) {
				case nonoise:
					maxgaussian = (int) Math.ceil((Math.max(psf[0], psf[1])) / ds);
					break;
				case noise:
					maxgaussian = (int) Math.ceil((Math.max(psf[0], psf[1])) / ds);
					break;
				default:	
					maxgaussian = 2;
					break;
				}

				ArrayList<double[]> startlist = new ArrayList<double[]>();
				ArrayList<double[]> endlist = new ArrayList<double[]>();
				final Noiseclassifier noiseinterval = new Noiseclassifier(inputimg, intimg);
				final long radius = (long) Math.ceil(Math.sqrt(psf[0] * psf[0] + psf[1] * psf[1]));
				final int numberofgaussians = maxgaussian;
				// for (int numberofgaussians = 0; numberofgaussians <
				// maxgaussian; numberofgaussians++) {
				startfit = peakFitter.GaussianMaskFit.sumofgaussianMaskFit(inputimg, intimg, startpos, psf, iterations,
						maxintensity, dx, dy, slope, numberofgaussians, Endfit.Start, label);

				if (startfit != null) {
					switch (noiseorno) {
					case noise:
						RandomAccess<FloatType> ranac = inputimg.randomAccess();
						for (int d = 0; d < ndims; ++d) {
							longstartfit[d] = (long) startfit[d];

						}
						ranac.setPosition(longstartfit);
						double[] noiseparam = noiseinterval.Getnoiseparam(ranac, radius);
						if (noiseparam[0] > 1.0 / SNR) {
							final double[] startandgauss = { startfit[0], startfit[1], numberofgaussians, slope };
							startlist.add(startandgauss);
						}
						
						break;
					case nonoise:
						final double[] startandgauss = { startfit[0], startfit[1], numberofgaussians, slope };
						startlist.add(startandgauss);
						break;
					}
				
				}
				

						
					
				
				// }
				// for (int numberofgaussians = 0; numberofgaussians <
				// maxgaussian; numberofgaussians++) {
				endfit = peakFitter.GaussianMaskFit.sumofgaussianMaskFit(inputimg, intimg, endpos, psf, iterations,
						maxintensity, dx, dy, slope, numberofgaussians, Endfit.End, label);
				if (endfit != null) {
					switch(noiseorno){
						
					case noise:
						
						RandomAccess<FloatType> ranac = inputimg.randomAccess();
						for (int d = 0; d < ndims; ++d) {
							longendfit[d] = (long) endfit[d];
						}
						ranac.setPosition(longendfit);
						double[] noiseparam = noiseinterval.Getnoiseparam(ranac, radius);

						if (noiseparam[0] > 1.0 / SNR) {

							final double[] endandgauss = { endfit[0], endfit[1], numberofgaussians, slope };
							endlist.add(endandgauss);
						}
						
						break;
						
					case nonoise:
						final double[] endandgauss = { endfit[0], endfit[1], numberofgaussians, slope };
						endlist.add(endandgauss);

						break;
					}

					

				}
				// }

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

					System.out.println(startpoint[0] + " " + startpoint[1]);
				}
				for (int secondindex = 0; secondindex < endlist.size(); ++secondindex) {
					double[] endpoint = { endlist.get(secondindex)[0], endlist.get(secondindex)[1] };

					seconddistance = Math.abs(endpoint[1]) + Math.abs(endpoint[0]);

					System.out.println(endpoint[0] + " " + endpoint[1]);

					if (seconddistance >= maxdistance) {

						maxdistance = seconddistance;

						finalendfit = endpoint;

						finalnumbersecond = (int) endlist.get(secondindex)[2];
					}

				}

				System.out.println("Number of gaussians used for start:" + (finalnumber) + " "
						+ "Number of gaussians used for end:" + (finalnumbersecond));
				final double[] refindedparam = { finalstartfit[0], finalstartfit[1], finalendfit[0], finalendfit[1], dx,
						dy, finalparam[5] };

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

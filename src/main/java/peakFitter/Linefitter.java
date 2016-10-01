package peakFitter;

import java.io.FileWriter;
import java.io.IOException;
import houghandWatershed.PerformWatershedding;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.Point;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import peakFitter.GaussianMaskFit.Endfit;
import preProcessing.GetLocalmaxmin;

public class Linefitter {

	private final RandomAccessibleInterval<FloatType> inputimg;
	private final RandomAccessibleInterval<IntType> intimg;
	private final int ndims;
	// LM solver iteration params
	final int maxiter = 500;
	final double lambda = 1e-3;
	final double termepsilon = 1e-1;
	//Mask fits iteration param
	final int iterations = 500;
	final double cutoffdistance = 6;

	public Linefitter(RandomAccessibleInterval<FloatType> inputimg, RandomAccessibleInterval<IntType> intimg) {

		this.inputimg = inputimg;
		this.intimg = intimg;
		this.ndims = inputimg.numDimensions();
		assert inputimg.numDimensions() == intimg.numDimensions();

	}

	private final double[] MakeimprovedLineguess(double slope, double intercept, double[] psf, int label,
			boolean offsetting) throws Exception {
		long[] newposition = new long[ndims];
		double[] minVal = { Double.MAX_VALUE, Double.MAX_VALUE };
		double[] maxVal = { -Double.MIN_VALUE, -Double.MIN_VALUE };

		RandomAccessibleInterval<FloatType> currentimg = PerformWatershedding.CurrentLabelImage(intimg, inputimg,
				label);

		double newintercept = intercept;

		final long[] minCorner = PerformWatershedding.GetMincorners(intimg, label);
		final long[] maxCorner = PerformWatershedding.GetMaxcorners(intimg, label);
		FinalInterval smallinterval = new FinalInterval(minCorner, maxCorner);
		currentimg = Views.interval(currentimg, smallinterval);
		if (offsetting) {

			currentimg = Views.offsetInterval(currentimg, smallinterval);

			newintercept = intercept - (smallinterval.realMin(1) - slope * smallinterval.realMin(0));

		}

		final Cursor<FloatType> outcursor = Views.iterable(currentimg).localizingCursor();

		final double maxintensityline = GetLocalmaxmin.computeMaxIntensity(currentimg);

		while (outcursor.hasNext()) {

			outcursor.fwd();

			if (outcursor.get().get() / maxintensityline > 0.5) {
				outcursor.localize(newposition);

				long pointonline = (long) (newposition[1] - slope * newposition[0] - newintercept);

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

		final double[] MinandMax = new double[2 * ndims + 3];

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

		// This parameter is guess estimate for spacing between the Gaussians
		MinandMax[2 * ndims] =   Math.min(psf[0], psf[1]) ;
		MinandMax[2 * ndims + 1] = maxintensityline; 
		// This parameter guess estimates the background noise level
		MinandMax[2 * ndims + 2] = 1.0; 
		
		if (offsetting == false)
		System.out.println("Label: " + label + " " + "Hough Detection: " + " StartX: " + MinandMax[0] + " StartY: "
				+ MinandMax[1] + " EndX: " + MinandMax[2] + " EndY: " + MinandMax[3]);

		if (offsetting) {
			System.out.println("Ofsetting on: " + "Label: " + label + " " + "Hough Detection: " + " StartX: "
					+ (MinandMax[0] + smallinterval.realMin(0)) + " StartY: "
					+ (MinandMax[1] + smallinterval.realMin(1)) + " EndX: " + (MinandMax[2] + smallinterval.realMin(0))
					+ " EndY: " + (MinandMax[3] + smallinterval.realMin(1)));

		}
		if (offsetting == false) {
			for (int d = 0; d < ndims; ++d) {

				if (MinandMax[d] == Double.MAX_VALUE || MinandMax[d + ndims] == -Double.MIN_VALUE)
					return null;
				if (MinandMax[d] >= inputimg.dimension(d) || MinandMax[d + ndims] >= inputimg.dimension(d))
					return null;
				if (MinandMax[d] <= 0 || MinandMax[d + ndims] <= 0)
					return null;

			}
		}

		if (offsetting) {

			for (int d = 0; d < ndims; ++d) {

				if (MinandMax[d] + smallinterval.realMin(d) == Double.MAX_VALUE
						|| MinandMax[d + ndims] + smallinterval.realMin(d) == -Double.MIN_VALUE)
					return null;
				if (MinandMax[d] + smallinterval.realMin(d) >= inputimg.dimension(d)
						|| MinandMax[d + ndims] + smallinterval.realMin(d) >= inputimg.dimension(d))
					return null;
				if (MinandMax[d] + smallinterval.realMin(d) <= 0
						|| MinandMax[d + ndims] + smallinterval.realMin(d) <= 0)
					return null;

			}

		}

		return MinandMax;

	}

	private final double[] MakerepeatedLineguess(double[] iniparam, int label, boolean offsetting) throws Exception {
		long[] newposition = new long[ndims];
		double[] minVal = { Double.MAX_VALUE, Double.MAX_VALUE };
		double[] maxVal = { -Double.MIN_VALUE, -Double.MIN_VALUE };

		RandomAccessibleInterval<FloatType> currentimg = PerformWatershedding.CurrentLabelImage(intimg, inputimg,
				label);
		final double[] cordone = { iniparam[0], iniparam[1] };
		final double[] cordtwo = { iniparam[2], iniparam[3] };

		double slope = (cordone[1] - cordtwo[1]) / (cordone[0] - cordtwo[0]);
		double intercept = cordone[1] - slope * cordone[0];
		double newintercept = intercept;

		final long[] minCorner = PerformWatershedding.GetMincorners(intimg, label);
		final long[] maxCorner = PerformWatershedding.GetMaxcorners(intimg, label);
		FinalInterval smallinterval = new FinalInterval(minCorner, maxCorner);
		currentimg = Views.interval(currentimg, smallinterval);
		if (offsetting) {

			currentimg = Views.offsetInterval(currentimg, smallinterval);

			newintercept = intercept - (smallinterval.realMin(1) - slope * smallinterval.realMin(0));

		}

		final Cursor<FloatType> outcursor = Views.iterable(currentimg).localizingCursor();

		final double maxintensityline = GetLocalmaxmin.computeMaxIntensity(currentimg);

		while (outcursor.hasNext()) {

			outcursor.fwd();

			if (outcursor.get().get() / maxintensityline > 0.5) {
				outcursor.localize(newposition);

				long pointonline = (long) (outcursor.getLongPosition(1) - slope * outcursor.getLongPosition(0) - newintercept);

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

		final double[] MinandMax = new double[2 * ndims + 3];

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

		MinandMax[2 * ndims] = iniparam[2 * ndims];
		MinandMax[2 * ndims + 1] = iniparam[2 * ndims + 1];
		MinandMax[2 * ndims + 2] = iniparam[2 * ndims + 2];
		System.out.println("Label: " + label + " " + "Initial guess: " + " StartX: " + MinandMax[0] + " StartY: "
				+ MinandMax[1] + " EndX: " + MinandMax[2] + " EndY: " + MinandMax[3]);

		if (offsetting) {
			System.out.println("Ofsetting on: " + "Label: " + label + " " + "Initial guess: " + " StartX: "
					+ (MinandMax[0] + smallinterval.realMin(0)) + " StartY: "
					+ (MinandMax[1] + smallinterval.realMin(1)) + " EndX: " + (MinandMax[2] + smallinterval.realMin(0))
					+ " EndY: " + (MinandMax[3] + smallinterval.realMin(1)));

		}
		if (offsetting == false) {
			for (int d = 0; d < ndims; ++d) {

				if (MinandMax[d] == Double.MAX_VALUE || MinandMax[d + ndims] == -Double.MIN_VALUE)
					return null;
				if (MinandMax[d] >= inputimg.dimension(d) || MinandMax[d + ndims] >= inputimg.dimension(d))
					return null;
				if (MinandMax[d] <= 0 || MinandMax[d + ndims] <= 0)
					return null;

			}
		}

		if (offsetting) {

			for (int d = 0; d < ndims; ++d) {

				if (MinandMax[d] + smallinterval.realMin(d) == Double.MAX_VALUE
						|| MinandMax[d + ndims] + smallinterval.realMin(d) == -Double.MIN_VALUE)
					return null;
				if (MinandMax[d] + smallinterval.realMin(d) >= inputimg.dimension(d)
						|| MinandMax[d + ndims] + smallinterval.realMin(d) >= inputimg.dimension(d))
					return null;
				if (MinandMax[d] + smallinterval.realMin(d) <= 0
						|| MinandMax[d + ndims] + smallinterval.realMin(d) <= 0)
					return null;

			}

		}

		return MinandMax;

	}

	// Get line parameters for fitting line to a line in a label

	public double[] Getfinallineparam(final int label, final double slope, final double intercept, final double[] psf,
			double minlength, final boolean offsetting) throws Exception {

		PointSampleList<FloatType> datalist = gatherfullData(label, offsetting);
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

		final double[] start_param = MakeimprovedLineguess(slope, intercept, psf, label, offsetting);
		if (start_param == null)
			return null;

		else {

			final double[] finalparamstart = start_param.clone();
			// LM solver part

			RandomAccessibleInterval<FloatType> currentimg = PerformWatershedding.CurrentLabelImage(intimg, inputimg,
					label);

			final long[] minCorner = PerformWatershedding.GetMincorners(intimg, label);
			final long[] maxCorner = PerformWatershedding.GetMaxcorners(intimg, label);
			FinalInterval smallinterval = new FinalInterval(minCorner, maxCorner);
			currentimg = Views.interval(currentimg, smallinterval);
			if (offsetting) {

				currentimg = Views.offsetInterval(currentimg, smallinterval);

			}

			final double[] fixed_param = new double[ndims + 1];

			for (int d = 0; d < ndims; ++d) {

				fixed_param[d] = 1.0 / Math.pow(psf[d], 2);
			}

			fixed_param[ndims] = GetLocalmaxmin.computeMaxIntensity(currentimg);
		
			double inicutoffdistance = 0;
			if (offsetting == false){
			final double[] inistartpos = { start_param[0], start_param[1] };
			final double[] iniendpos = { start_param[2], start_param[3] };
			inicutoffdistance = Distance(inistartpos, iniendpos);
			
			
			}
			
			if (offsetting){
				final double[] inistartpos = { start_param[0] + smallinterval.realMin(0), start_param[1] + smallinterval.realMin(1)};
			    final double[] iniendpos = { start_param[2]+ smallinterval.realMin(0), start_param[3] + smallinterval.realMin(1)};
			    inicutoffdistance = Distance(inistartpos, iniendpos);
			}
			
			if (inicutoffdistance > minlength) {
				LevenbergMarquardtSolverLine.solve(X, finalparamstart, fixed_param, I, new GaussianLineds(), lambda,
						termepsilon, maxiter);

				final double[] startpos = { finalparamstart[0], finalparamstart[1] };
				final double[] endpos = { finalparamstart[2], finalparamstart[3] };
				// NaN protection: we prefer returning the crude estimate than
				// NaN
				for (int j = 0; j < finalparamstart.length; j++) {
					if (Double.isNaN(finalparamstart[j]))
						finalparamstart[j] = start_param[j];
				}


				double newslope = (endpos[1] - startpos[1]) / (endpos[0] - startpos[0]);
				double newintercept = (endpos[1] - newslope * endpos[0]);
				double dx = finalparamstart[4]/ Math.sqrt(1 + newslope * newslope);
				double dy = newslope * dx;
				double ds = finalparamstart[4];
				final double LMdist = sqDistance(startpos, endpos);

				double[] dxvector = { dx, dy };

				double[] startfit = new double[ndims];
				double[] endfit = new double[ndims];
				final double maxintensityline = fixed_param[ndims];

				final int numberofgaussians = 2;

				startfit = peakFitter.GaussianMaskFit.sumofgaussianMaskFit(currentimg, intimg, startpos.clone(), psf,
						iterations, dxvector, newslope, newintercept, maxintensityline, numberofgaussians, Endfit.Start,
						label);

				endfit = peakFitter.GaussianMaskFit.sumofgaussianMaskFit(currentimg, intimg, endpos.clone(), psf,
						iterations, dxvector, newslope, newintercept, maxintensityline, numberofgaussians, Endfit.End,
						label);

				final double Maskdist = sqDistance(startfit, endfit);
				// If mask fits fail, return LM solver results, very crucial for
				// noisy data

				double[] returnparam = new double[2 * ndims + 3];
				double[] LMparam = new double[2 * ndims];

				if (offsetting) {
					for (int d = 0; d < ndims; ++d) {
						startpos[d] += smallinterval.realMin(d);
						endpos[d] += smallinterval.realMin(d);
						startfit[d] += smallinterval.realMin(d);
						endfit[d] += smallinterval.realMin(d);
					}

				}

				for (int d = 0; d < ndims; ++d) {
					LMparam[d] = startpos[d];
					LMparam[ndims + d] = endpos[d];
				}

				for (int d = 0; d < ndims; ++d) {
					returnparam[d] = startfit[d];
					returnparam[ndims + d] = endfit[d];
				}
				
				
				if (Math.abs(Math.sqrt(Maskdist)) - Math.sqrt(LMdist) >= cutoffdistance){
				if (Math.abs(startpos[0] - startfit[0]) >= cutoffdistance / 2 && Math.abs(startpos[1] - startfit[1]) >= cutoffdistance / 2
						|| Math.abs(endpos[0] - endfit[0]) >= cutoffdistance / 2 && Math.abs(endpos[1] - endfit[1]) >= cutoffdistance / 2 ){
					System.out.println("Mask fits fail, returning LM solver results!");

					for (int d = 0; d < ndims; ++d) {
						returnparam[d] = startpos[d];
						returnparam[ndims + d] = endpos[d];
					}
					}
				if (Math.abs(startpos[0] - startfit[0]) >= cutoffdistance || Math.abs(startpos[1] - startfit[1]) >= cutoffdistance 
						|| Math.abs(endpos[0] - endfit[0]) >= cutoffdistance  || Math.abs(endpos[1] - endfit[1]) >= cutoffdistance  ){
					System.out.println("Mask fits fail, returning LM solver results!");
					for (int d = 0; d < ndims; ++d) {
						returnparam[d] = startpos[d];
						returnparam[ndims + d] = endpos[d];
					}
					
				}
				
				
				
				
				}

				for (int d = 0; d < ndims; ++d) {
					if (Double.isNaN(startfit[d]) || Double.isNaN(endfit[d])) {
						System.out.println("Mask fits fail, returning LM solver results!");
						returnparam[d] = startpos[d];
						returnparam[ndims + d] = endpos[d];

					}
				}
				Testerrorthree(returnparam, label, Distance(new double[] { returnparam[0], returnparam[1] },
						new double[] { returnparam[2], returnparam[3] }));

				returnparam[2 * ndims] = finalparamstart[4];
				returnparam[2 * ndims + 1] = finalparamstart[5];
				returnparam[2 * ndims + 2] = finalparamstart[6];
				return returnparam;

			}

			else
				return null;

		}
	}

	public double[] Getfinaltrackparam(final double[] iniparam, final int label, final double[] psf, final int rate,
			boolean offsetting) throws Exception {

		if (iniparam == null)
			return null;

		else {

			PointSampleList<FloatType> datalist = gatherfullData(label, offsetting);
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

			final double[] finalparamstart = MakerepeatedLineguess(iniparam, label, offsetting);
			if (finalparamstart == null)
				return null;

			else {
				RandomAccessibleInterval<FloatType> currentimg = PerformWatershedding.CurrentLabelImage(intimg,
						inputimg, label);

				final long[] minCorner = PerformWatershedding.GetMincorners(intimg, label);
				final long[] maxCorner = PerformWatershedding.GetMaxcorners(intimg, label);
				FinalInterval smallinterval = new FinalInterval(minCorner, maxCorner);
				currentimg = Views.interval(currentimg, smallinterval);
				if (offsetting) {

					currentimg = Views.offsetInterval(currentimg, smallinterval);

				}

				final double[] fixed_param = new double[ndims + 1];

				for (int d = 0; d < ndims; ++d) {

					fixed_param[d] = 1.0 / Math.pow(psf[d], 2);
				}
				fixed_param[ndims] = GetLocalmaxmin.computeMaxIntensity(currentimg);

				final double[] inistartpos = { finalparamstart[0], finalparamstart[1] };
				final double[] iniendpos = { finalparamstart[2], finalparamstart[3] };

				double inicutoffdistance = Distance(inistartpos, iniendpos);

				// LM solver part
				if (inicutoffdistance > 1) {
					LevenbergMarquardtSolverLine.solve(X, finalparamstart, fixed_param, I, new GaussianLineds(), lambda,
							termepsilon, maxiter);

					final double[] startpos = { finalparamstart[0], finalparamstart[1] };
					final double[] endpos = { finalparamstart[2], finalparamstart[3] };
					// NaN protection: we prefer returning the crude estimate
					// than
					// NaN
					for (int j = 0; j < finalparamstart.length; j++) {
						if (Double.isNaN(finalparamstart[j]))
							finalparamstart[j] = iniparam[j];
					}

					final double LMdist = sqDistance(startpos, endpos);

					double[] returnparam = new double[2 * ndims + 3];

					final double maxintensityline = fixed_param[ndims];

					

					double newslope = (endpos[1] - startpos[1]) / (endpos[0] - startpos[0]);
					double newintercept = (endpos[1] - newslope * endpos[0]);
					double dx = finalparamstart[4] / Math.sqrt(1 + newslope * newslope);
					double dy = newslope * dx;
					double ds = finalparamstart[4];
					double[] dxvector = { dx,  dy };

					double[] startfit = new double[ndims];
					double[] endfit = new double[ndims];

					final int numberofgaussians = 2;

					startfit = peakFitter.GaussianMaskFit.sumofgaussianMaskFit(currentimg, intimg, startpos.clone(),
							psf, iterations, dxvector, newslope, newintercept, maxintensityline, numberofgaussians,
							Endfit.Start, label);

					endfit = peakFitter.GaussianMaskFit.sumofgaussianMaskFit(currentimg, intimg, endpos.clone(), psf,
							iterations, dxvector, newslope, newintercept, maxintensityline, numberofgaussians,
							Endfit.End, label);

					final double Maskdist = sqDistance(startfit, endfit);
					// If mask fits fail, return LM solver results, very crucial
					// for
					// noisy data

					if (offsetting) {
						for (int d = 0; d < ndims; ++d) {
							startpos[d] += smallinterval.realMin(d);
							endpos[d] += smallinterval.realMin(d);
							startfit[d] += smallinterval.realMin(d);
							endfit[d] += smallinterval.realMin(d);
						}

					}

					for (int d = 0; d < ndims; ++d) {
						returnparam[d] = startfit[d];
						returnparam[ndims + d] = endfit[d];
					}

					
					if (Math.abs(Math.sqrt(Maskdist)) - Math.sqrt(LMdist) > cutoffdistance){
						if (Math.abs(startpos[0] - startfit[0]) >= cutoffdistance / 2 && Math.abs(startpos[1] - startfit[1]) >= cutoffdistance / 2
								|| Math.abs(endpos[0] - endfit[0]) >= cutoffdistance / 2 && Math.abs(endpos[1] - endfit[1]) >= cutoffdistance / 2 ){
							System.out.println("Mask fits fail, returning LM solver results!");
						
							for (int d = 0; d < ndims; ++d) {
							returnparam[d] = startpos[d];
							returnparam[ndims + d] = endpos[d];
						}
						}
					
						if (Math.abs(startpos[0] - startfit[0]) >= cutoffdistance || Math.abs(startpos[1] - startfit[1]) >= cutoffdistance 
								|| Math.abs(endpos[0] - endfit[0]) >= cutoffdistance  || Math.abs(endpos[1] - endfit[1]) >= cutoffdistance  ){
							System.out.println("Mask fits fail, returning LM solver results!");
							for (int d = 0; d < ndims; ++d) {
								returnparam[d] = startpos[d];
								returnparam[ndims + d] = endpos[d];
							}
							
						}
					
					
					}
					
					
					
					for (int d = 0; d < ndims; ++d) {
						if (Double.isNaN(startfit[d]) || Double.isNaN(endfit[d])) {
							System.out.println("Mask fits fail, returning LM solver results!");
							returnparam[d] = startpos[d];
							returnparam[ndims + d] = endpos[d];

						}

					}

				//	Testmovingline(returnparam, label, Distance(new double[] { returnparam[0], returnparam[1] },
				//			new double[] { returnparam[2], returnparam[3] }), rate);

					returnparam[2 * ndims] = finalparamstart[4];
					returnparam[2 * ndims + 1] = finalparamstart[5];
					returnparam[2 * ndims + 2] = finalparamstart[6];
					return returnparam;
				} else
					return null;
			}

		}

	}

	public int Getlabel(final Point linepoint) {

		RandomAccess<IntType> intranac = intimg.randomAccess();

		intranac.setPosition(linepoint);
		int currentlabel = intranac.get().get();

		return currentlabel;
	}

	private void Testerrorone(final double[] point, int label, double length) {

		// Errorlist 
		try {
			FileWriter writer = new FileWriter("../res/error-Pnoise1snr35.txt", true);

			if (label == 1)
				writer.write((point[0] - 55.55217459454749) + " " + (point[1] - 48.30096130273662) + " "
						+ (point[2] - 66.32106541633595) + " " + (point[3] - 20.62906649269944) + " "
						+ (length - 29.69348029297614));
			writer.write("\r\n");
			if (label == 2)
				writer.write((point[0] - 213.71705870001662) + " " + (point[1] - 22.18575371852285) + " "
						+ (point[2] - 224.17901898271776) + " " + (point[3] - 49.02869036045246) + " "
						+ (length - 28.809648739952788));
			writer.write("\r\n");
			if (label == 3)
				writer.write((point[0] - 382.40921141961365) + " " + (point[1] - 25.47946702808274) + " "
						+ (point[2] - 385.43282134647404) + " " + (point[3] - 39.50964398952846) + " "
						+ (length - 14.352284924683282));
			writer.write("\r\n");
			if (label == 4)
				writer.write((point[0] - 102.44291073766455) + " " + (point[1] - 33.58242883236636) + " "
						+ (point[2] - 129.6755712883426) + " " + (point[3] - 45.01986810740948) + " "
						+ (length - 29.53697374205347));
			writer.write("\r\n");
			
			if (label == 6)
				writer.write((point[0] - 432.8011066191138) + " " + (point[1] - 104.44891402543196 ) + " "
						+ (point[2] - 448.33893858877803) + " " + (point[3] - 82.55405768010598) + " "
						+ (length - 26.847885516367604));
			writer.write("\r\n");

			if (label == 7)
				writer.write((point[0] - 183.32000789300656) + " " + (point[1] - 158.22256763937276) + " "
						+ (point[2] - 189.96789844665068) + " " + (point[3] - 149.1111230957096) + " "
						+ (length - 11.278868315814314));
			writer.write("\r\n");
			if (label == 8)
				writer.write((point[0] - 446.1653019733676) + " " + (point[1] - 252.3134572376292) + " "
						+ (point[2] - 449.2993340411474) + " " + (point[3] - 227.84043364472504) + " "
						+ (length - 24.67288067455269));
			writer.write("\r\n");
			
			if (label == 9)
				writer.write((point[0] - 119.80528923622467) + " " + (point[1] - 248.9289843415045) + " "
						+ (point[2] - 129.17793025086274) + " " + (point[3] - 265.8217915876105) + " "
						+ (length - 19.318730192312522));
			writer.write("\r\n");
			if (label == 10)
				writer.write((point[0] - 276.8721035260036) + " " + (point[1] - 273.3940355635414) + " "
						+ (point[2] - 287.7085268575912) + " " + (point[3] - 274.71972909508486) + " "
						+ (length - 10.917212737734472));
			writer.write("\r\n");
			if (label == 11)
				writer.write((point[0] - 13.700439658553767) + " " + (point[1] - 310.5532653119629) + " "
						+ (point[2] - 31.153396540922103) + " " + (point[3] - 286.6858019889502) + " "
						+ (length - 29.56791351132449));
			writer.write("\r\n");
			if (label == 12)
				writer.write((point[0] - 309.35892754128355 ) + " " + (point[1] - 358.1480590965386) + " "
						+ (point[2] - 319.39022531594753) + " " + (point[3] - 371.67661027458706 ) + " "
						+ (length - 16.841871393080265));
			writer.write("\r\n");
			

			if (label == 14)
				writer.write((point[0] - 111.18093055644358) + " " + (point[1] - 376.8726796271507) + " "
						+ (point[2] - 122.24121775748232) + " " + (point[3] - 391.0255834923995) + " "
						+ (length - 17.962033314422822));
			writer.write("\r\n");

			if (label == 15)
				writer.write((point[0] -328.2849742729377) + " " + (point[1] - 404.7499831894944) + " "
						+ (point[2] - 345.6947463117291) + " " + (point[3] - 387.60846978547454) + " "
						+ (length - 24.432184597838884));
			writer.write("\r\n");

			if (label == 16)
				writer.write((point[0] - 14.71094251606357) + " " + (point[1] - 418.7655285832544) + " "
						+ (point[2] - 35.60154650513873) + " " + (point[3] - 416.25106429920544) + " "
						+ (length - 21.041384594748536));
			writer.write("\r\n");
			if (label == 17)
				writer.write((point[0] - 333.67781551268257) + " " + (point[1] - 426.36288843507003) + " "
						+ (point[2] - 359.2263637922085) + " " + (point[3] - 417.7641149277784) + " "
						+ (length - 26.95676584868753));
			writer.write("\r\n");
			if (label == 18)
				writer.write((point[0] - 384.84869591121316) + " " + (point[1] - 443.35195442039355) + " "
						+ (point[2] - 400.9015584864585) + " " + (point[3] - 431.87031173913215) + " "
						+ (length - 19.736324772355065));
			writer.write("\r\n");
			if (label == 19)
				writer.write((point[0] - 223.07334445046993) + " " + (point[1] - 442.6587430230545) + " "
						+ (point[2] - 250.90911483273447) + " " + (point[3] - 441.94146302339226) + " "
						+ (length - 27.845010385562283));
			writer.write("\r\n");
			
			

			writer.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void Testerrorthree(final double[] point, int label, double length) {

		// Errorlist for Fake_big_file.tif
		try {
			FileWriter writer = new FileWriter("../res/error-Pnonoise3.txt", true);

			if (label == 1)
				writer.write((point[0] - 115.47053635052531) + " " + (point[1] - 9.023194352830702) + " "
						+ (point[2] - 122.54159964459754 ) + " " + (point[3] - 17.518385880585186) + " "
						+ (length - 11.05297313857259));
			writer.write("\r\n");
			if (label == 2)
				writer.write((point[0] - 357.6402383747052) + " " + (point[1] - 12.579259572881716) + " "
						+ (point[2] - 373.54103649310053) + " " + (point[3] - 12.561613441134531 ) + " "
						+ (length - 15.900807909912343));
			writer.write("\r\n");
			if (label == 3)
				writer.write((point[0] - 422.5586733459164) + " " + (point[1] - 16.290977908274005) + " "
						+ (point[2] - 435.4336281855708) + " " + (point[3] - 26.26657858066185) + " "
						+ (length - 16.287328537795403));
			writer.write("\r\n");
			if (label == 4)
				writer.write((point[0] - 117.57934569226755) + " " + (point[1] - 55.570567093715496) + " "
						+ (point[2] - 140.5646361981157) + " " + (point[3] - 39.696394685092) + " "
						+ (length - 27.93408185885006));
			writer.write("\r\n");
			if (label == 5)
				writer.write((point[0] - 58.987715024176104) + " " + (point[1] - 57.55851541154964) + " "
						+ (point[2] - 74.94713093352168) + " " + (point[3] - 49.67797967614733) + " "
						+ (length - 17.799039289928697));
			writer.write("\r\n");
			if (label == 6)
				writer.write((point[0] - 328.44586313581783) + " " + (point[1] - 100.24824829548898) + " "
						+ (point[2] - 338.0883235637351 ) + " " + (point[3] - 77.82862502979053) + " "
						+ (length - 24.405256615733382));
			writer.write("\r\n");

			
			if (label == 8)
				writer.write((point[0] - 168.95040819988822) + " " + (point[1] - 132.82573384052787) + " "
						+ (point[2] - 178.95745390755877) + " " + (point[3] - 128.7433252784837) + " "
						+ (length - 10.807729801529058));
			writer.write("\r\n");
			
			if (label == 11)
				writer.write((point[0] - 272.2561724772123) + " " + (point[1] - 267.4288436625433) + " "
						+ (point[2] - 284.8912784517835) + " " + (point[3] - 242.9391645537842) + " "
						+ (length - 27.55703695680351));
			writer.write("\r\n");
			if (label == 12)
				writer.write((point[0] - 343.6518316899644) + " " + (point[1] - 260.0891884124472) + " "
						+ (point[2] - 370.71778225694993) + " " + (point[3] - 258.9261650453071) + " "
						+ (length - 27.09092658893412));
			writer.write("\r\n");
			if (label == 14)
				writer.write((point[0] - 77.39669991423358) + " " + (point[1] - 354.2685283707751) + " "
						+ (point[2] - 82.76833457785416) + " " + (point[3] - 374.0734877135922) + " "
						+ (length - 20.52049886162737));
			writer.write("\r\n");

			if (label == 15)
				writer.write((point[0] - 243.19267085703873) + " " + (point[1] - 366.9101410699944) + " "
						+ (point[2] - 250.7120614406193) + " " + (point[3] - 375.48033078965346) + " "
						+ (length - 11.401288812208477));
			writer.write("\r\n");

			

			if (label == 16)
				writer.write((point[0] - 96.81917869639119 ) + " " + (point[1] - 384.6837984735284) + " "
						+ (point[2] - 115.26948850758293) + " " + (point[3] - 381.69979905505596 ) + " "
						+ (length - 18.690055769269446));
			writer.write("\r\n");
			if (label == 17)
				writer.write((point[0] - 49.3170184804688) + " " + (point[1] - 425.767490948651) + " "
						+ (point[2] - 63.63103519112512) + " " + (point[3] - 413.79617785393265) + " "
						+ (length - 18.660209312993352));
			writer.write("\r\n");
			if (label == 18)
				writer.write((point[0] - 211.52755886980339) + " " + (point[1] - 433.8468820584074) + " "
						+ (point[2] - 229.50961952188607) + " " + (point[3] - 434.6956128998127 ) + " "
						+ (length - 18.002079028165962));
			writer.write("\r\n");
			if (label == 19)
				writer.write((point[0] - 295.2103675652863) + " " + (point[1] - 433.9056000528416) + " "
						+ (point[2] - 323.01968348294497) + " " + (point[3] - 443.9300793954226) + " "
						+ (length - 29.560924171919567));
			writer.write("\r\n");
			if (label == 20)
				writer.write((point[0] - 86.8988927325537) + " " + (point[1] - 445.970800094747) + " "
						+ (point[2] - 106.58885746641593) + " " + (point[3] - 467.5851100482377) + " "
						+ (length - 29.238213112061167));
			writer.write("\r\n");

			writer.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}


	private void Testerrortwo(final double[] point, int label, double length) {

		// Errorlist for Fake_big_file.tif
		try {
			FileWriter writer = new FileWriter("../res/error-Pnonoise2.txt", true);

			if (label == 1)
				writer.write((point[0] - 362.19642734312356) + " " + (point[1] - 13.388136845729605) + " "
						+ (point[2] - 391.6473283989393) + " " + (point[3] - 15.641411094308705) + " "
						+ (length - 29.5369737420535));
			writer.write("\r\n");
			if (label == 2)
				writer.write((point[0] - 328.24475874783906) + " " + (point[1] - 31.875423603617193) + " "
						+ (point[2] - 346.7232577746829) + " " + (point[3] - 15.891784533061049 ) + " "
						+ (length - 24.432184597838933));
			writer.write("\r\n");
			if (label == 3)
				writer.write((point[0] - 246.54081005002695) + " " + (point[1] - 42.00209390745403) + " "
						+ (point[2] - 270.9480464426979) + " " + (point[3] - 38.39129725671395) + " "
						+ (length - 24.672880674552708));
			writer.write("\r\n");
			if (label == 4)
				writer.write((point[0] - 188.3005499067052) + " " + (point[1] - 56.30346049666297) + " "
						+ (point[2] - 217.0825161243029) + " " + (point[3] - 49.00269392332098) + " "
						+ (length - 29.693480292976147));
			writer.write("\r\n");
			
			if (label == 6)
				writer.write((point[0] - 438.9030794375588) + " " + (point[1] - 71.13858404163925) + " "
						+ (point[2] - 441.8158510394937) + " " + (point[3] - 82.03485092899298) + " "
						+ (length - 11.278868315814298));
			writer.write("\r\n");

			if (label == 7)
				writer.write((point[0] - 256.680078177105) + " " + (point[1] - 114.09772454324029) + " "
						+ (point[2] - 266.8128414214157) + " " + (point[3] - 110.03414939373252) + " "
						+ (length - 10.917212737734483));
			writer.write("\r\n");
			if (label == 8)
				writer.write((point[0] - 149.86939978852126) + " " + (point[1] - 116.74320280616583) + " "
						+ (point[2] - 164.7378527951561) + " " + (point[3] - 124.65381196163204) + " "
						+ (length - 16.84187139308024));
			writer.write("\r\n");
			if (label == 9)
				writer.write((point[0] - 199.08227189244099) + " " + (point[1] - 185.3749862504838) + " "
						+ (point[2] - 213.2208981598052) + " " + (point[3] - 210.47668156755876) + " "
						+ (length - 28.8096487399528));
			writer.write("\r\n");
			if (label == 10)
				writer.write((point[0] - 75.95400961802716) + " " + (point[1] - 202.68483041367028) + " "
						+ (point[2] - 91.89085485376414) + " " + (point[3] - 191.76556761000376) + " "
						+ (length - 19.318730192312536));
			writer.write("\r\n");
			if (label == 11)
				writer.write((point[0] - 335.7895014938657) + " " + (point[1] - 210.3897819816849) + " "
						+ (point[2] - 349.1797000191111) + " " + (point[3] - 194.15887168316385) + " "
						+ (length - 21.04138459474854));
			writer.write("\r\n");
			if (label == 12)
				writer.write((point[0] - 127.07526194900117) + " " + (point[1] - 224.99882903292337) + " "
						+ (point[2] - 133.11537907728987) + " " + (point[3] - 233.91590943312383) + " "
						+ (length - 10.770206023428045));
			writer.write("\r\n");
			if (label == 14)
				writer.write((point[0] - 448.6166745061395) + " " + (point[1] - 242.45669201167456) + " "
						+ (point[2] - 451.3954104349838) + " " + (point[3] - 256.5374132004859) + " "
						+ (length - 14.352284924683275));
			writer.write("\r\n");

			if (label == 15)
				writer.write((point[0] - 102.07326282107839) + " " + (point[1] - 296.48350789270467) + " "
						+ (point[2] - 113.23772070041582) + " " + (point[3] - 280.2084534909407) + " "
						+ (length - 19.736324772355058));
			writer.write("\r\n");

			

			if (label == 16)
				writer.write((point[0] - 30.136884216334302 ) + " " + (point[1] - 344.90248585105746) + " "
						+ (point[2] - 55.51175123493876) + " " + (point[3] - 356.36811828717765) + " "
						+ (length - 27.845010385562286));
			writer.write("\r\n");
			if (label == 17)
				writer.write((point[0] - 68.57738854955394) + " " + (point[1] - 408.71403968183364) + " "
						+ (point[2] - 83.39252020656194) + " " + (point[3] - 383.1254996023628) + " "
						+ (length - 29.56791351132448));
			writer.write("\r\n");
			if (label == 18)
				writer.write((point[0] - 218.86847126822028) + " " + (point[1] - 423.4828607049002) + " "
						+ (point[2] - 232.1580617020211) + " " + (point[3] - 443.0346629756175 ) + " "
						+ (length - 23.640773801451516));
			writer.write("\r\n");
			if (label == 19)
				writer.write((point[0] - 364.21538301871306) + " " + (point[1] - 424.49188501126974) + " "
						+ (point[2] - 390.83182460376787) + " " + (point[3] - 428.76181034672146) + " "
						+ (length - 26.956765848687557));
			writer.write("\r\n");
			if (label == 20)
				writer.write((point[0] - 310.1599594899791) + " " + (point[1] - 431.5838448697581) + " "
						+ (point[2] - 328.10291270960886) + " " + (point[3] - 432.41153475529126) + " "
						+ (length - 17.96203331442285));
			writer.write("\r\n");

			writer.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	private void Testerrorfour(final double[] point, int label, double length) {

		// Errorlist for Fake_big_file.tif
		try {
			FileWriter writer = new FileWriter("../res/error-Pnoise4snr45.txt", true);

			if (label == 1)
				writer.write((point[0] - 327.76209215919863) + " " + (point[1] - 47.380628285633605) + " "
						+ (point[2] - 329.4156263605521) + " " + (point[3] - 15.258050424764825) + " "
						+ (length - 32.16510817302847));
			writer.write("\r\n");
			if (label == 2)
				writer.write((point[0] - 4.875453301405078) + " " + (point[1] - 77.8324862348326) + " "
						+ (point[2] - 9.864253108433854) + " " + (point[3] - 49.79207039132902 ) + " "
						+ (length - 28.4807486592473));
			writer.write("\r\n");
			if (label == 3)
				writer.write((point[0] - 408.89938816919187) + " " + (point[1] - 69.82732619025738) + " "
						+ (point[2] - 409.99222448920546 ) + " " + (point[3] - 85.44862419856435) + " "
						+ (length - 15.65947772713624));
			writer.write("\r\n");
			if (label == 4)
				writer.write((point[0] - 447.7387612495996) + " " + (point[1] - 96.75498327965913) + " "
						+ (point[2] - 469.40903018475075) + " " + (point[3] - 98.45695235464969) + " "
						+ (length - 21.73700196563459));
			writer.write("\r\n");
			if (label == 5)
				writer.write((point[0] - 434.2810603432197) + " " + (point[1] - 152.33342406129356) + " "
						+ (point[2] - 455.5560388521431) + " " + (point[3] - 130.8829197976539) + " "
						+ (length - 30.21173354376701));
			writer.write("\r\n");
			if (label == 6)
				writer.write((point[0] - 345.10627168817774) + " " + (point[1] -152.1834118866404) + " "
						+ (point[2] - 367.64539043884594) + " " + (point[3] - 144.05117734498975 ) + " "
						+ (length - 23.96132535351788));
			writer.write("\r\n");

			if (label == 7)
				writer.write((point[0] - 3.908286620541236) + " " + (point[1] - 190.1178679247746) + " "
						+ (point[2] - 38.819415473125794) + " " + (point[3] - 173.35128294593778) + " "
						+ (length - 38.728610736951545));
			writer.write("\r\n");
			if (label == 8)
				writer.write((point[0] - 266.1954610919443) + " " + (point[1] - 189.91181955051073) + " "
						+ (point[2] - 282.29275670870487) + " " + (point[3] - 193.15782758353862) + " "
						+ (length - 16.42131219859912));
			writer.write("\r\n");
			if (label == 9)
				writer.write((point[0] - 96.04691366056015) + " " + (point[1] - 212.26565365789554) + " "
						+ (point[2] - 102.66535594773856) + " " + (point[3] - 194.1139650837576 ) + " "
						+ (length - 19.320651552191546));
			writer.write("\r\n");
			if (label == 10)
				writer.write((point[0] - 336.23861322738924) + " " + (point[1] - 220.056673149817) + " "
						+ (point[2] - 356.2924875226813) + " " + (point[3] - 241.45973560520414) + " "
						+ (length - 29.330000966937924));
			writer.write("\r\n");
			if (label == 11)
				writer.write((point[0] - 278.10070833534184) + " " + (point[1] - 239.64339233587503) + " "
						+ (point[2] - 290.55485014203896) + " " + (point[3] - 259.7640723616255 ) + " "
						+ (length - 23.663207999761035));
			writer.write("\r\n");
			if (label == 12)
				writer.write((point[0] - 300.7680991925678 ) + " " + (point[1] - 330.0621629416804) + " "
						+ (point[2] - 314.6509543673136) + " " + (point[3] - 339.5818090159117) + " "
						+ (length - 16.833220998418394));
			writer.write("\r\n");
			if (label == 14)
				writer.write((point[0] - 427.46299173976337) + " " + (point[1] - 359.22598602730744) + " "
						+ (point[2] - 439.1154146009716) + " " + (point[3] - 367.4235431977808) + " "
						+ (length - 14.247066438379065));
			writer.write("\r\n");

			if (label == 15)
				writer.write((point[0] - 362.8999580440501) + " " + (point[1] - 373.13074741573) + " "
						+ (point[2] - 374.34870157838213) + " " + (point[3] - 361.8824868380733) + " "
						+ (length - 16.049831604654788));
			writer.write("\r\n");

			

			if (label == 16)
				writer.write((point[0] - 147.73224268059857 ) + " " + (point[1] - 403.88078930974433) + " "
						+ (point[2] - 177.54144738240447) + " " + (point[3] - 399.40758592299727) + " "
						+ (length - 30.142963250041827));
			writer.write("\r\n");
			
			if (label == 17)
				writer.write((point[0] - 281.66617970983623) + " " + (point[1] - 402.03989076001506) + " "
						+ (point[2] - 293.8247670408476) + " " + (point[3] - 418.9879654976899) + " "
						+ (length - 20.858295309052522));
			writer.write("\r\n");
			if (label == 18)
				writer.write((point[0] - 331.0027901598171) + " " + (point[1] - 442.34452697768046) + " "
						+ (point[2] - 345.36552081520773) + " " + (point[3] - 435.74263790951096) + " "
						+ (length - 15.807370785418865));
			writer.write("\r\n");
			if (label == 19)
				writer.write((point[0] - 14.375693702864483) + " " + (point[1] - 447.54202575273393) + " "
						+ (point[2] - 46.73404009764429) + " " + (point[3] - 467.5466930496815) + " "
						+ (length - 38.042729858228135));
			writer.write("\r\n");
			

			writer.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	private void Testerrorfive(final double[] point, int label, double length) {

		// Errorlist for Fake_big_file.tif
		try {
			FileWriter writer = new FileWriter("../res/error-Pnonoise5.txt", true);

			if (label == 1)
				writer.write((point[0] - 213.23931042861454) + " " + (point[1] - 28.768853068698554) + " "
						+ (point[2] - 231.3531997384953) + " " + (point[3] - 15.92967030628381) + " "
						+ (length - 22.20264848925237));
			writer.write("\r\n");
			if (label == 2)
				writer.write((point[0] - 262.5684393922506) + " " + (point[1] - 20.050913217219318) + " "
						+ (point[2] - 272.2177907332039) + " " + (point[3] - 38.723725344178916 ) + " "
						+ (length - 21.018655856878368));
			writer.write("\r\n");
			if (label == 3)
				writer.write((point[0] - 299.3825946207502) + " " + (point[1] - 94.80260795693071) + " "
						+ (point[2] - 330.2531548768725 ) + " " + (point[3] - 70.57277447839192) + " "
						+ (length - 39.24380614727115));
			writer.write("\r\n");
			if (label == 4)
				writer.write((point[0] - 45.675470197175386) + " " + (point[1] - 87.0860702029411) + " "
						+ (point[2] - 73.43981930651455) + " " + (point[3] - 77.35201386219319 ) + " "
						+ (length - 29.421266701318572));
			writer.write("\r\n");
			if (label == 5)
				writer.write((point[0] - 358.51322277498105) + " " + (point[1] - 78.69665600584534) + " "
						+ (point[2] - 379.0787208734989) + " " + (point[3] - 98.10097272030382) + " "
						+ (length - 28.274851355845446));
			writer.write("\r\n");
			if (label == 6)
				writer.write((point[0] - 49.76148246466852) + " " + (point[1] -124.09604469489028) + " "
						+ (point[2] - 84.2592830009717) + " " + (point[3] - 115.61717182871169 ) + " "
						+ (length - 35.52449193054533));
			writer.write("\r\n");

			if (label == 7)
				writer.write((point[0] - 223.4539445808136) + " " + (point[1] - 132.57748562595273) + " "
						+ (point[2] - 236.38132216426166) + " " + (point[3] - 117.365369486839) + " "
						+ (length - 19.963105184688068));
			writer.write("\r\n");
			if (label == 8)
				writer.write((point[0] - 370.0214563609312) + " " + (point[1] - 142.97032407028482) + " "
						+ (point[2] - 398.32465146782636) + " " + (point[3] - 128.4247044367045) + " "
						+ (length - 31.82209772758186));
			writer.write("\r\n");
			if (label == 9)
				writer.write((point[0] - 327.4805366491584) + " " + (point[1] -131.2586576595218) + " "
						+ (point[2] - 345.48254113777244) + " " + (point[3] - 157.9396065221859 ) + " "
						+ (length - 32.186102557162414));
			writer.write("\r\n");
			if (label == 10)
				writer.write((point[0] - 153.68916577838593) + " " + (point[1] - 138.1918303951238) + " "
						+ (point[2] - 186.84440518521032) + " " + (point[3] - 140.40210558701762) + " "
						+ (length - 33.22883110414422));
			writer.write("\r\n");
			if (label == 11)
				writer.write((point[0] - 117.3400876334884) + " " + (point[1] - 199.9830506291008) + " "
						+ (point[2] - 125.52613128612809) + " " + (point[3] - 189.40155469545007  ) + " "
						+ (length - 13.378317042019507));
			writer.write("\r\n");
			if (label == 12)
				writer.write((point[0] - 355.9814722371587) + " " + (point[1] - 212.54284215973115) + " "
						+ (point[2] - 369.88081057901996) + " " + (point[3] - 189.87436359988837) + " "
						+ (length - 26.59044051082243));
			writer.write("\r\n");
			if (label == 13)
				writer.write((point[0] - 16.519680837702623) + " " + (point[1] - 208.35089834691382) + " "
						+ (point[2] - 44.309044477124075) + " " + (point[3] - 220.29558511100268 ) + " "
						+ (length - 30.247715176128573));
			writer.write("\r\n");
			if (label == 14)
				writer.write((point[0] - 322.0314627356028) + " " + (point[1] - 230.46786868249026) + " "
						+ (point[2] - 344.2062978545556) + " " + (point[3] - 211.85766750777054 ) + " "
						+ (length - 28.9493160595597));
			writer.write("\r\n");

			if (label == 15)
				writer.write((point[0] - 140.77986413780962) + " " + (point[1] - 222.04591247544067) + " "
						+ (point[2] - 169.1431556561999) + " " + (point[3] - 214.4000255926368) + " "
						+ (length - 29.375770491713432));
			writer.write("\r\n");

			

			if (label == 16)
				writer.write((point[0] - 150.2821689064564) + " " + (point[1] - 292.25054302689387) + " "
						+ (point[2] - 178.12273121820238) + " " + (point[3] - 269.81238539563384) + " "
						+ (length - 35.757066822091026));
			writer.write("\r\n");
			
			if (label == 17)
				writer.write((point[0] - 405.3441866021064) + " " + (point[1] -283.53349253317094) + " "
						+ (point[2] - 423.8493855332569) + " " + (point[3] - 304.5534970808495) + " "
						+ (length - 28.00505273456707));
			writer.write("\r\n");
			if (label == 18)
				writer.write((point[0] -359.5585400563441) + " " + (point[1] - 324.0239934624044) + " "
						+ (point[2] - 372.9916314475815) + " " + (point[3] - 336.8607944237224) + " "
						+ (length - 18.58040374281006));
			writer.write("\r\n");
			if (label == 19)
				writer.write((point[0] - 46.42913218018552) + " " + (point[1] - 388.9428935558981) + " "
						+ (point[2] - 60.6517127858811 ) + " " + (point[3] - 353.42924286364405) + " "
						+ (length - 38.25573400912526));
			writer.write("\r\n");
			if (label == 20)
				writer.write((point[0] - 181.55286568938038) + " " + (point[1] - 427.8287110953466) + " "
						+ (point[2] - 193.34300560337857) + " " + (point[3] - 415.2158643110275) + " "
						+ (length - 17.265320825179078));
			writer.write("\r\n");

			writer.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private PointSampleList<FloatType> gatherfullData(final int label, final boolean offsetting) {
		final PointSampleList<FloatType> datalist = new PointSampleList<FloatType>(ndims);

		RandomAccessibleInterval<FloatType> currentimg = PerformWatershedding.CurrentLabelImage(intimg, inputimg,
				label);

		boolean outofbounds = false;
		final long[] minCorner = PerformWatershedding.GetMincorners(intimg, label);
		final long[] maxCorner = PerformWatershedding.GetMaxcorners(intimg, label);
		FinalInterval smallinterval = new FinalInterval(minCorner, maxCorner);

		if (offsetting) {

			currentimg = Views.offsetInterval(currentimg, smallinterval);

		}

		Cursor<FloatType> localcursor = Views.iterable(currentimg).localizingCursor();

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

			Point newpoint = new Point(localcursor);
			datalist.add(newpoint, localcursor.get().copy());

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

	public double sqDistance(final double[] cordone, final double[] cordtwo) {

		double distance = 0;

		for (int d = 0; d < ndims; ++d) {

			distance += Math.pow((cordone[d] - cordtwo[d]), 2);

		}
		return (distance);
	}

	private void Testmovingline(final double[] point, int label, double length, int rate) {
		try {
			FileWriter writer = new FileWriter("error-Fake_movingline.txt", true);

			if (rate == 0) {

				writer.write(" Rate :" + rate + " ");
				if (label == 1)
					writer.write((point[0] - 169.73720446120004) + " " + (point[1] - 104.95635903304128) + " "
							+ (point[2] - 181.43458804585626) + " " + (point[3] - 60.78333591544249) + " "
							+ (length - 45.69556602203874));
				writer.write("\r\n");
				if (label == 2)
					writer.write((point[0] - 294.6604728889232) + " " + (point[1] - 109.07711973280686) + " "
							+ (point[2] - 305.2140902032087) + " " + (point[3] - 82.10337364262949) + " "
							+ (length - 28.964837588941272));
				writer.write("\r\n");
				if (label == 3)
					writer.write((point[0] - 391.2990529717664) + " " + (point[1] - 124.0347325138523) + " "
							+ (point[2] - 401.71221291284917) + " " + (point[3] - 82.7441754252086) + " "
							+ (length - 42.58337709305259));
				writer.write("\r\n");
				if (label == 4)
					writer.write((point[0] - 203.44041398569848) + " " + (point[1] - 130.38479493829118) + " "
							+ (point[2] - 214.9881816669676) + " " + (point[3] - 99.63148998257181) + " "
							+ (length - 32.84991178253051));
				writer.write("\r\n");
				if (label == 5)
					writer.write((point[0] - 132.08290714831938) + " " + (point[1] - 221.96874365412077) + " "
							+ (point[2] - 143.00136725366892) + " " + (point[3] - 177.4220897013742) + " "
							+ (length - 45.86520630562846));
				writer.write("\r\n");
				if (label == 6)
					writer.write((point[0] - 4.561985449277813) + " " + (point[1] - 228.0886059139005) + " "
							+ (point[2] - 15.44169613220083) + " " + (point[3] - 198.4982244474968) + " "
							+ (length - 31.52711182254722));
				writer.write("\r\n");

				if (label == 7)
					writer.write((point[0] - 360.6948730664424) + " " + (point[1] - 246.8371168164922) + " "
							+ (point[2] - 370.8155639541149) + " " + (point[3] - 212.5197210031867) + " "
							+ (length - 35.778653404661476));
				writer.write("\r\n");
				if (label == 8)
					writer.write((point[0] - 303.8536184864458) + " " + (point[1] - 341.40556128106505) + " "
							+ (point[2] - 315.5894071390266) + " " + (point[3] - 303.15203819534213) + " "
							+ (length - 40.0132573501331));
				writer.write("\r\n");
				if (label == 9)
					writer.write((point[0] - 104.56770264695977) + " " + (point[1] - 377.49548375649414) + " "
							+ (point[2] - 115.11670119664012) + " " + (point[3] - 339.5842911894899) + " "
							+ (length - 39.35149161408799));
				writer.write("\r\n");
				if (label == 10)
					writer.write((point[0] - 181.75417277692245) + " " + (point[1] - 377.86831136750124) + " "
							+ (point[2] - 193.23751455228975) + " " + (point[3] - 351.97062628010985) + " "
							+ (length - 28.329441067828885));
				writer.write("\r\n");
			}

			if (rate == 1) {

				writer.write(" Rate :" + rate + " ");
				if (label == 1)
					writer.write((point[0] - 171.1859449967859) + " " + (point[1] - 99.48545618483672) + " "
							+ (point[2] - 181.43458804585626) + " " + (point[3] - 60.78333591544249) + " "
							+ (length - 40.036093686746135));
				writer.write("\r\n");
				if (label == 2)
					writer.write((point[0] - 298.63632549289247) + " " + (point[1] - 98.91533013614685) + " "
							+ (point[2] - 305.2140902032087) + " " + (point[3] - 82.10337364262949) + " "
							+ (length - 18.05294628929588));
				writer.write("\r\n");
				if (label == 3)
					writer.write((point[0] - 392.2658523406384) + " " + (point[1] - 120.20115224047345) + " "
							+ (point[2] - 401.71221291284917) + " " + (point[3] - 82.7441754252086) + " "
							+ (length - 38.62976624572697));
				writer.write("\r\n");
				if (label == 4)
					writer.write((point[0] - 206.0628373691796) + " " + (point[1] - 123.40091820460714) + " "
							+ (point[2] - 214.9881816669676) + " " + (point[3] - 99.63148998257181) + " "
							+ (length - 25.38990919315283));
				writer.write("\r\n");
				if (label == 5)
					writer.write((point[0] - 135.10455586538396) + " " + (point[1] - 209.6406004482957) + " "
							+ (point[2] - 143.00136725366892) + " " + (point[3] - 177.4220897013742) + " "
							+ (length - 33.17215797700901));
				writer.write("\r\n");
				if (label == 6)
					writer.write((point[0] - 4.798655006524217) + " " + (point[1] - 227.44491761331008) + " "
							+ (point[2] - 15.44169613220083) + " " + (point[3] - 198.4982244474968) + " "
							+ (length - 30.841293254962338));
				writer.write("\r\n");

				if (label == 7)
					writer.write((point[0] - 363.616204272503) + " " + (point[1] - 236.93142160133107) + " "
							+ (point[2] - 370.8155639541149) + " " + (point[3] - 212.5197210031867) + " "
							+ (length - 25.45116708362626));
				writer.write("\r\n");
				if (label == 8)
					writer.write((point[0] - 307.02690788497586) + " " + (point[1] - 331.06202981011825) + " "
							+ (point[2] - 315.5894071390266) + " " + (point[3] - 303.15203819534213) + " "
							+ (length - 29.19390390839315));
				writer.write("\r\n");
				if (label == 9)
					writer.write((point[0] - 105.70099297153607) + " " + (point[1] - 373.42264333216434) + " "
							+ (point[2] - 115.11670119664012) + " " + (point[3] - 339.5842911894899) + " "
							+ (length - 35.123918305222354));
				writer.write("\r\n");
				if (label == 10)
					writer.write((point[0] - 182.7739990645031) + " " + (point[1] - 375.56835891844486) + " "
							+ (point[2] - 193.23751455228975) + " " + (point[3] - 351.97062628010985) + " "
							+ (length - 25.813526338597992));
				writer.write("\r\n");
			}

			if (rate == 2) {

				writer.write(" Rate :" + rate + " ");
				if (label == 1)
					writer.write((point[0] - 172.43796447094846) + " " + (point[1] - 94.75743433485547) + " "
							+ (point[2] - 181.43458804585626) + " " + (point[3] - 60.78333591544249) + " "
							+ (length - 35.14510775571112));
				writer.write("\r\n");
				if (label == 2)
					writer.write((point[0] - 298.63632549289247) + " " + (point[1] - 98.91533013614685) + " "
							+ (point[2] - 305.2140902032087) + " " + (point[3] - 82.10337364262949) + " "
							+ (length - 18.05294628929588));
				writer.write("\r\n");
				if (label == 3)
					writer.write((point[0] - 393.1753225631553) + " " + (point[1] - 116.59489512392793) + " "
							+ (point[2] - 401.71221291284917) + " " + (point[3] - 82.7441754252086) + " "
							+ (length - 34.91059611298497));
				writer.write("\r\n");
				if (label == 4)
					writer.write((point[0] - 207.40081502209952) + " " + (point[1] - 119.83769837546845) + " "
							+ (point[2] - 214.9881816669676) + " " + (point[3] - 99.63148998257181) + " "
							+ (length - 21.583766821869236));
				writer.write("\r\n");
				if (label == 5)
					writer.write((point[0] - 136.04277722443922) + " " + (point[1] - 205.81271429863318) + " "
							+ (point[2] - 143.00136725366892) + " " + (point[3] - 177.4220897013742) + " "
							+ (length - 29.23096885526344));
				writer.write("\r\n");
				if (label == 6)
					writer.write((point[0] - 5.034495309332331) + " " + (point[1] - 226.8034846994553) + " "
							+ (point[2] - 15.44169613220083) + " " + (point[3] - 198.4982244474968) + " "
							+ (length - 30.157877692215305));
				writer.write("\r\n");

				if (label == 7)
					writer.write((point[0] - 364.68594610462594) + " " + (point[1] - 233.3041243863031) + " "
							+ (point[2] - 370.8155639541149) + " " + (point[3] - 212.5197210031867) + " "
							+ (length - 21.669417135051745));
				writer.write("\r\n");
				if (label == 8)
					writer.write((point[0] - 307.7174271751884) + " " + (point[1] - 328.81123982137103) + " "
							+ (point[2] - 315.5894071390266) + " " + (point[3] - 303.15203819534213) + " "
							+ (length - 26.839573331859704));
				writer.write("\r\n");
				if (label == 9)
					writer.write((point[0] - 106.74140940290457) + " " + (point[1] - 369.6835748799574) + " "
							+ (point[2] - 115.11670119664012) + " " + (point[3] - 339.5842911894899) + " "
							+ (length - 31.24279743091932));
				writer.write("\r\n");
				if (label == 10)
					writer.write((point[0] - 184.54848516438193) + " " + (point[1] - 371.56646790922724) + " "
							+ (point[2] - 193.23751455228975) + " " + (point[3] - 351.97062628010985) + " "
							+ (length - 21.435863426915553));
				writer.write("\r\n");
			}

			writer.close();

		} catch (IOException e) {
			e.printStackTrace();
		}

	}

}

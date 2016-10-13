package peakFitter;

import houghandWatershed.Boundingboxes;
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
	final int iterations = 1500;
	final double cutoffdistance = 15;
	final boolean halfgaussian = false;
	final double Intensityratio = 0.5;
	public Linefitter(RandomAccessibleInterval<FloatType> inputimg, RandomAccessibleInterval<IntType> intimg) {

		this.inputimg = inputimg;
		this.intimg = intimg;
		this.ndims = inputimg.numDimensions();
		assert (inputimg.numDimensions() == intimg.numDimensions());

	}

	private final double[] MakeimprovedLineguess(double slope, double intercept, double[] psf, int label,
			boolean offsetting) throws Exception {
		long[] newposition = new long[ndims];
		double[] minVal = { Double.MAX_VALUE, Double.MAX_VALUE };
		double[] maxVal = { -Double.MIN_VALUE, -Double.MIN_VALUE };

		RandomAccessibleInterval<FloatType> currentimg = Boundingboxes.CurrentLabelImage(intimg, inputimg,
				label);

		double newintercept = intercept;

		final long[] minCorner = Boundingboxes.GetMincorners(intimg, label);
		final long[] maxCorner = Boundingboxes.GetMaxcorners(intimg, label);
		FinalInterval smallinterval = new FinalInterval(minCorner, maxCorner);
		
		if (offsetting == false)
		currentimg = Views.interval(currentimg, smallinterval);
		
		if (offsetting) {

			currentimg = Views.offsetInterval(currentimg, smallinterval);

			newintercept = intercept - (smallinterval.realMin(1) - slope * smallinterval.realMin(0));

		}

		final Cursor<FloatType> outcursor = Views.iterable(currentimg).localizingCursor();

		final double maxintensityline = GetLocalmaxmin.computeMaxIntensity(currentimg);

		while (outcursor.hasNext()) {

			outcursor.fwd();

			if (outcursor.get().get() / maxintensityline > Intensityratio) {
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
		MinandMax[2 * ndims] =  0.5 * Math.min(psf[0], psf[1]);
		MinandMax[2 * ndims + 1] = maxintensityline; 
		// This parameter guess estimates the background noise level
		MinandMax[2 * ndims + 2] = 1; 
		
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

		RandomAccessibleInterval<FloatType> currentimg = Boundingboxes.CurrentLabelImage(intimg, inputimg,
				label);
		final double[] cordone = { iniparam[0], iniparam[1] };
		final double[] cordtwo = { iniparam[2], iniparam[3] };

		double slope = (cordone[1] - cordtwo[1]) / (cordone[0] - cordtwo[0]);
		double intercept = cordone[1] - slope * cordone[0];
		double newintercept = intercept;

		final long[] minCorner = Boundingboxes.GetMincorners(intimg, label);
		final long[] maxCorner = Boundingboxes.GetMaxcorners(intimg, label);
		FinalInterval smallinterval = new FinalInterval(minCorner, maxCorner);
		if (offsetting == false)
		currentimg = Views.interval(currentimg, smallinterval);
		
		
		if (offsetting) {

			currentimg = Views.offsetInterval(currentimg, smallinterval);

			newintercept = intercept - (smallinterval.realMin(1) - slope * smallinterval.realMin(0));

		}

		final Cursor<FloatType> outcursor = Views.iterable(currentimg).localizingCursor();

		final double maxintensityline = GetLocalmaxmin.computeMaxIntensity(currentimg);

		while (outcursor.hasNext()) {

			outcursor.fwd();

			if (outcursor.get().get() / maxintensityline > Intensityratio) {
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

			RandomAccessibleInterval<FloatType> currentimg = Boundingboxes.CurrentLabelImage(intimg, inputimg,
					label);

			final long[] minCorner = Boundingboxes.GetMincorners(intimg, label);
			final long[] maxCorner = Boundingboxes.GetMaxcorners(intimg, label);
			FinalInterval smallinterval = new FinalInterval(minCorner, maxCorner);
			
			if (offsetting == false)
			currentimg = Views.interval(currentimg, smallinterval);
			
			if (offsetting) {

				currentimg = Views.offsetInterval(currentimg, smallinterval);

			}

			final double[] fixed_param = new double[ndims];

			for (int d = 0; d < ndims; ++d) {

				fixed_param[d] = 1.0 / Math.pow(psf[d], 2);
			}

		
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
				final double LMdist = sqDistance(startpos, endpos);

				double[] dxvector = { dx, dy };

				double[] startfit = new double[ndims];
				double[] endfit = new double[ndims];
				final double maxintensityline = GetLocalmaxmin.computeMaxIntensity(currentimg);

				

				System.out.println("Doing Mask Fits: ");
				startfit = peakFitter.GaussianMaskFit.sumofgaussianMaskFit(currentimg, startpos.clone(), psf,
						iterations, dxvector, newslope, newintercept, maxintensityline, halfgaussian, Endfit.Start,
						label);

				endfit = peakFitter.GaussianMaskFit.sumofgaussianMaskFit(currentimg, endpos.clone(), psf,
						iterations, dxvector, newslope, newintercept, maxintensityline,  halfgaussian, Endfit.End,
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
					System.out.println("Mask fits fail, both cords move far, returning LM solver results!");

					for (int d = 0; d < ndims; ++d) {
						returnparam[d] = startpos[d];
						returnparam[ndims + d] = endpos[d];
					}
					}
				if (Math.abs(startpos[0] - startfit[0]) >= cutoffdistance || Math.abs(startpos[1] - startfit[1]) >= cutoffdistance 
						|| Math.abs(endpos[0] - endfit[0]) >= cutoffdistance  || Math.abs(endpos[1] - endfit[1]) >= cutoffdistance  ){
					System.out.println("Mask fits fail, one cord moves too much, returning LM solver results!");
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
				RandomAccessibleInterval<FloatType> currentimg = Boundingboxes.CurrentLabelImage(intimg,
						inputimg, label);

				final long[] minCorner = Boundingboxes.GetMincorners(intimg, label);
				final long[] maxCorner = Boundingboxes.GetMaxcorners(intimg, label);
				FinalInterval smallinterval = new FinalInterval(minCorner, maxCorner);
				if (offsetting == false)
				currentimg = Views.interval(currentimg, smallinterval);
				
				if (offsetting) {

					currentimg = Views.offsetInterval(currentimg, smallinterval);

				}

				final double[] fixed_param = new double[ndims];

				for (int d = 0; d < ndims; ++d) {

					fixed_param[d] = 1.0 / Math.pow(psf[d], 2);
				}

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

					final double maxintensityline = GetLocalmaxmin.computeMaxIntensity(currentimg);

					

					double newslope = (endpos[1] - startpos[1]) / (endpos[0] - startpos[0]);
					double newintercept = (endpos[1] - newslope * endpos[0]);
					double dx = finalparamstart[4] / Math.sqrt(1 + newslope * newslope);
					double dy = newslope * dx;
					double[] dxvector = { dx,  dy };

					double[] startfit = new double[ndims];
					double[] endfit = new double[ndims];

					
					

					startfit = peakFitter.GaussianMaskFit.sumofgaussianMaskFit(currentimg,  startpos.clone(),
							psf, iterations, dxvector, newslope, newintercept, maxintensityline,  halfgaussian,
							Endfit.Start, label);

					endfit = peakFitter.GaussianMaskFit.sumofgaussianMaskFit(currentimg,  endpos.clone(), psf,
							iterations, dxvector, newslope, newintercept, maxintensityline,  halfgaussian,
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



	private PointSampleList<FloatType> gatherfullData(final int label, final boolean offsetting) {
		final PointSampleList<FloatType> datalist = new PointSampleList<FloatType>(ndims);

		RandomAccessibleInterval<FloatType> currentimg = Boundingboxes.CurrentLabelImage(intimg, inputimg,
				label);

		boolean outofbounds = false;
		final long[] minCorner = Boundingboxes.GetMincorners(intimg, label);
		final long[] maxCorner = Boundingboxes.GetMaxcorners(intimg, label);
		FinalInterval smallinterval = new FinalInterval(minCorner, maxCorner);
		
		if (offsetting == false)
			currentimg = Views.interval(currentimg, smallinterval);
			
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

	

}

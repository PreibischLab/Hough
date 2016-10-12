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

		RandomAccessibleInterval<FloatType> currentimg = PerformWatershedding.CurrentLabelImage(intimg, inputimg,
				label);

		double newintercept = intercept;

		final long[] minCorner = PerformWatershedding.GetMincorners(intimg, label);
		final long[] maxCorner = PerformWatershedding.GetMaxcorners(intimg, label);
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

			RandomAccessibleInterval<FloatType> currentimg = PerformWatershedding.CurrentLabelImage(intimg, inputimg,
					label);

			final long[] minCorner = PerformWatershedding.GetMincorners(intimg, label);
			final long[] maxCorner = PerformWatershedding.GetMaxcorners(intimg, label);
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
				Testerrorone(returnparam, label, Distance(new double[] { returnparam[0], returnparam[1] },
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
					double ds = finalparamstart[4];
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

		try {
			FileWriter writer = new FileWriter("../res/error-Pnonoise1mask.txt", true);

			if (label == 1)
				writer.write((point[0] - 331.5364522249255) + " " + (point[1] - 19.893563661868367) + " "
						+ (point[2] - 336.0544738270745) + " " + (point[3] - 4.544704394726786 ) + " "
						+ (length - 16));
			writer.write("\r\n");
			if (label == 2)
				writer.write((point[0] - 358.564293824347) + " " + (point[1] - 134.43573787129) + " "
						+ (point[2] - 376.66554377520197) + " " + (point[3] - 120.24613779881884) + " "
						+ (length - 23));
			writer.write("\r\n");
			if (label == 3)
				writer.write((point[0] - 30.497372093540314) + " " + (point[1] - 134.31506624504095) + " "
						+ (point[2] - 44.831239399479294) + " " + (point[3] - 138.7355015240304) + " "
						+ (length - 15));
			writer.write("\r\n");
			if (label == 4)
				writer.write((point[0] - 244.98542147699598) + " " + (point[1] - 137.38834701507116) + " "
						+ (point[2] - 268.95343229596267) + " " + (point[3] - 138.62707918220516) + " "
						+ (length - 23.999));
			writer.write("\r\n");
			if (label == 5)
				writer.write((point[0] - 350.99495980590706) + " " + (point[1] - 148.89732806584144) + " "
						+ (point[2] - 365.03886102752364) + " " + (point[3] - 141.23124441918836) + " "
						+ (length - 16));
			writer.write("\r\n");
			if (label == 6)
				writer.write((point[0] - 399.07801551170985) + " " + (point[1] - 215.6920990720688 ) + " "
						+ (point[2] - 408.899852747816) + " " + (point[3] - 203.06152056958365) + " "
						+ (length - 16));
			writer.write("\r\n");

			if (label == 7)
				writer.write((point[0] - 108.78034293666624) + " " + (point[1] - 216.68419139490626) + " "
						+ (point[2] - 128.92832030710323) + " " + (point[3] - 210.7631266358985) + " "
						+ (length - 21));
			writer.write("\r\n");
			if (label == 8)
				writer.write((point[0] - 149.75242608356956) + " " + (point[1] - 258.88973927032725) + " "
						+ (point[2] - 164.23612622612032) + " " + (point[3] - 238.51273188011726) + " "
						+ (length - 24.999));
			writer.write("\r\n");
			
			if (label == 9)
				writer.write((point[0] - 4.779577587106245) + " " + (point[1] - 242.66186450394378) + " "
						+ (point[2] - 19.72450661678264) + " " + (point[3] - 260.1446873966559) + " "
						+ (length - 23));
			writer.write("\r\n");
			if (label == 10)
				writer.write((point[0] - 407.3501916855185) + " " + (point[1] - 284.8041109947676) + " "
						+ (point[2] - 423.02540981055705) + " " + (point[3] - 272.3828575663706) + " "
						+ (length - 20));
			writer.write("\r\n");
			if (label == 11)
				writer.write((point[0] - 188.9178728935583) + " " + (point[1] - 360.25743858520144) + " "
						+ (point[2] - 201.19360936039828) + " " + (point[3] - 343.21906148535084) + " "
						+ (length - 21));
			writer.write("\r\n");
			if (label == 13)
				writer.write((point[0] - 242.9832929183822 ) + " " + (point[1] - 401.7231934614111) + " "
						+ (point[2] - 257.94695655884016) + " " + (point[3] - 416.4571260815654 ) + " "
						+ (length - 21));
			writer.write("\r\n");
			

			if (label == 14)
				writer.write((point[0] - 128.9864315967172) + " " + (point[1] - 423.679581484538) + " "
						+ (point[2] - 149.15996670775016) + " " + (point[3] - 438.44536900582943) + " "
						+ (length - 25));
			writer.write("\r\n");

			if (label == 15)
				writer.write((point[0] -384.1527617353957) + " " + (point[1] - 429.634979413752) + " "
						+ (point[2] - 395.76151966251473) + " " + (point[3] - 447.13460053506606) + " "
						+ (length - 20.999));
			writer.write("\r\n");
			if (label == 16)
				writer.write((point[0] - 452.69591108476527) + " " + (point[1] - 431.5437825330961) + " "
						+ (point[2] - 480.50980808810493) + " " + (point[3] - 439.7527491847318) + " "
						+ (length - 29));
			writer.write("\r\n");
			if (label == 17)
				writer.write((point[0] - 212.21926128834326) + " " + (point[1] - 431.3941152486386) + " "
						+ (point[2] - 218.0397374553543) + " " + (point[3] - 447.36665597293324 ) + " "
						+ (length - 16.999));
			writer.write("\r\n");
			
			if (label == 18)
				writer.write((point[0] - 132.78463690710356) + " " + (point[1] -481.99168690306874) + " "
						+ (point[2] - 143.55251539520583) + " " + (point[3] - 503.4405281973079) + " "
						+ (length - 23.999));
			writer.write("\r\n");
			if (label == 19)
				writer.write((point[0] - 16.622432112826097) + " " + (point[1] -506.4063877219946) + " "
						+ (point[2] - 32.41110113306011) + " " + (point[3] - 503.81449062343876) + " "
						+ (length - 15.999));
			writer.write("\r\n");
			if (label == 20)
				writer.write((point[0] - 108.30819273568439) + " " + (point[1] - 524.9816812046726) + " "
						+ (point[2] - 129.18791411562972) + " " + (point[3] - 515.3361000870906 ) + " "
						+ (length - 23));
			writer.write("\r\n");
			
			

			writer.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	
	}

	private void Testerrorthree(final double[] point, int label, double length) {

		// Errorlist for Fake_big_file.tif
		try {
			FileWriter writer = new FileWriter("../res/error-Pnonoise3mask.txt", true);

			if (label == 1)
				writer.write((point[0] - 121.11739066464122) + " " + (point[1] - 34.945911412857576) + " "
						+ (point[2] - 128.3358379749917) + " " + (point[3] - 11.010706829907168) + " "
						+ (length - 24.99999));
			writer.write("\r\n");
			if (label == 2)
				writer.write((point[0] - 383.9612604373775) + " " + (point[1] - 23.687713327920708) + " "
						+ (point[2] - 397.58964339022407) + " " + (point[3] - 15.305161324469086 ) + " "
						+ (length - 16));
			writer.write("\r\n");
			if (label == 3)
				writer.write((point[0] - 461.6064166590004) + " " + (point[1] - 32.56142533075916) + " "
						+ (point[2] - 469.6329721344923) + " " + (point[3] - 19.889629593814597) + " "
						+ (length - 14.99999999));
			writer.write("\r\n");
			if (label == 4)
				writer.write((point[0] - 131.3943799424747) + " " + (point[1] - 68.63943112705725) + " "
						+ (point[2] - 145.10413223838316) + " " + (point[3] - 81.79400037066326) + " "
						+ (length - 18.99999999));
			writer.write("\r\n");
			if (label == 5)
				writer.write((point[0] - 65.9818008542298) + " " + (point[1] - 71.28283726838842) + " "
						+ (point[2] - 76.18895202493668) + " " + (point[3] - 93.0041186527437) + " "
						+ (length - 24));
			writer.write("\r\n");
			if (label == 6)
				writer.write((point[0] - 365.0494095855029) + " " + (point[1] - 121.57565242259585) + " "
						+ (point[2] - 365.5585088819804 ) + " " + (point[3] - 103.58285336557503) + " "
						+ (length - 18));
			writer.write("\r\n");
			if (label == 7)
				writer.write((point[0] - 456.2305259667868) + " " + (point[1] - 140.8247176153394) + " "
						+ (point[2] - 457.9747823668729 ) + " " + (point[3] - 158.74000632668202) + " "
						+ (length - 18));
			writer.write("\r\n");
			
			if (label == 8)
				writer.write((point[0] - 188.7537677873778) + " " + (point[1] - 162.262291516167) + " "
						+ (point[2] - 205.57157794640466) + " " + (point[3] - 159.78010267263954) + " "
						+ (length - 16.999999));
			writer.write("\r\n");
			if (label == 9)
				writer.write((point[0] - 227.5717307456629) + " " + (point[1] - 210.79546138178526) + " "
						+ (point[2] - 237.03830634648833) + " " + (point[3] - 193.17774454918006) + " "
						+ (length - 20));
			writer.write("\r\n");
			if (label == 10)
				writer.write((point[0] - 365.3800199779222) + " " + (point[1] - 226.63329912891996) + " "
						+ (point[2] - 385.708753416097) + " " + (point[3] - 239.39035979349333) + " "
						+ (length - 23.99999));
			writer.write("\r\n");
			if (label == 11)
				writer.write((point[0] - 303.2754986723827) + " " + (point[1] - 326.276934516481) + " "
						+ (point[2] - 315.1970359599219) + " " + (point[3] - 314.157654475909) + " "
						+ (length -16.99999999));
			writer.write("\r\n");
			if (label == 12)
				writer.write((point[0] - 368.77628705423024) + " " + (point[1] - 323.9212578102178) + " "
						+ (point[2] - 382.4586555012279) + " " + (point[3] - 317.77367343254616) + " "
						+ (length - 14.9999999));
			writer.write("\r\n");
			if (label == 13)
				writer.write((point[0] - 2.836157079979025) + " " + (point[1] - 389.02094856916733) + " "
						+ (point[2] - 23.40689682082701) + " " + (point[3] - 368.57969942645775) + " "
						+ (length - 29));
			writer.write("\r\n");
			if (label == 14)
				writer.write((point[0] - 86.67754892229159) + " " + (point[1] - 432.4207567280493) + " "
						+ (point[2] - 99.25674261766153) + " " + (point[3] - 418.1812382679907) + " "
						+ (length - 18.99999999));
			writer.write("\r\n");

			if (label == 15)
				writer.write((point[0] - 256.28327884809266) + " " + (point[1] - 428.06197407612234) + " "
						+ (point[2] - 270.3455272306394) + " " + (point[3] - 448.73207737599887) + " "
						+ (length - 25));
			writer.write("\r\n");

			

			if (label == 16)
				writer.write((point[0] - 107.94540348847204 ) + " " + (point[1] - 469.4207201182959) + " "
						+ (point[2] - 117.59361709892363) + " " + (point[3] - 443.1355156464434 ) + " "
						+ (length - 27.999));
			writer.write("\r\n");
			if (label == 17)
				writer.write((point[0] - 55.20799587346018) + " " + (point[1] - 519.6763073216334) + " "
						+ (point[2] - 63.19008609304531) + " " + (point[3] - 501.33819848439543) + " "
						+ (length - 19.999999));
			writer.write("\r\n");
			if (label == 18)
				writer.write((point[0] - 328.2849954773905) + " " + (point[1] - 529.5641087001159) + " "
						+ (point[2] - 333.4959718939446) + " " + (point[3] - 504.0916583902992 ) + " "
						+ (length - 25.9999999));
			writer.write("\r\n");
			if (label == 19)
				writer.write((point[0] - 235.2562094745425) + " " + (point[1] - 529.4871437456533) + " "
						+ (point[2] - 238.0281000059525) + " " + (point[3] - 513.7290784176422) + " "
						+ (length - 15.999999));
			writer.write("\r\n");
			

			writer.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}


	private void Testerrortwo(final double[] point, int label, double length) {

		// Errorlist for Fake_big_file.tif
		try {
			FileWriter writer = new FileWriter("../res/error-Pnonoise2mask.txt", true);

			if (label == 1)
				writer.write((point[0] - 33.38153502099174) + " " + (point[1] - 86.98929700953852) + " "
						+ (point[2] - 39.082232864228075) + " " + (point[3] - 107.2007305089249) + " "
						+ (length - 20.999));
			writer.write("\r\n");
			if (label == 2)
				writer.write((point[0] - 374.21126941619156) + " " + (point[1] - 92.18925779701563) + " "
						+ (point[2] - 397.04372678806374 ) + " " + (point[3] - 89.41817671526472 ) + " "
						+ (length - 22.999999));
			writer.write("\r\n");
			if (label == 3)
				writer.write((point[0] - 121.17237814080222) + " " + (point[1] - 96.11726976873652) + " "
						+ (point[2] - 139.29907998413208) + " " + (point[3] - 106.72022597363654 ) + " "
						+ (length - 21));
			writer.write("\r\n");
			if (label == 4)
				writer.write((point[0] - 338.3957301359977) + " " + (point[1] - 125.33087245209133) + " "
						+ (point[2] - 348.06166101122255) + " " + (point[3] - 110.14634808300484) + " "
						+ (length - 18));
			writer.write("\r\n");
			if (label == 5)
				writer.write((point[0] -364.5130007285168) + " " + (point[1] - 140.53413207134295) + " "
						+ (point[2] - 365.3049150278223) + " " + (point[3] - 112.54533303060956) + " "
						+ (length - 28));
			writer.write("\r\n");
			
			if (label == 6)
				writer.write((point[0] - 177.8170530463071) + " " + (point[1] - 146.18901103211232) + " "
						+ (point[2] - 188.48711985814245) + " " + (point[3] - 122.47931524682264) + " "
						+ (length - 26));
			writer.write("\r\n");

			
			if (label == 8)
				writer.write((point[0] - 299.2682363380778) + " " + (point[1] - 153.59115200001608) + " "
						+ (point[2] - 320.4984566918457) + " " + (point[3] - 138.58189673099972) + " "
						+ (length - 26));
			writer.write("\r\n");
			if (label == 9)
				writer.write((point[0] - 78.33703419838355) + " " + (point[1] - 142.06250076854613) + " "
						+ (point[2] - 98.98319910448109) + " " + (point[3] -145.9012342406997) + " "
						+ (length - 21));
			writer.write("\r\n");
			if (label == 10)
				writer.write((point[0] - 357.6922528755881) + " " + (point[1] - 161.16774810283246) + " "
						+ (point[2] - 365.3377270145437 ) + " " + (point[3] - 178.56161843667596  ) + " "
						+ (length - 18.999999));
			writer.write("\r\n");
			if (label == 11)
				writer.write((point[0] - 382.5648799245133) + " " + (point[1] - 245.67237981995828) + " "
						+ (point[2] - 400.79781881539213) + " " + (point[3] - 227.1370098640744) + " "
						+ (length - 26.000));
			writer.write("\r\n");
			if (label == 12)
				writer.write((point[0] - 379.3914908119228) + " " + (point[1] - 275.35463578297276) + " "
						+ (point[2] - 404.292742464066) + " " + (point[3] - 277.57447049677717) + " "
						+ (length - 24.999999));
			writer.write("\r\n");
			

			if (label == 15)
				writer.write((point[0] - 432.3282882044739) + " " + (point[1] - 291.5621787855454) + " "
						+ (point[2] - 440.15589000771365) + " " + (point[3] - 315.30514932442173) + " "
						+ (length - 25));
			writer.write("\r\n");

			

			if (label == 16)
				writer.write((point[0] - 122.32450802270037 ) + " " + (point[1] - 321.3581919493962) + " "
						+ (point[2] - 137.31874376945635) + " " + (point[3] - 331.31675079948616) + " "
						+ (length - 17.99999999));
			writer.write("\r\n");
			if (label == 17)
				writer.write((point[0] - 404.4161512154085) + " " + (point[1] - 366.3849784182234) + " "
						+ (point[2] - 407.22494910278147 ) + " " + (point[3] - 340.53714186604935) + " "
						+ (length -26));
			writer.write("\r\n");
			if (label == 18)
				writer.write((point[0] - 249.99348275739862) + " " + (point[1] - 350.27744895704296) + " "
						+ (point[2] - 269.689292267354) + " " + (point[3] - 368.745667272505) + " "
						+ (length - 27));
			writer.write("\r\n");
			if (label == 19)
				writer.write((point[0] - 83.82342536303035) + " " + (point[1] - 509.29204677455294) + " "
						+ (point[2] - 93.71071227776532) + " " + (point[3] - 486.3303056065647) + " "
						+ (length -24.999));
			writer.write("\r\n");
			if (label == 20)
				writer.write((point[0] - 32.3923087857011) + " " + (point[1] - 543.8528582018996 ) + " "
						+ (point[2] - 59.607648991293395) + " " + (point[3] - 533.8366085297988) + " "
						+ (length - 28.99999));
			writer.write("\r\n");

			writer.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	private void Testerrorfour(final double[] point, int label, double length) {

		// Errorlist for Fake_big_file.tif
		try {
			FileWriter writer = new FileWriter("../res/error-Pnoise4secsnr15.txt", true);

			if (label == 1)
				writer.write((point[0] - 327.76209215919863) + " " + (point[1] - 47.380628285633605) + " "
						+ (point[2] - 329.4174200690935) + " " + (point[3] - 15.223204733509242) + " "
						+ (length - 32.199999999999996));
			writer.write("\r\n");
			if (label == 2)
				writer.write((point[0] - 4.875453301405078) + " " + (point[1] - 77.8324862348326) + " "
						+ (point[2] - 10.02527278408304) + " " + (point[3] - 48.88703133607818 ) + " "
						+ (length - 29.399999999999928));
			writer.write("\r\n");
			if (label == 3)
				writer.write((point[0] - 408.89938816919187) + " " + (point[1] - 69.82732619025738) + " "
						+ (point[2] - 410.0718187250601 ) + " " + (point[3] - 86.5863657388612) + " "
						+ (length - 16.80000000000002));
			writer.write("\r\n");
			if (label == 4)
				writer.write((point[0] - 447.7387612495996) + " " + (point[1] - 96.75498327965913) + " "
						+ (point[2] - 470.06999280189945) + " " + (point[3] - 98.50886393624926) + " "
						+ (length - 22.39999999999987));
			writer.write("\r\n");
			if (label == 5)
				writer.write((point[0] - 434.2810603432197) + " " + (point[1] - 152.33342406129356) + " "
						+ (point[2] - 455.9702936712432) + " " + (point[3] - 130.46524723645155) + " "
						+ (length - 30.800000000000065));
			writer.write("\r\n");
			if (label == 6)
				writer.write((point[0] - 345.10627168817774) + " " + (point[1] -152.1834118866404) + " "
						+ (point[2] - 368.8105444751357) + " " + (point[3] - 143.6307835396648 ) + " "
						+ (length - 25.200000000000138));
			writer.write("\r\n");

			if (label == 7)
				writer.write((point[0] - 3.908286620541236) + " " + (point[1] - 190.1178679247746) + " "
						+ (point[2] - 39.24433986335979) + " " + (point[3] - 173.14720673492457) + " "
						+ (length - 39.199999999999875));
			writer.write("\r\n");
			if (label == 8)
				writer.write((point[0] - 266.1954610919443) + " " + (point[1] - 189.91181955051073) + " "
						+ (point[2] - 282.6639724442005) + " " + (point[3] - 193.23268296844893 ) + " "
						+ (length - 16.79999999999975));
			writer.write("\r\n");
			if (label == 9)
				writer.write((point[0] - 96.04691366056015) + " " + (point[1] - 212.26565365789554) + " "
						+ (point[2] - 102.76104897210875) + " " + (point[3] - 193.85151813601408 ) + " "
						+ (length - 19.600000000000012));
			writer.write("\r\n");
			if (label == 10)
				writer.write((point[0] - 336.23861322738924) + " " + (point[1] - 220.056673149817) + " "
						+ (point[2] - 356.3403481351445) + " " + (point[3] - 241.51081619255652) + " "
						+ (length - 29.39999999999992));
			writer.write("\r\n");
			if (label == 11)
				writer.write((point[0] - 278.10070833534184) + " " + (point[1] - 239.64339233587503) + " "
						+ (point[2] - 290.62684490154123) + " " + (point[3] - 259.88038575695714 ) + " "
						+ (length - 23.799999999999958));
			writer.write("\r\n");
			if (label == 12)
				writer.write((point[0] - 300.7680991925678 ) + " " + (point[1] - 330.0621629416804) + " "
						+ (point[2] - 315.7781774303643) + " " + (point[3] - 340.35475985865105) + " "
						+ (length - 18.20000000000016));
			writer.write("\r\n");
			if (label == 14)
				writer.write((point[0] - 427.46299173976337) + " " + (point[1] - 359.22598602730744) + " "
						+ (point[2] - 440.0583785050841) + " " + (point[3] - 368.08692459113803) + " "
						+ (length - 15.399999999999876));
			writer.write("\r\n");

			if (label == 15)
				writer.write((point[0] - 362.8999580440501) + " " + (point[1] - 373.13074741573) + " "
						+ (point[2] - 374.8838153276987) + " " + (point[3] - 361.356743652985) + " "
						+ (length - 16.799999999999734));
			writer.write("\r\n");

			

			if (label == 16)
				writer.write((point[0] - 147.73224268059857 ) + " " + (point[1] - 403.88078930974433) + " "
						+ (point[2] - 178.19120908013963) + " " + (point[3] - 399.3100819387086) + " "
						+ (length - 30.800000000000185));
			writer.write("\r\n");
			
			if (label == 17)
				writer.write((point[0] - 281.66617970983623) + " " + (point[1] - 402.03989076001506) + " "
						+ (point[2] - 293.9073686549574) + " " + (point[3] - 419.1031053812364) + " "
						+ (length - 21.000000000000433));
			writer.write("\r\n");
			if (label == 18)
				writer.write((point[0] - 331.0027901598171) + " " + (point[1] - 442.34452697768046) + " "
						+ (point[2] - 346.26743336241316) + " " + (point[3] - 435.3280700443392) + " "
						+ (length - 16.799999999999844));
			writer.write("\r\n");
			if (label == 19)
				writer.write((point[0] - 14.375693702864483) + " " + (point[1] - 447.54202575273393) + " "
						+ (point[2] - 47.718389755498606) + " " + (point[3] - 468.1552404448105) + " "
						+ (length - 39.1999999999998));
			writer.write("\r\n");
			

			writer.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	private void Testerrorfive(final double[] point, int label, double length) {

		// Errorlist for Fake_big_file.tif
		try {
			FileWriter writer = new FileWriter("../res/error-Pnonoise5sec.txt", true);

			if (label == 1)
				writer.write((point[0] - 213.23931042861454) + " " + (point[1] - 28.768853068698554) + " "
						+ (point[2] - 231.51420771115022) + " " + (point[3] - 15.815547339550626) + " "
						+ (length - 22.40000000000014));
			writer.write("\r\n");
			if (label == 2)
				writer.write((point[0] - 262.5684393922506) + " " + (point[1] - 20.050913217219318) + " "
						+ (point[2] - 272.85194522198293) + " " + (point[3] - 39.95090085163391 ) + " "
						+ (length - 22.399999999999842));
			writer.write("\r\n");
			if (label == 3)
				writer.write((point[0] - 299.3825946207502) + " " + (point[1] - 94.80260795693071) + " "
						+ (point[2] - 331.31998473258363) + " " + (point[3] - 69.7354359416322) + " "
						+ (length - 40.60000000000002));
			writer.write("\r\n");
			if (label == 4)
				writer.write((point[0] - 45.675470197175386) + " " + (point[1] - 87.0860702029411) + " "
						+ (point[2] - 74.74090647661748 ) + " " + (point[3] - 76.89585853263253 ) + " "
						+ (length - 30.80000000000003));
			writer.write("\r\n");
			if (label == 5)
				writer.write((point[0] - 358.51322277498105) + " " + (point[1] - 78.69665600584534) + " "
						+ (point[2] - 379.8970891828737) + " " + (point[3] - 98.87313383150646 ) + " "
						+ (length - 29.399999999999903));
			writer.write("\r\n");
			if (label == 6)
				writer.write((point[0] - 49.76148246466852) + " " + (point[1] -124.09604469489028) + " "
						+ (point[2] - 85.1094880596964) + " " + (point[3] - 115.40820834434334 ) + " "
						+ (length - 36.400000000000105));
			writer.write("\r\n");

			if (label == 7)
				writer.write((point[0] - 223.4539445808136) + " " + (point[1] - 132.57748562595273) + " "
						+ (point[2] - 237.05277736463475 ) + " " + (point[3] - 116.57524369184648) + " "
						+ (length - 20.999999999999932));
			writer.write("\r\n");
			if (label == 8)
				writer.write((point[0] - 370.0214563609312) + " " + (point[1] - 142.97032407028482) + " "
						+ (point[2] - 398.6607651281875) + " " + (point[3] - 128.2519684103018) + " "
						+ (length - 32.19999999999996));
			writer.write("\r\n");
			if (label == 9)
				writer.write((point[0] - 327.4805366491584) + " " + (point[1] -131.2586576595218) + " "
						+ (point[2] - 345.49031411456883) + " " + (point[3] - 157.9511269279574 ) + " "
						+ (length - 32.19999999999976));
			writer.write("\r\n");
			if (label == 10)
				writer.write((point[0] - 153.68916577838593) + " " + (point[1] - 138.1918303951238) + " "
						+ (point[2] - 187.21475205548984) + " " + (point[3] - 140.4267945506451 ) + " "
						+ (length - 33.60000000000005));
			writer.write("\r\n");
			if (label == 11)
				writer.write((point[0] - 117.3400876334884) + " " + (point[1] - 199.9830506291008) + " "
						+ (point[2] - 125.90653218529553) + " " + (point[3] - 188.90983846004988  ) + " "
						+ (length - 14.000000000000126));
			writer.write("\r\n");
			if (label == 12)
				writer.write((point[0] - 355.9814722371587) + " " + (point[1] - 212.54284215973115) + " "
						+ (point[2] - 369.88580750919107) + " " + (point[3] - 189.86621408933854) + " "
						+ (length - 26.60000000000003));
			writer.write("\r\n");
			if (label == 13)
				writer.write((point[0] - 16.519680837702623) + " " + (point[1] - 208.35089834691382) + " "
						+ (point[2] - 44.81644292734521) + " " + (point[3] - 220.51367990551142 ) + " "
						+ (length - 30.79999999999993));
			writer.write("\r\n");
			if (label == 14)
				writer.write((point[0] - 322.0314627356028) + " " + (point[1] - 230.46786868249026) + " "
						+ (point[2] - 344.5515164445036) + " " + (point[3] - 211.56794326072625  ) + " "
						+ (length - 29.40000000000029));
			writer.write("\r\n");

			if (label == 15)
				writer.write((point[0] - 140.77986413780962) + " " + (point[1] - 222.04591247544067) + " "
						+ (point[2] - 169.16655005898042) + " " + (point[3] - 214.39371916810785) + " "
						+ (length - 29.399999999999753));
			writer.write("\r\n");

			

			if (label == 16)
				writer.write((point[0] - 150.2821689064564) + " " + (point[1] - 292.25054302689387) + " "
						+ (point[2] - 178.62332096989914) + " " + (point[3] - 269.4089340704993) + " "
						+ (length - 36.4000000000004));
			writer.write("\r\n");
			
			if (label == 17)
				writer.write((point[0] - 405.3441866021064) + " " + (point[1] -283.53349253317094) + " "
						+ (point[2] - 424.7711397927458) + " " + (point[3] - 305.6005152079386) + " "
						+ (length - 29.400000000000098));
			writer.write("\r\n");
			if (label == 18)
				writer.write((point[0] -359.5585400563441) + " " + (point[1] - 324.0239934624044) + " "
						+ (point[2] - 373.7287698808781) + " " + (point[3] - 337.56521152873336) + " "
						+ (length - 19.599999999999635));
			writer.write("\r\n");
			if (label == 19)
				writer.write((point[0] - 46.42913218018552) + " " + (point[1] - 388.9428935558981) + " "
						+ (point[2] - 61.00276863453459) + " " + (point[3] - 352.55265969437255) + " "
						+ (length - 39.20000000000066));
			writer.write("\r\n");
			if (label == 20)
				writer.write((point[0] - 181.55286568938038) + " " + (point[1] - 427.8287110953466) + " "
						+ (point[2] - 193.98127918991804) + " " + (point[3] - 414.53305248587316) + " "
						+ (length - 18.200000000000326));
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

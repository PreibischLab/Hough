package peakFitter;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
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
	final int maxiter = 2500;
	final double lambda = 1e-4;
	final double termepsilon = 1e-3;
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
		
		RandomAccessibleInterval<FloatType> currentimg = PerformWatershedding.CurrentLabelImage(intimg, inputimg, label);
		
		double newintercept = intercept;
		
		
		final long[] minCorner = PerformWatershedding.GetMincorners(intimg, label);
		final long[] maxCorner = PerformWatershedding.GetMaxcorners(intimg, label);
		FinalInterval smallinterval = new FinalInterval(minCorner, maxCorner);
		currentimg = Views.interval(currentimg, smallinterval);
		if(offsetting){
		
			currentimg = Views.offsetInterval(currentimg, smallinterval);
			
			newintercept = intercept - (smallinterval.realMin(1) - slope * smallinterval.realMin(0));
			
		}
		
		final Cursor<FloatType> outcursor = Views.iterable(currentimg).localizingCursor();
		
		final double maxintensityline = GetLocalmaxmin.computeMaxIntensity(currentimg);
		
		while (outcursor.hasNext()) {

			outcursor.fwd();
			
				if (outcursor.get().get()/maxintensityline  > 0.5) {
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

		MinandMax[2 * ndims] =   Math.min(psf[0], psf[1]);
		MinandMax[2 * ndims + 1] = 0.01 * maxintensityline;
		
		System.out.println("Label: " + label + " " + "Hough Detection: " + " StartX: " + MinandMax[0]  + " StartY: "
				+ MinandMax[1] + " EndX: " + MinandMax[2] + " EndY: " + MinandMax[3]);
		
		if (offsetting){
			System.out.println("Ofsetting on: " + "Label: " + label + " " + "Hough Detection: " 
		                                                      + " StartX: "
					+ (MinandMax[0] + smallinterval.realMin(0))   + " StartY: "
					+ (MinandMax[1] + smallinterval.realMin(1)) + " EndX: " 
		            + (MinandMax[2] + smallinterval.realMin(0)) + " EndY: " 
					+ (MinandMax[3] + smallinterval.realMin(1)));
			
		}
		for (int d = 0; d < ndims; ++d) {

			if (MinandMax[d] == Double.MAX_VALUE || MinandMax[d + ndims] == -Double.MIN_VALUE)
				return null;
			if (MinandMax[d] >= inputimg.dimension(d) || MinandMax[d + ndims] >= inputimg.dimension(d))
				return null;
			if (MinandMax[d] <= 0 || MinandMax[d + ndims] <= 0)
				return null;

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
			
			
			 RandomAccessibleInterval<FloatType> currentimg =  PerformWatershedding.CurrentLabelImage(intimg, inputimg, label);
			 
				final long[] minCorner = PerformWatershedding.GetMincorners(intimg, label);
				final long[] maxCorner = PerformWatershedding.GetMaxcorners(intimg, label);
				FinalInterval smallinterval = new FinalInterval(minCorner, maxCorner);
				currentimg = Views.interval(currentimg, smallinterval);
				if(offsetting){
				
					currentimg = Views.offsetInterval(currentimg, smallinterval);
					
				}
				
				final double[] fixed_param = new double[ndims + 1];

				for (int d = 0; d < ndims; ++d) {

					fixed_param[d] = 1.0 / Math.pow(psf[d], 2);
				}
				
				fixed_param[ndims] =  GetLocalmaxmin.computeMaxIntensity(currentimg);

			final double[] inistartpos = { start_param[0], start_param[1] };
			final double[] iniendpos = { start_param[2], start_param[3] };

			double inicutoffdistance = Distance(inistartpos, iniendpos);

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

				int iterations = 1800;

				double newslope = (endpos[1] - startpos[1]) / (endpos[0] - startpos[0]);
				double newintercept = (endpos[1]  -newslope *endpos[0]);
				double dx = finalparamstart[4]/ Math.sqrt( 1 + slope * slope) ;
				double dy = newslope * dx;
				double ds = finalparamstart[4];
				final double LMdist = sqDistance(startpos, endpos);
				
				double[] dxvector = {dx, dy};
				

				double[] startfit = new double[ndims];
				double[] endfit = new double[ndims];
				final double maxintensityline = fixed_param[ndims];

				final int numberofgaussians = 2; 
                
                
                
				startfit = peakFitter.GaussianMaskFit.sumofgaussianMaskFit(currentimg, intimg, startpos.clone(), psf, iterations,
						 dxvector, newslope,newintercept, maxintensityline, numberofgaussians, Endfit.Start, label);

				endfit = peakFitter.GaussianMaskFit.sumofgaussianMaskFit(currentimg, intimg, endpos.clone(), psf, iterations,
						 dxvector, newslope,newintercept, maxintensityline, numberofgaussians, Endfit.End, label);

				final double Maskdist = sqDistance(startfit, endfit);
				// If mask fits fail, return LM solver results, very crucial for
				// noisy data
				
				
				double[] returnparam = new double[2 * ndims + 2];
				double[] LMparam = new double[2 * ndims];
				
				if (offsetting){
					for (int d = 0; d < ndims; ++d) {
						startpos[d] += smallinterval.realMin(d);
						endpos[d] += smallinterval.realMin(d);
						startfit[d] += smallinterval.realMin(d);
						endfit[d] += smallinterval.realMin(d);
					}
					
					
				}
				System.out.println("ds: " + ds );
				System.out.println("LM solver : " + " StartX: " + startpos[0] + " StartY:  " + startpos[1]);
				System.out.println("LM solver : " + " EndX: " + endpos[0] + " EndY:  " + endpos[1]);
				System.out.println(" Length:  " + Math.sqrt(LMdist));
				
				for (int d = 0; d < ndims; ++d) {
					LMparam[d] = startpos[d];
					LMparam[ndims + d] = endpos[d];
				}
				
				for (int d = 0; d < ndims; ++d) {
						returnparam[d] = startfit[d];
						returnparam[ndims + d] = endfit[d];
				} 
				if (Math.abs(Math.sqrt(Maskdist) - Math.sqrt(LMdist)) > 10) {
                        System.out.println("Mask fits fail, returning LM solver results!");
						for (int d = 0; d < ndims; ++d) {
							returnparam[d] = startpos[d];
							returnparam[ndims + d] = endpos[d];
						}

					}
				
				//Testmovingline(returnparam, label, Distance(new double[] {returnparam[0], returnparam[1]},
				//		new double[] {returnparam[2], returnparam[3]}), 0);
				Testerrorfive(returnparam, label, Distance(new double[] {returnparam[0], returnparam[1]},
						new double[] {returnparam[2], returnparam[3]}));
				
				
			//	System.out.println("Number of gaussians for mask fit:" + (numberofgaussians) );

				
				returnparam[2* ndims] = finalparamstart[4];
				returnparam[2* ndims + 1] = finalparamstart[5];
				return returnparam;

			}

			else
				return null;

		}
	}
	public double[] Getfinaltrackparam(final double[] iniparam,
			final int label, final double[] psf, final int rate, boolean offsetting) throws Exception {

		
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
			
			
			final double[] finalparamstart = iniparam.clone();
			
			 RandomAccessibleInterval<FloatType> currentimg =  PerformWatershedding.CurrentLabelImage(intimg, inputimg, label);
			 
				final long[] minCorner = PerformWatershedding.GetMincorners(intimg, label);
				final long[] maxCorner = PerformWatershedding.GetMaxcorners(intimg, label);
				FinalInterval smallinterval = new FinalInterval(minCorner, maxCorner);
				currentimg = Views.interval(currentimg, smallinterval);
				if(offsetting){
				
					currentimg = Views.offsetInterval(currentimg, smallinterval);
					
					for (int d = 0; d < ndims; ++d){
						
						finalparamstart[d] -=smallinterval.realMin(d);
						finalparamstart[d + ndims] -= smallinterval.realMin(d);
						
					}
					
				}

				final double[] fixed_param = new double[ndims];

				for (int d = 0; d < ndims; ++d) {

					fixed_param[d] = 1.0 / Math.pow(psf[d], 2);
				}
				fixed_param[ndims] =  GetLocalmaxmin.computeMaxIntensity(currentimg);
			// LM solver part
						
				LevenbergMarquardtSolverLine.solve(X, finalparamstart, fixed_param, I, new GaussianLineds(), lambda,
						termepsilon, maxiter);

				final double[] startpos = { finalparamstart[0], finalparamstart[1] };
				final double[] endpos = { finalparamstart[2], finalparamstart[3] };
				// NaN protection: we prefer returning the crude estimate than
				// NaN
				for (int j = 0; j < finalparamstart.length; j++) {
					if (Double.isNaN(finalparamstart[j]))
						finalparamstart[j] = iniparam[j];
				}

				
				
				final double LMdist = sqDistance(startpos, endpos);
				
				double[] returnparam = new double[2 * ndims + 2];
				
				final double maxintensityline = fixed_param[ndims];
				
				
			
				int iterations = 1800;

				double newslope = (endpos[1] - startpos[1]) / (endpos[0] - startpos[0]);
				double newintercept = (endpos[1] - newslope *endpos[0]);
				double dx = finalparamstart[4]/ Math.sqrt( 1 + newslope * newslope);
				double dy = newslope * dx;
				double ds = finalparamstart[4];
				double[] dxvector = {dx, dy};
				System.out.println("Label: " +label);
				System.out.println("ds: " + ds );
				System.out.println("Initial solver : " + " StartX: " + startpos[0] + " StartY:  " + startpos[1]);
				System.out.println("Initial solver : " + " EndX: " + endpos[0] + " EndY:  " + endpos[1]);
				System.out.println(" Length:  " + Math.sqrt(LMdist));

				double[] startfit = new double[ndims];
				double[] endfit = new double[ndims];


				final int numberofgaussians = 2; 
                
                
                
				startfit = peakFitter.GaussianMaskFit.sumofgaussianMaskFit(currentimg, intimg, startpos.clone(), psf, iterations,
						 dxvector, newslope,newintercept, maxintensityline, numberofgaussians, Endfit.Start, label);

				endfit = peakFitter.GaussianMaskFit.sumofgaussianMaskFit(currentimg, intimg, endpos.clone(), psf, iterations,
						 dxvector, newslope,newintercept, maxintensityline, numberofgaussians, Endfit.End, label);

				final double Maskdist = sqDistance(startfit, endfit);
				// If mask fits fail, return LM solver results, very crucial for
				// noisy data
				
				
				
				
				if (offsetting){
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
				
				if (Math.abs(Math.sqrt(Maskdist) - Math.sqrt(LMdist)) > 10) {
                    System.out.println("Mask fits fail, returning LM solver results!");
					for (int d = 0; d < ndims; ++d) {
						returnparam[d] = startpos[d];
						returnparam[ndims + d] = endpos[d];
					}

				}
			Testmovingline(returnparam, label, Distance(new double[] {returnparam[0], returnparam[1]},
								new double[] {returnparam[2], returnparam[3]}), rate);
				
				
				returnparam[2* ndims] = finalparamstart[4];
				returnparam[2* ndims + 1] = finalparamstart[5];
				
				
				return returnparam;

			}

			

		
	}
	
	
	
	public int Getlabel(final Point linepoint ){
		
		
		RandomAccess<IntType> intranac = intimg.randomAccess();
		
		intranac.setPosition(linepoint);
		int currentlabel = intranac.get().get();
		
		return currentlabel;
	}

	
	
	private void Testerrorone(final double[] point, int label, double length) {

		// Errorlist for Fake_big_file.tif
		try {
			FileWriter writer = new FileWriter("error-Fake_one.txt", true);

			if (label == 1)
				writer.write((point[0] - 101.84356229708258) + " " + (point[1] - 34.72269000621467) + " "
						+ (point[2] - 113.01169863486068) + " " + (point[3] - 20.744421134813877) + " "
						+ (length - 17.891877204485056));
			writer.write("\r\n");
			if (label == 2)
				writer.write((point[0] - 383.70959453295944) + " " + (point[1] - 44.04000691401485) + " "
						+ (point[2] - 395.02257842025983 ) + " " + (point[3] - 30.02979900503253) + " "
						+ (length - 18.00748539044886));
			writer.write("\r\n");
			if (label == 3)
				writer.write((point[0] - 222.0517680211812) + " " + (point[1] - 54.27136909965083) + " "
						+ (point[2] - 233.22065561451288) + " " + (point[3] - 37.39675789813177) + " "
						+ (length - 20.23602118191496));
			writer.write("\r\n");
			if (label == 4)
				writer.write((point[0] - 264.86622286943935) + " " + (point[1] - 69.50041336082712) + " "
						+ (point[2] - 276.4634176659475) + " " + (point[3] - 55.14114238921126) + " "
						+ (length - 18.457616042827553));
			writer.write("\r\n");
			if (label == 5)
				writer.write((point[0] - 91.88527561650841) + " " + (point[1] - 89.36916692652001) + " "
						+ (point[2] - 103.54241296505866) + " " + (point[3] - 72.66055721968836) + " "
						+ (length - 20.37318064265357));
			writer.write("\r\n");
			if (label == 6)
				writer.write((point[0] - 78.72868102476444) + " " + (point[1] - 92.35128849766473) + " "
						+ (point[2] - 90.53310270073669) + " " + (point[3] - 77.39996918842077) + " "
						+ (length - 19.04957532836716));
			writer.write("\r\n");
       
			if (label == 7)
				writer.write((point[0] - 258.3028741382653) + " " + (point[1] - 169.0627765834802) + " "
						+ (point[2] - 270.22003738935496) + " " + (point[3] - 151.52609599285373 ) + " "
						+ (length - 21.20268723748897));
			writer.write("\r\n");
			if (label == 8)
				writer.write((point[0] - 300.5852204836343) + " " + (point[1] - 165.6067364109946) + " "
						+ (point[2] - 312.44765714375006) + " " + (point[3] - 154.44886993675695) + " "
						+ (length - 16.28543483521932));
			writer.write("\r\n");
			if (label == 9)
				writer.write((point[0] - 9.135689496597063) + " " + (point[1] - 220.6474589796404) + " "
						+ (point[2] - 20.269038650282194) + " " + (point[3] - 202.82555475178765) + " "
						+ (length - 21.01360829759175));
			writer.write("\r\n");
			if (label == 10)
				writer.write((point[0] - 322.347195774028) + " " + (point[1] - 251.36431346661595) + " "
						+ (point[2] - 334.23841952254406) + " " + (point[3] - 236.13478121016215) + " "
						+ (length - 19.322004424687385));
			writer.write("\r\n");
			if (label == 11)
				writer.write((point[0] - 182.5844087344585) + " " + (point[1] - 267.4427608509476) + " "
						+ (point[2] - 194.29386588004672) + " " + (point[3] - 250.758719045703) + " "
						+ (length - 20.383047799667164));
			writer.write("\r\n");
			if (label == 12)
				writer.write((point[0] - 43.13133334955034) + " " + (point[1] - 295.9839748756092) + " "
						+ (point[2] - 55.11228730824559) + " " + (point[3] -285.7189935156261) + " "
						+ (length - 15.776980068478753));
			writer.write("\r\n");
			if (label == 13)
				writer.write((point[0] - 87.90943527029843) + " " + (point[1] - 353.00924787494216) + " "
						+ (point[2] - 99.23449433610313) + " " + (point[3] - 334.4341802148607) + " "
						+ (length - 21.75523158738075));
			writer.write("\r\n");
        
			if (label == 14)
				writer.write((point[0] - 225.48251307414588) + " " + (point[1] - 366.9829333722683) + " "
						+ (point[2] - 236.9539633568867) + " " + (point[3] - 358.40376907956824) + " "
						+ (length - 14.32467212715649));
			writer.write("\r\n");

			if (label == 15)
				writer.write((point[0] - 35.62599947258589) + " " + (point[1] - 380.9522282309819) + " "
						+ (point[2] - 46.97940312410058) + " " + (point[3] - 371.014896388969) + " "
						+ (length - 15.08808598240713));
			writer.write("\r\n");
			
			if (label == 16)
				writer.write((point[0] - 318.79809608678465) + " " + (point[1] - 394.58485865054877) + " "
						+ (point[2] - 330.20918012825047) + " " + (point[3] - 382.22653982045375) + " "
						+ (length - 16.820846688192503));
			writer.write("\r\n");
			if (label == 17)
				writer.write((point[0] - 182.74199545178487) + " " + (point[1] - 399.2372805842788) + " "
						+ (point[2] - 193.745540933255) + " " + (point[3] - 386.7172540902593) + " "
						+ (length - 16.668205559499537));
			writer.write("\r\n");
			if (label == 18)
				writer.write((point[0] - 277.96411377430826) + " " + (point[1] - 414.01094489491777) + " "
						+ (point[2] - 289.02938093556145) + " " + (point[3] - 395.9148965060625) + " "
						+ (length - 21.211013757991367));
			writer.write("\r\n");
			if (label == 19)
				writer.write((point[0] - 256.14415851312833) + " " + (point[1] - 414.56338299609394 ) + " "
						+ (point[2] - 267.69268843741247) + " " + (point[3] - 402.5256330267631) + " "
						+ (length - 16.681605670204853));
			writer.write("\r\n");
			if (label == 20)
				writer.write((point[0] - 307.0221124229872) + " " + (point[1] - 425.74080251312057 ) + " "
						+ (point[2] - 318.90660598829754) + " " + (point[3] - 405.76253015063395 ) + " "
						+ (length - 23.2459147785928));
			writer.write("\r\n");
			
			writer.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private void Testerrortwo(final double[] point, int label, double length) {

		// Errorlist for Fake_big_file.tif
		try {
			FileWriter writer = new FileWriter("error-Fake_two.txt", true);

			if (label == 1)
				writer.write((point[0] - 338.3140365464876) + " " + (point[1] - 34.6062682241022) + " "
						+ (point[2] - 351.2873412307025) + " " + (point[3] - 64.56050672718627) + " "
						+ (length - 32.642963081330656));
			writer.write("\r\n");
			if (label == 2)
				writer.write((point[0] - 50.20240148337895) + " " + (point[1] - 74.70738084651842) + " "
						+ (point[2] - 62.899870601047446) + " " + (point[3] - 103.90557393663119 ) + " "
						+ (length - 31.83960115519122));
			writer.write("\r\n");
			if (label == 3)
				writer.write((point[0] - 3.328776397632022) + " " + (point[1] - 98.4960584135745) + " "
						+ (point[2] - 16.215625484255337) + " " + (point[3] - 125.7216261047921) + " "
						+ (length - 30.121461045083375));
			writer.write("\r\n");
			if (label == 4)
				writer.write((point[0] - 100.0364535301454) + " " + (point[1] - 112.27002417021258) + " "
						+ (point[2] - 112.66355873407907) + " " + (point[3] - 133.2351369102248) + " "
						+ (length - 24.474062556768793));
			writer.write("\r\n");
			if (label == 5)
				writer.write((point[0] - 316.8291502433898) + " " + (point[1] - 123.45464493609717) + " "
						+ (point[2] - 328.94284213631136) + " " + (point[3] - 149.3305924501906) + " "
						+ (length - 28.571072626500253));
			writer.write("\r\n");
			if (label == 6)
				writer.write((point[0] - 181.52203055786154) + " " + (point[1] - 136.73303593280056) + " "
						+ (point[2] - 194.30693230431007) + " " + (point[3] - 169.33350692240825) + " "
						+ (length - 35.01777293619049));
			writer.write("\r\n");
       
			if (label == 7)
				writer.write((point[0] - 220.3920284710936) + " " + (point[1] - 148.73466196724948) + " "
						+ (point[2] - 233.00881614719128) + " " + (point[3] - 170.67066765449704 ) + " "
						+ (length - 25.305566122390697));
			writer.write("\r\n");
			if (label == 8)
				writer.write((point[0] - 269.28105661270115) + " " + (point[1] - 154.26000633281953) + " "
						+ (point[2] - 281.92465988583575) + " " + (point[3] - 177.51878290536143) + " "
						+ (length - 26.4732202684117));
			writer.write("\r\n");
			if (label == 9)
				writer.write((point[0] - 313.6767573400796) + " " + (point[1] - 167.3048565453715) + " "
						+ (point[2] - 325.68360640963596) + " " + (point[3] - 190.86201905757824) + " "
						+ (length - 26.440581124582437));
			writer.write("\r\n");
			if (label == 10)
				writer.write((point[0] - 142.4113061517725) + " " + (point[1] - 197.12919949402738) + " "
						+ (point[2] - 155.26615799558937) + " " + (point[3] - 219.90907353258663) + " "
						+ (length - 26.156641166998263));
			writer.write("\r\n");
			if (label == 11)
				writer.write((point[0] - 181.0146684838801) + " " + (point[1] - 219.85082013938518) + " "
						+ (point[2] - 193.8977198168756) + " " + (point[3] - 246.95167976994097) + " "
						+ (length - 30.00715921848794));
			writer.write("\r\n");
			if (label == 12)
				writer.write((point[0] - 333.4436700399002) + " " + (point[1] - 226.89744562705945) + " "
						+ (point[2] - 346.18855535668945) + " " + (point[3] -250.2439212453266) + " "
						+ (length - 26.598684657938477));
			writer.write("\r\n");
			if (label == 13)
				writer.write((point[0] - 2.683434785065986) + " " + (point[1] - 227.39808836230733) + " "
						+ (point[2] - 14.895082826291112 ) + " " + (point[3] - 256.4683615614573) + " "
						+ (length - 31.531018565152223));
			writer.write("\r\n");
        
			if (label == 14)
				writer.write((point[0] - 139.40878293040046) + " " + (point[1] - 233.38981153080937) + " "
						+ (point[2] - 151.94475303337396) + " " + (point[3] - 263.7918526538246) + " "
						+ (length - 32.88517372415957));
			writer.write("\r\n");

			if (label == 15)
				writer.write((point[0] - 368.51992893964234) + " " + (point[1] - 246.49559242908282) + " "
						+ (point[2] - 380.8947083962499) + " " + (point[3] - 276.48119906526307) + " "
						+ (length - 32.43873875383279));
			writer.write("\r\n");
			
			if (label == 16)
				writer.write((point[0] - 322.6084356045666) + " " + (point[1] - 321.9972718564485) + " "
						+ (point[2] - 335.17857634720383) + " " + (point[3] - 353.9862176619938) + " "
						+ (length - 34.37006098394101));
			writer.write("\r\n");
			if (label == 17)
				writer.write((point[0] - 87.553951066869) + " " + (point[1] - 355.2475055449888) + " "
						+ (point[2] - 99.7736370435026) + " " + (point[3] - 376.7766785020526) + " "
						+ (length - 24.755322934324706));
			writer.write("\r\n");
			if (label == 18)
				writer.write((point[0] - 372.52471316888864) + " " + (point[1] - 359.3735279729181) + " "
						+ (point[2] - 384.7637643372141) + " " + (point[3] - 380.8555137213732) + " "
						+ (length - 24.723876823785485));
			writer.write("\r\n");
			if (label == 19)
				writer.write((point[0] - 345.9443907957086) + " " + (point[1] - 362.44175401675085 ) + " "
						+ (point[2] - 358.67385324684767) + " " + (point[3] - 387.1525254286868) + " "
						+ (length - 27.796788268933316));
			writer.write("\r\n");
			if (label == 20)
				writer.write((point[0] - 46.073101992551415) + " " + (point[1] - 406.08466530896965 ) + " "
						+ (point[2] - 58.58157266303186) + " " + (point[3] - 433.89196099959645 ) + " "
						+ (length - 30.491105787429575));
			writer.write("\r\n");
			
			writer.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private void Testerrorthree(final double[] point, int label, double length) {

		// Errorlist for Fake_big_file.tif
		try {
			FileWriter writer = new FileWriter("error-Fake_three.txt", true);

			if (label == 1)
				writer.write((point[0] - 102.61190201305033) + " " + (point[1] - 8.259538706043248) + " "
						+ (point[2] - 116.20495429979734) + " " + (point[3] - 47.49788822451655) + " "
						+ (length - 41.52612603415002));
			writer.write("\r\n");
			if (label == 2)
				writer.write((point[0] - 375.50314179291905) + " " + (point[1] - 14.912231448330996) + " "
						+ (point[2] - 387.8478932401985 ) + " " + (point[3] - 35.43179969630119) + " "
						+ (length - 23.94672356248794));
			writer.write("\r\n");
			if (label == 3)
				writer.write((point[0] - 166.2656617870702) + " " + (point[1] - 20.445405423494307) + " "
						+ (point[2] - 178.8210673581664) + " " + (point[3] - 57.85779380052519) + " "
						+ (length - 39.462957480763));
			writer.write("\r\n");
			if (label == 4)
				writer.write((point[0] - 386.5182277115242) + " " + (point[1] -20.701455890397558) + " "
						+ (point[2] - 399.21558816608757) + " " + (point[3] - 37.05412446917972) + " "
						+ (length - 20.703447349670988));
			writer.write("\r\n");
			if (label == 5)
				writer.write((point[0] - 104.48587735237139) + " " + (point[1] - 50.86749014591775) + " "
						+ (point[2] - 116.7486288669557) + " " + (point[3] - 87.72350007001484) + " "
						+ (length - 38.842509473945974));
			writer.write("\r\n");
			if (label == 6)
				writer.write((point[0] - 332.37937496471335) + " " + (point[1] - 81.67518017637789) + " "
						+ (point[2] - 345.7315293833449) + " " + (point[3] - 116.86750733938015) + " "
						+ (length - 37.640137071572525));
			writer.write("\r\n");
       
			if (label == 7)
				writer.write((point[0] - 64.54917058081423) + " " + (point[1] - 91.03967657791708) + " "
						+ (point[2] - 77.4242600784504) + " " + (point[3] - 108.99030738377144 ) + " "
						+ (length - 22.090565314183987));
			writer.write("\r\n");
			if (label == 8)
				writer.write((point[0] - 291.87059998038154) + " " + (point[1] - 91.7639867470957) + " "
						+ (point[2] - 304.30186580976414) + " " + (point[3] - 110.46284394495282) + " "
						+ (length - 22.45403372729783));
			writer.write("\r\n");
			if (label == 9)
				writer.write((point[0] - 150.13633156292963) + " " + (point[1] - 121.58435770257674) + " "
						+ (point[2] - 163.9645224360096) + " " + (point[3] - 154.54032141131216) + " "
						+ (length - 35.73953562644342));
			writer.write("\r\n");
			if (label == 10)
				writer.write((point[0] - 189.51595497473548) + " " + (point[1] - 144.65370360371278) + " "
						+ (point[2] - 201.73562989001596) + " " + (point[3] - 161.47891372002843) + " "
						+ (length - 20.79442594767418));
			writer.write("\r\n");
			if (label == 11)
				writer.write((point[0] - 40.7412533021802) + " " + (point[1] - 145.5782263378174) + " "
						+ (point[2] - 53.803573917493964) + " " + (point[3] - 160.40210181361812) + " "
						+ (length - 19.757821336860424));
			writer.write("\r\n");
			if (label == 12)
				writer.write((point[0] - 308.05934208691633) + " " + (point[1] - 179.0082922914368) + " "
						+ (point[2] - 321.2720644587525) + " " + (point[3] -203.26903587480268 ) + " "
						+ (length - 27.625345458347645));
			writer.write("\r\n");
			if (label == 13)
				writer.write((point[0] - 245.65286113426797) + " " + (point[1] - 180.99605440460334) + " "
						+ (point[2] - 259.12208398328875) + " " + (point[3] - 194.5249762317515) + " "
						+ (length - 19.09061785175298));
			writer.write("\r\n");
        
			if (label == 14)
				writer.write((point[0] - 287.22729116591654) + " " + (point[1] - 204.18458891357025) + " "
						+ (point[2] - 299.49883391642254) + " " + (point[3] - 238.21243491453822) + " "
						+ (length - 36.172988056602236));
			writer.write("\r\n");

			if (label == 15)
				writer.write((point[0] - 387.79159541310145) + " " + (point[1] - 220.4005654006134) + " "
						+ (point[2] - 399.8627095318793) + " " + (point[3] - 238.55108154821437) + " "
						+ (length - 21.798005241142636));
			writer.write("\r\n");
			
			if (label == 16)
				writer.write((point[0] - 305.3832535507701) + " " + (point[1] - 238.07718582965657) + " "
						+ (point[2] - 317.7703310239664) + " " + (point[3] - 253.47017090067237) + " "
						+ (length - 19.75812940851136));
			writer.write("\r\n");
			if (label == 17)
				writer.write((point[0] - 83.30628723439102) + " " + (point[1] - 271.7887197677545) + " "
						+ (point[2] - 97.2035177400693) + " " + (point[3] - 311.15070927231375) + " "
						+ (length - 41.74325374818051));
			writer.write("\r\n");
			if (label == 18)
				writer.write((point[0] - 196.8083178212894) + " " + (point[1] - 307.58702537538244) + " "
						+ (point[2] - 209.2698192157919) + " " + (point[3] - 330.102325095052) + " "
						+ (length - 25.73378593351043));
			writer.write("\r\n");
			if (label == 19)
				writer.write((point[0] - 187.97215142327738) + " " + (point[1] - 397.12932856571365 ) + " "
						+ (point[2] - 201.05541499056932) + " " + (point[3] - 415.7467803201413) + " "
						+ (length - 22.75480818200121));
			writer.write("\r\n");
			if (label == 20)
				writer.write((point[0] - 262.3361618230495) + " " + (point[1] - 397.1830771085031 ) + " "
						+ (point[2] - 274.34468615662405) + " " + (point[3] - 420.133780264653 ) + " "
						+ (length - 25.90249856735373));
			writer.write("\r\n");
			
			writer.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private void Testerrorfour(final double[] point, int label, double length) {

		// Errorlist for Fake_big_file.tif
		try {
			FileWriter writer = new FileWriter("error-Fake_four.txt", true);

			if (label == 1)
				writer.write((point[0] - 246.9961427022165) + " " + (point[1] - 19.804594837699458) + " "
						+ (point[2] - 261.36977964030376) + " " + (point[3] - 1.0418674544315567) + " "
						+ (length - 23.63559556446119));
			writer.write("\r\n");
			if (label == 2)
				writer.write((point[0] - 372.56166620886665) + " " + (point[1] - 21.720349759718815) + " "
						+ (point[2] - 387.25534099606835 ) + " " + (point[3] - 13.708919626473548) + " "
						+ (length - 16.735802685617312));
			writer.write("\r\n");
			if (label == 3)
				writer.write((point[0] - 128.50544705924705) + " " + (point[1] - 61.96424673453504) + " "
						+ (point[2] - 142.0507088984846) + " " + (point[3] - 44.49777923994937) + " "
						+ (length - 22.103203501595033));
			writer.write("\r\n");
			if (label == 4)
				writer.write((point[0] - 196.31944223919538) + " " + (point[1] -69.59397621060532) + " "
						+ (point[2] - 210.9692701132741) + " " + (point[3] - 45.69409751268013 ) + " "
						+ (length - 28.032510742273367));
			writer.write("\r\n");
			if (label == 5)
				writer.write((point[0] - 58.61112320382826) + " " + (point[1] - 63.42905477642745) + " "
						+ (point[2] - 71.66576978169842) + " " + (point[3] - 51.00748898087748) + " "
						+ (length - 18.01996376484353));
			writer.write("\r\n");
			if (label == 6)
				writer.write((point[0] - 91.67834782432345) + " " + (point[1] - 62.678064701487045) + " "
						+ (point[2] - 104.36211216953461) + " " + (point[3] - 53.8672435624161) + " "
						+ (length - 15.443718694328426));
			writer.write("\r\n");
       
			if (label == 7)
				writer.write((point[0] - 277.7761745378519) + " " + (point[1] - 93.02698136834665) + " "
						+ (point[2] - 291.64732993963435) + " " + (point[3] - 87.76837114864867 ) + " "
						+ (length - 14.834484609284909));
			writer.write("\r\n");
			if (label == 8)
				writer.write((point[0] - 182.9069448917257) + " " + (point[1] - 103.61775776279539) + " "
						+ (point[2] - 196.05993089823494) + " " + (point[3] - 94.13698786731403) + " "
						+ (length - 16.213760782079333));
			writer.write("\r\n");
			if (label == 9)
				writer.write((point[0] - 255.4078812743308) + " " + (point[1] - 169.86605171478112) + " "
						+ (point[2] - 269.5747729176865) + " " + (point[3] - 171.21638730628098) + " "
						+ (length - 14.231100626594317));
			writer.write("\r\n");
			if (label == 10)
				writer.write((point[0] - 181.72534344407936) + " " + (point[1] - 185.81542903477995) + " "
						+ (point[2] - 194.43942497020035) + " " + (point[3] - 171.47225665119402) + " "
						+ (length - 19.167014975684488));
			writer.write("\r\n");
			if (label == 11)
				writer.write((point[0] - 51.82675589676167) + " " + (point[1] - 180.1851579327744) + " "
						+ (point[2] - 64.3276628519142) + " " + (point[3] - 178.97461153417208) + " "
						+ (length - 12.559382838521563));
			writer.write("\r\n");
			if (label == 12)
				writer.write((point[0] - 351.0942761263851) + " " + (point[1] - 214.98359386963045) + " "
						+ (point[2] - 364.46938077544723) + " " + (point[3] -192.65919780274118 ) + " "
						+ (length - 26.02445165848273));
			writer.write("\r\n");
			if (label == 13)
				writer.write((point[0] - 277.53462196685587) + " " + (point[1] - 299.24612843893607) + " "
						+ (point[2] - 292.1273078410969) + " " + (point[3] - 295.3689064514161) + " "
						+ (length - 15.098984448127062));
			writer.write("\r\n");
        
			if (label == 14)
				writer.write((point[0] - 226.92311949430476) + " " + (point[1] - 344.7937059034321) + " "
						+ (point[2] - 241.3425032845104 ) + " " + (point[3] - 313.98844907831386) + " "
						+ (length - 34.01297512642403));
			writer.write("\r\n");

			if (label == 15)
				writer.write((point[0] - 89.01673780411032) + " " + (point[1] - 334.288777809216) + " "
						+ (point[2] - 102.59180581181907) + " " + (point[3] - 327.3442287581981) + " "
						+ (length - 15.248253438866751));
			writer.write("\r\n");
			
			if (label == 16)
				writer.write((point[0] - 380.13035754507115) + " " + (point[1] - 364.8437584545042) + " "
						+ (point[2] - 393.67672482889003) + " " + (point[3] - 339.75576791755964) + " "
						+ (length - 28.511600021218364));
			writer.write("\r\n");
			if (label == 17)
				writer.write((point[0] -79.08785656321496) + " " + (point[1] - 379.1919061648394) + " "
						+ (point[2] -93.38053218888024) + " " + (point[3] - 367.83353345808314) + " "
						+ (length - 18.256319647893765));
			writer.write("\r\n");
			if (label == 18)
				writer.write((point[0] - 28.394015551227405) + " " + (point[1] - 410.32826030425895) + " "
						+ (point[2] - 42.87688838030825 ) + " " + (point[3] - 386.760186564061) + " "
						+ (length - 27.662387915845905));
			writer.write("\r\n");
			if (label == 19)
				writer.write((point[0] - 107.48363735749044) + " " + (point[1] - 422.44065142624805 ) + " "
						+ (point[2] - 121.08890896190134) + " " + (point[3] - 393.68727069233506) + " "
						+ (length - 31.809751949035192));
			writer.write("\r\n");
			if (label == 20)
				writer.write((point[0] - 364.52057197326985) + " " + (point[1] - 425.57696520366466 ) + " "
						+ (point[2] - 378.4745312914248) + " " + (point[3] - 398.09246755793447 ) + " "
						+ (length - 30.823863993517612));
			writer.write("\r\n");
			
			writer.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	
	private void Testerrorfive(final double[] point, int label, double length) {

		// Errorlist for Fake_big_file.tif
		try {
			FileWriter writer = new FileWriter("error-Fake_five.txt", true);

			if (label == 1)
				writer.write((point[0] - 222.90371643469152) + " " + (point[1] - 16.028490341799035) + " "
						+ (point[2] - 255.2312724055033) + " " + (point[3] - 36.02343425047382) + " "
						+ (length - 38.011427978399055));
			writer.write("\r\n");
			if (label == 2)
				writer.write((point[0] - 360.91877489118144) + " " + (point[1] - 17.45868075904963) + " "
						+ (point[2] - 391.7462103358474 ) + " " + (point[3] - 63.491803165883596) + " "
						+ (length - 55.401977713955375));
			writer.write("\r\n");
			if (label == 3)
				writer.write((point[0] - 128.96895079354505) + " " + (point[1] - 55.02041872625349) + " "
						+ (point[2] - 159.71930546208145) + " " + (point[3] - 59.17680554453315) + " "
						+ (length - 31.02998330041361));
			writer.write("\r\n");
			if (label == 4)
				writer.write((point[0] - 198.8351823841789) + " " + (point[1] -136.2169214947499) + " "
						+ (point[2] - 230.26320248688194) + " " + (point[3] - 186.4924501678817) + " "
						+ (length - 59.29038059364104));
			writer.write("\r\n");
			if (label == 5)
				writer.write((point[0] - 106.15907400122536 ) + " " + (point[1] - 151.76688612246622) + " "
						+ (point[2] - 136.3486806598561) + " " + (point[3] - 201.56365209669153) + " "
						+ (length - 58.23341181568021));
			writer.write("\r\n");
			if (label == 6)
				writer.write((point[0] - 258.23494092128317) + " " + (point[1] - 202.2428912845651) + " "
						+ (point[2] - 289.6033517149683) + " " + (point[3] - 239.63541831163295 ) + " "
						+ (length - 48.80756368834017));
			writer.write("\r\n");
       
			if (label == 7)
				writer.write((point[0] - 161.12270856208303) + " " + (point[1] - 219.44944752857754) + " "
						+ (point[2] - 192.30906160583731) + " " + (point[3] - 232.04958503613773 ) + " "
						+ (length -33.63557761328185));
			writer.write("\r\n");
			if (label == 8)
				writer.write((point[0] - 388.78833042563593) + " " + (point[1] - 351.8343450878364) + " "
						+ (point[2] - 418.9224848972848) + " " + (point[3] - 375.98872415617024 ) + " "
						+ (length - 38.61995978633273));
			writer.write("\r\n");
			if (label == 9)
				writer.write((point[0] - 24.202688921682622) + " " + (point[1] - 362.7494073097436) + " "
						+ (point[2] - 54.487901283178786) + " " + (point[3] - 405.33575976290547) + " "
						+ (length - 52.25697564006019));
			writer.write("\r\n");
			if (label == 10)
				writer.write((point[0] - 314.24729092693474) + " " + (point[1] - 379.8908139375993) + " "
						+ (point[2] - 344.34962327362956) + " " + (point[3] - 415.1772672850986) + " "
						+ (length - 46.38193832254224));
			writer.write("\r\n");
			
			
			writer.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private PointSampleList<FloatType> gatherfullData(final int label, final boolean offsetting) {
		final PointSampleList<FloatType> datalist = new PointSampleList<FloatType>(ndims);

		RandomAccessibleInterval<FloatType> currentimg = PerformWatershedding.CurrentLabelImage(intimg, inputimg, label);
		
		

		boolean outofbounds = false;
		final long[] minCorner = PerformWatershedding.GetMincorners(intimg, label);
		final long[] maxCorner = PerformWatershedding.GetMaxcorners(intimg, label);
		FinalInterval smallinterval = new FinalInterval(minCorner, maxCorner);
		
		if(offsetting){
		
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
	
	private void Testmovingline (final double[] point, int label, double length, int rate){
		try {
			FileWriter writer = new FileWriter("error-Fake_movingline.txt", true);

		if (rate == 0){
			
			
			if (label == 1)
				writer.write((point[0] - 169.73720446120004) + " " + (point[1] - 104.95635903304128 ) + " "
						+ (point[2] - 181.43458804585626) + " " + (point[3] - 60.78333591544249) + " "
						+ (length - 45.69556602203874));
			writer.write("\r\n");
			if (label == 2)
				writer.write((point[0] - 294.6604728889232) + " " + (point[1] - 109.07711973280686) + " "
						+ (point[2] - 305.2140902032087 ) + " " + (point[3] - 82.10337364262949) + " "
						+ (length - 28.964837588941272));
			writer.write("\r\n");
			if (label == 3)
				writer.write((point[0] - 391.2990529717664) + " " + (point[1] - 124.0347325138523) + " "
						+ (point[2] - 401.71221291284917) + " " + (point[3] - 82.7441754252086) + " "
						+ (length - 42.58337709305259));
			writer.write("\r\n");
			if (label == 4)
				writer.write((point[0] - 203.44041398569848) + " " + (point[1] -130.38479493829118) + " "
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
				writer.write((point[0] - 360.6948730664424 ) + " " + (point[1] - 246.8371168164922) + " "
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
						+ (length -39.35149161408799));
			writer.write("\r\n");
			if (label == 10)
				writer.write((point[0] - 181.75417277692245) + " " + (point[1] - 377.86831136750124) + " "
						+ (point[2] -193.23751455228975) + " " + (point[3] - 351.97062628010985) + " "
						+ (length - 28.329441067828885));
			writer.write("\r\n");
		}
		
         if (rate == 1){
			
			
			if (label == 1)
				writer.write((point[0] - 171.1859449967859) + " " + (point[1] - 99.48545618483672) + " "
						+ (point[2] - 181.43458804585626) + " " + (point[3] - 60.78333591544249) + " "
						+ (length - 40.036093686746135));
			writer.write("\r\n");
			if (label == 2)
				writer.write((point[0] - 298.63632549289247) + " " + (point[1] - 98.91533013614685) + " "
						+ (point[2] - 305.2140902032087 ) + " " + (point[3] - 82.10337364262949) + " "
						+ (length - 18.05294628929588));
			writer.write("\r\n");
			if (label == 3)
				writer.write((point[0] -392.2658523406384) + " " + (point[1] - 120.20115224047345) + " "
						+ (point[2] - 401.71221291284917) + " " + (point[3] - 82.7441754252086) + " "
						+ (length - 38.62976624572697));
			writer.write("\r\n");
			if (label == 4)
				writer.write((point[0] - 206.0628373691796) + " " + (point[1] - 123.40091820460714) + " "
						+ (point[2] - 214.9881816669676) + " " + (point[3] - 99.63148998257181 ) + " "
						+ (length - 25.38990919315283));
			writer.write("\r\n");
			if (label == 5)
				writer.write((point[0] -135.10455586538396 ) + " " + (point[1] - 209.6406004482957) + " "
						+ (point[2] - 143.00136725366892) + " " + (point[3] - 177.4220897013742) + " "
						+ (length - 33.17215797700901));
			writer.write("\r\n");
			if (label == 6)
				writer.write((point[0] - 4.798655006524217 ) + " " + (point[1] - 227.44491761331008) + " "
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
		
		
         if (rate == 2){
 			
 			
 			if (label == 1)
 				writer.write((point[0] -172.43796447094846) + " " + (point[1] - 94.75743433485547) + " "
 						+ (point[2] - 181.43458804585626) + " " + (point[3] - 60.78333591544249) + " "
 						+ (length - 35.14510775571112));
 			writer.write("\r\n");
 			if (label == 2)
 				writer.write((point[0] - 298.63632549289247) + " " + (point[1] - 98.91533013614685) + " "
 						+ (point[2] - 305.2140902032087 ) + " " + (point[3] - 82.10337364262949) + " "
 						+ (length - 18.05294628929588));
 			writer.write("\r\n");
 			if (label == 3)
 				writer.write((point[0] -393.1753225631553) + " " + (point[1] - 116.59489512392793) + " "
 						+ (point[2] - 401.71221291284917) + " " + (point[3] - 82.7441754252086) + " "
 						+ (length -34.91059611298497));
 			writer.write("\r\n");
 			if (label == 4)
 				writer.write((point[0] - 207.40081502209952) + " " + (point[1] - 119.83769837546845) + " "
 						+ (point[2] - 214.9881816669676) + " " + (point[3] - 99.63148998257181 ) + " "
 						+ (length - 21.583766821869236));
 			writer.write("\r\n");
 			if (label == 5)
 				writer.write((point[0] - 136.04277722443922) + " " + (point[1] - 205.81271429863318) + " "
 						+ (point[2] - 143.00136725366892) + " " + (point[3] - 177.4220897013742) + " "
 						+ (length - 29.23096885526344));
 			writer.write("\r\n");
 			if (label == 6)
 				writer.write((point[0] - 5.034495309332331 ) + " " + (point[1] - 226.8034846994553) + " "
 						+ (point[2] - 15.44169613220083) + " " + (point[3] - 198.4982244474968) + " "
 						+ (length - 30.157877692215305));
 			writer.write("\r\n");
        
 			if (label == 7)
 				writer.write((point[0] - 364.68594610462594) + " " + (point[1] - 233.3041243863031 ) + " "
 						+ (point[2] - 370.8155639541149) + " " + (point[3] - 212.5197210031867) + " "
 						+ (length - 21.669417135051745));
 			writer.write("\r\n");
 			if (label == 8)
 				writer.write((point[0] - 307.7174271751884) + " " + (point[1] -328.81123982137103) + " "
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

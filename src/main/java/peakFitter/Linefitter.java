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
		if(offsetting){
		
			currentimg = Views.offsetInterval(currentimg, smallinterval);
			
			newintercept = intercept - (smallinterval.realMin(1) - slope * smallinterval.realMin(0));
			
		}
		
		final Cursor<FloatType> outcursor = Views.iterable(currentimg).localizingCursor();
		
		final double maxintensityline = GetLocalmaxmin.computeMaxIntensityinlabel(inputimg, intimg, label);
		
		
		while (outcursor.hasNext()) {

			outcursor.fwd();
			
				if (outcursor.get().get() / maxintensityline >= 0.5) {
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

		MinandMax[2 * ndims] = 1.0;
		MinandMax[2 * ndims + 1] = 1.0;
		MinandMax[2 * ndims + 2] = 0;
		
		System.out.println("Label: " + label + " " + "Hough Detection: " + " StartX: " + MinandMax[0]  + " StartY: "
				+ MinandMax[1] + " EndX: " + MinandMax[2] + " EndY: " + MinandMax[3]);
		
		if (offsetting){
			System.out.println("Ofsetting on: " + "Label: " + label + " " + "Hough Detection: " 
		                                                      + " StartX: "
					+ (MinandMax[0] + smallinterval.realMin(0))   + " StartY: "
					+ (MinandMax[1] + smallinterval.realMin(1)) + " EndX: " 
		            + (MinandMax[2] + smallinterval.realMin(0)) + " EndY: " 
					+ (MinandMax[3] + +smallinterval.realMin(1)));
			
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
			final double[] fixed_param = new double[ndims + 1];

			for (int d = 0; d < ndims; ++d) {

				fixed_param[d] = 1.0 / Math.pow(psf[d], 2);
			}
			final double[] finalparamstart = start_param.clone();
			// LM solver part
			int maxiter = 1500;
			double lambda = 1e-3;
			double termepsilon = 1e-2;
			
			 RandomAccessibleInterval<FloatType> currentimg =  PerformWatershedding.CurrentLabelImage(intimg, inputimg, label);
			 
				final long[] minCorner = PerformWatershedding.GetMincorners(intimg, label);
				final long[] maxCorner = PerformWatershedding.GetMaxcorners(intimg, label);
				FinalInterval smallinterval = new FinalInterval(minCorner, maxCorner);
				
				if(offsetting){
				
					currentimg = Views.offsetInterval(currentimg, smallinterval);
					
				}

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

				int iterations = 800;

				double newslope = (endpos[1] - startpos[1]) / (endpos[0] - startpos[0]);
				double newintercept = (endpos[1] +(1/ newslope) *endpos[0]);
				double dx = finalparamstart[4] / Math.sqrt(1 + newslope * newslope);
				double dy = newslope * dx;
				double ds = finalparamstart[4];
				final double LMdist = sqDistance(startpos, endpos);
				/*
				double[] dxvector = {dx, dy};
				

				double[] startfit = new double[ndims];
				double[] endfit = new double[ndims];


				final int numberofgaussians = 2; 
                final double lineintensity = finalparamstart[5];
                
                
				startfit = peakFitter.GaussianMaskFit.sumofgaussianMaskFit(currentimg, intimg, startpos.clone(), psf, iterations,
						lineintensity, dxvector, newslope,newintercept, numberofgaussians, Endfit.Start, label);

				endfit = peakFitter.GaussianMaskFit.sumofgaussianMaskFit(currentimg, intimg, endpos.clone(), psf, iterations,
						lineintensity, dxvector, newslope,newintercept, numberofgaussians, Endfit.End, label);

				final double Maskdist = sqDistance(startfit, endfit);
				// If mask fits fail, return LM solver results, very crucial for
				// noisy data
				
				*/
				double[] returnparam = new double[2 * ndims];
				double[] LMparam = new double[2 * ndims];
				
				if (offsetting){
					for (int d = 0; d < ndims; ++d) {
						startpos[d] += smallinterval.realMin(d);
						endpos[d] += smallinterval.realMin(d);
					//	startfit[d] += smallinterval.realMin(d);
					//	endfit[d] += smallinterval.realMin(d);
					}
					
					
				}
				System.out.println("ds: " + ds + " " + "Const:" + finalparamstart[6]);
				System.out.println("LM solver : " + " StartX: " + startpos[0] + " StartY:  " + startpos[1]);
				System.out.println("LM solver : " + " EndX: " + endpos[0] + " EndY:  " + endpos[1]);
				System.out.println(" Length:  " + Math.sqrt(LMdist));
				
				for (int d = 0; d < ndims; ++d) {
					LMparam[d] = startpos[d];
					LMparam[ndims + d] = endpos[d];
				}
				/*
				for (int d = 0; d < ndims; ++d) {
						returnparam[d] = startfit[d];
						returnparam[ndims + d] = endfit[d];
				} */
			//		if (Math.abs(Math.sqrt(Maskdist) - Math.sqrt(LMdist)) > 10) {
             //           System.out.println("Mask fits fail, returning LM solver results!");
						for (int d = 0; d < ndims; ++d) {
							returnparam[d] = startpos[d];
							returnparam[ndims + d] = endpos[d];
						}

				//	}
				
				Testerroreight(LMparam, label, Distance(new double[] {LMparam[0], LMparam[1]},
						new double[] {LMparam[2], LMparam[3]}));
				
			//	System.out.println("Number of gaussians for mask fit:" + (numberofgaussians) );

				return returnparam;

			}

			else
				return null;

		}
	}

	private void Testerror(final double[] point, int label, double length) {

		// Errorlist for Fake_big_file.tif
		try {
			FileWriter writer = new FileWriter("error-Fake_big_file-offsettingLM.txt", true);

			if (label == 1)
				writer.write((point[0] - 146.94592464242646) + " " + (point[1] - 55.42312239408533) + " "
						+ (point[2] - 153.54928295839164) + " " + (point[3] - 49.42739020885564) + " "
						+ (length - 8.919257003025853));
			writer.write("\r\n");
			if (label == 2)
				writer.write((point[0] - 229.8239785496787) + " " + (point[1] - 64.27520704680428) + " "
						+ (point[2] - 234.1218573487757 ) + " " + (point[3] - 54.69958671952257) + " "
						+ (length - 10.495916673829786));
			writer.write("\r\n");
			if (label == 3)
				writer.write((point[0] - 302.9465004497457) + " " + (point[1] - 58.47639954792287) + " "
						+ (point[2] - 313.16619975599116) + " " + (point[3] - 73.52451716171643) + " "
						+ (length - 18.19032978339419));
			writer.write("\r\n");
			if (label == 4)
				writer.write((point[0] - 240.1193313521409) + " " + (point[1] - 116.67602295228242) + " "
						+ (point[2] - 245.5693164385891) + " " + (point[3] - 107.23807809464897) + " "
						+ (length - 10.898492582840879));
			writer.write("\r\n");
			if (label == 5)
				writer.write((point[0] - 133.35913522709873) + " " + (point[1] - 140.7563693740396) + " "
						+ (point[2] - 136.5441108511015) + " " + (point[3] - 155.84081987255448) + " "
						+ (length - 15.417026839427818));
			writer.write("\r\n");
		
       
			if (label == 7)
				writer.write((point[0] - 215.28576496031062) + " " + (point[1] - 196.00276945260316) + " "
						+ (point[2] - 218.68751554101746) + " " + (point[3] - 211.05663130841134) + " "
						+ (length - 15.433426832272081));
			writer.write("\r\n");
			if (label == 8)
				writer.write((point[0] - 121.13389234638807) + " " + (point[1] - 196.582277391324) + " "
						+ (point[2] - 130.5276503416368) + " " + (point[3] - 214.46479753020282) + " "
						+ (length - 20.19968355669726));
			writer.write("\r\n");
			if (label == 9)
				writer.write((point[0] - 340.37580237095966) + " " + (point[1] - 219.09631742654946) + " "
						+ (point[2] - 351.0314266781306) + " " + (point[3] - 234.43943863070862) + " "
						+ (length - 18.68030239747446));
			writer.write("\r\n");
			if (label == 10)
				writer.write((point[0] - 385.5135453587466) + " " + (point[1] - 230.40541983841212) + " "
						+ (point[2] - 392.6022931853405) + " " + (point[3] - 245.46467445749235) + " "
						+ (length - 16.64426313873128));
			writer.write("\r\n");
			if (label == 11)
				writer.write((point[0] - 5.556310385317046) + " " + (point[1] - 272.71517939626614) + " "
						+ (point[2] - 12.726024523817513) + " " + (point[3] - 263.9600050091784) + " "
						+ (length - 11.316266141096644));
			writer.write("\r\n");
			if (label == 12)
				writer.write((point[0] - 308.05785927861723) + " " + (point[1] - 275.08437728781024) + " "
						+ (point[2] - 314.6860883859905) + " " + (point[3] - 292.793305456796) + " "
						+ (length - 18.908716455490264));
			writer.write("\r\n");
			if (label == 13)
				writer.write((point[0] - 60.969008493931426) + " " + (point[1] - 282.3347616677545) + " "
						+ (point[2] - 62.4978725178862) + " " + (point[3] - 294.59885719136463) + " "
						+ (length - 12.359023594765786));
			writer.write("\r\n");
        
			if (label == 14)
				writer.write((point[0] - 79.27699219076632) + " " + (point[1] - 294.24744896193636) + " "
						+ (point[2] - 85.54813185659465) + " " + (point[3] - 283.7722508049237) + " "
						+ (length - 12.208888939498404));
			writer.write("\r\n");

			
			
			if (label == 16)
				writer.write((point[0] - 174.4791492931807) + " " + (point[1] - 300.8682377344869) + " "
						+ (point[2] - 184.98436572151002) + " " + (point[3] - 322.86506735085936) + " "
						+ (length - 24.376629901972915));
			writer.write("\r\n");
			if (label == 17)
				writer.write((point[0] - 387.5044473398028) + " " + (point[1] - 353.3233494741606) + " "
						+ (point[2] - 392.26957599273436) + " " + (point[3] - 345.6137121862541) + " "
						+ (length - 9.063385581010413));
			writer.write("\r\n");
			if (label == 18)
				writer.write((point[0] - 147.68645636726484) + " " + (point[1] - 355.0991002838538) + " "
						+ (point[2] - 154.95312760941732) + " " + (point[3] - 368.27466301377046) + " "
						+ (length - 15.0465931091226));
			writer.write("\r\n");
			if (label == 19)
				writer.write((point[0] - 31.239230737346844) + " " + (point[1] - 362.97213316931663 ) + " "
						+ (point[2] - 34.14047219898989) + " " + (point[3] - 377.24877146243585) + " "
						+ (length - 14.568445454862875));
			writer.write("\r\n");
			if (label == 20)
				writer.write((point[0] - 97.97120006079818) + " " + (point[1] - 369.5266828584096 ) + " "
						+ (point[2] - 105.10103519287271) + " " + (point[3] - 387.0954284840828) + " "
						+ (length - 18.960363178752036));
			writer.write("\r\n");
			if (label == 21)
				writer.write((point[0] - 284.6932383894) + " " + (point[1] - 375.30728038281813 ) + " "
						+ (point[2] - 291.6624758784747 ) + " " + (point[3] - 394.13153147321896) + " "
						+ (length - 20.07293452122985));
			if (label == 22)
				writer.write((point[0] - 304.75019499065377 ) + " " + (point[1] - 387.85799639315786 ) + " "
						+ (point[2] - 306.99650103631984) + " " + (point[3] - 377.82242936768705) + " "
						+ (length - 10.283895000120967));
			writer.write("\r\n");
			if (label == 23)
				writer.write((point[0] - 110.8275331566493 ) + " " + (point[1] - 376.97411872343633 ) + " "
						+ (point[2] - 116.2795632060759  ) + " " + (point[3] - 391.63032118810133 ) + " "
						+ (length - 15.637419938887078));
			writer.write("\r\n");
			if (label == 24)
				writer.write((point[0] - 39.57235449465267 ) + " " + (point[1] - 383.13064284394886 ) + " "
						+ (point[2] - 50.77250286652348 ) + " " + (point[3] - 404.26314211550664  ) + " "
						+ (length - 23.917061880889765));
			writer.write("\r\n");
			if (label == 25)
				writer.write((point[0] - 56.282438519576644 ) + " " + (point[1] - 395.16307287277664 ) + " "
						+ (point[2] - 62.250265839645685   ) + " " + (point[3] - 411.8090857339708 ) + " "
						+ (length - 17.683458572835903));
			writer.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	private void Testerroreight(final double[] point, int label, double length) {

		// Errorlist for Fake_big_file.tif
		try {
			FileWriter writer = new FileWriter("error-Fake_eight-offsettingLM.txt", true);

			if (label == 1)
				writer.write((point[0] - 103.42410325530419) + " " + (point[1] - 20.066334908307294) + " "
						+ (point[2] - 113.01169863486068) + " " + (point[3] - 20.744421134813877) + " "
						+ (length - 9.61154441766103));
			writer.write("\r\n");
			if (label == 2)
				writer.write((point[0] - 385.0149251470887) + " " + (point[1] - 28.8351210942548) + " "
						+ (point[2] - 395.02257842025983 ) + " " + (point[3] - 30.02979900503253) + " "
						+ (length - 10.078709210336052));
			writer.write("\r\n");
			if (label == 3)
				writer.write((point[0] - 223.63088159385103) + " " + (point[1] - 46.662735533349625) + " "
						+ (point[2] - 233.22065561451288) + " " + (point[3] - 37.39675789813177) + " "
						+ (length - 13.334995586940344));
			writer.write("\r\n");
			if (label == 4)
				writer.write((point[0] - 265.63155275607386) + " " + (point[1] - 54.06721342595662) + " "
						+ (point[2] - 276.4634176659475) + " " + (point[3] - 55.14114238921126) + " "
						+ (length - 10.884972248190095));
			writer.write("\r\n");
			if (label == 6)
				writer.write((point[0] - 300.84659082941437) + " " + (point[1] - 139.85962648478687) + " "
						+ (point[2] - 312.44765714375006) + " " + (point[3] - 154.44886993675695) + " "
						+ (length - 18.639494739141032));
			writer.write("\r\n");
		
       
			if (label == 7)
				writer.write((point[0] - 258.4602639611949) + " " + (point[1] - 161.8029504673369) + " "
						+ (point[2] - 270.22003738935496) + " " + (point[3] - 151.52609599285373) + " "
						+ (length - 15.61749048251234));
			writer.write("\r\n");
			if (label == 8)
				writer.write((point[0] - 10.782326104595313) + " " + (point[1] - 215.40174428753127) + " "
						+ (point[2] - 20.269038650282194) + " " + (point[3] - 202.82555475178765) + " "
						+ (length - 15.753039648380254));
			writer.write("\r\n");
			
			if (label == 10)
				writer.write((point[0] - 183.13644015784087) + " " + (point[1] - 258.1201972778982) + " "
						+ (point[2] - 194.29386588004672) + " " + (point[3] - 250.758719045703) + " "
						+ (length - 13.367105539705484));
			writer.write("\r\n");
			if (label == 11)
				writer.write((point[0] - 43.16752082802935) + " " + (point[1] - 266.8965342235055) + " "
						+ (point[2] - 55.11228730824559) + " " + (point[3] - 285.7189935156261) + " "
						+ (length - 22.292653948564176));
			writer.write("\r\n");
			if (label == 12)
				writer.write((point[0] - 89.19182304526953) + " " + (point[1] - 349.37985704768585) + " "
						+ (point[2] - 99.23449433610313) + " " + (point[3] - 334.4341802148607) + " "
						+ (length - 18.00634617702824));
			writer.write("\r\n");
			if (label == 13)
				writer.write((point[0] - 226.4867575369383) + " " + (point[1] - 336.7196803069454) + " "
						+ (point[2] - 236.9539633568867) + " " + (point[3] - 358.40376907956824) + " "
						+ (length - 24.078249595353665));
			writer.write("\r\n");
        
			if (label == 14)
				writer.write((point[0] - 36.854532534707985) + " " + (point[1] - 355.11968257966225) + " "
						+ (point[2] - 46.97940312410058) + " " + (point[3] - 371.014896388969) + " "
						+ (length - 18.845976400694227));
			writer.write("\r\n");

			if (label == 15)
				writer.write((point[0] - 319.91703640799966) + " " + (point[1] - 374.6571062489546) + " "
						+ (point[2] - 330.20918012825047) + " " + (point[3] - 382.22653982045375) + " "
						+ (length - 12.775936245600017));
			writer.write("\r\n");
			
			if (label == 16)
				writer.write((point[0] - 184.6352590369916) + " " + (point[1] - 381.94535949289735) + " "
						+ (point[2] - 193.745540933255) + " " + (point[3] - 386.7172540902593) + " "
						+ (length - 10.284367470958854));
			writer.write("\r\n");
			if (label == 17)
				writer.write((point[0] - 257.00195165698847) + " " + (point[1] - 393.0422893229684) + " "
						+ (point[2] - 267.69268843741247) + " " + (point[3] - 402.5256330267631) + " "
						+ (length - 14.290754378709815));
			writer.write("\r\n");
			if (label == 18)
				writer.write((point[0] - 279.74010616792725) + " " + (point[1] - 409.52143870609723) + " "
						+ (point[2] - 289.02938093556145) + " " + (point[3] - 395.9148965060625) + " "
						+ (length - 16.47509078426977));
			writer.write("\r\n");
			if (label == 19)
				writer.write((point[0] - 307.24157464889754) + " " + (point[1] - 425.6772245139631 ) + " "
						+ (point[2] - 318.90660598829754) + " " + (point[3] - 405.7625301506339) + " "
						+ (length - 23.079601550590034));
			writer.write("\r\n");
			
			
			writer.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	
	private void Testerrorfourth(final double[] point, int label, double length) {

		// Errorlist for Fake_big_file_nonoise.tif
		try {
			FileWriter writer = new FileWriter("error-Fake_fourth-offsettingLM.txt", true);
			
			
			if (label == 2)
				writer.write((point[0] - 92.9218278288335 ) + " " + (point[1] - 27.17362370877159) + " "
						+ (point[2] - 102.61190201305033  ) + " " + (point[3] - 8.259538706043248 ) + " "
						+ (length - 21.2518269611358));
			writer.write("\r\n");
			if (label == 3)
				writer.write((point[0] - 92.86686704851863) + " " + (point[1] - 64.90960165022722 ) + " "
						+ (point[2] -  104.48587735237139 ) + " " + (point[3] - 50.86749014591775 ) + " "
						+ (length - 18.22586886654507));
			writer.write("\r\n");
			if (label == 4)
				writer.write((point[0] - 280.4959354329863 ) + " " + (point[1] - 72.78270243124555) + " "
						+ (point[2] - 291.87059998038154) + " " + (point[3] - 91.7639867470957) + " "
						+ (length - 22.128536956719685));
			writer.write("\r\n");
			if (label == 5)
				writer.write((point[0] - 140.78720832889556 ) + " " + (point[1] - 132.0595306364064) + " "
						+ (point[2] - 150.13633156292963) + " " + (point[3] - 121.58435770257674) + " "
						+ (length - 14.04048977916345));
			writer.write("\r\n");
		
			if (label == 6)
				writer.write((point[0] - 177.83448360189217) + " " + (point[1] - 120.34554541214051 ) + " "
						+ (point[2] - 189.51595497473548) + " " + (point[3] - 144.65370360371278) + " "
						+ (length - 26.969303441154352));
			writer.write("\r\n");
			if (label == 7)
				writer.write((point[0] - 297.81778952607874) + " " + (point[1] - 174.6084120361677 ) + " "
						+ (point[2] - 308.05934208691633) + " " + (point[3] - 179.0082922914368) + " "
						+ (length - 11.146674172913903));
			writer.write("\r\n");
		
			if (label == 8)
				writer.write((point[0] - 293.94451588690464) + " " + (point[1] - 212.63092313007522) + " "
						+ (point[2] - 305.3832535507701) + " " + (point[3] -238.07718582965657) + " "
						+ (length - 27.89905024761309));
			writer.write("\r\n");
			if (label == 9)
				writer.write((point[0] - 177.5428835958507) + " " + (point[1] - 383.039314046178) + " "
						+ (point[2] - 187.97215142327738) + " " + (point[3] -397.12932856571365 ) + " "
						+ (length - 17.529921179997444));
			writer.write("\r\n");
			if (label == 10)
				writer.write((point[0] - 250.34852210673256) + " " + (point[1] - 383.14378596776106 ) + " "
						+ (point[2] - 262.3361618230495) + " " + (point[3] - 397.1830771085031) + " "
						+ (length - 18.460910099524806));
			writer.write("\r\n");
			
			
			
			writer.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	private void Testerrorsixth(final double[] point, int label, double length) {

		// Errorlist for Fake_big_file_nonoise.tif
		try {
			FileWriter writer = new FileWriter("error-Fake_sixth-offsettingLM.txt", true);
			
			if (label == 1)
				writer.write((point[0] - 169.64118749921616 ) + " " + (point[1] - 11.613782570008082) + " "
						+ (point[2] - 178.99307294259404  ) + " " + (point[3] - 15.89183603164018 ) + " "
						+ (length - 10.283943930547535));
			writer.write("\r\n");
			if (label == 2)
				writer.write((point[0] - 115.75116431067165 ) + " " + (point[1] - 55.6534204178164) + " "
						+ (point[2] - 125.81034889143847  ) + " " + (point[3] - 43.11496689247307 ) + " "
						+ (length - 16.074825387453867));
			writer.write("\r\n");
			if (label == 3)
				writer.write((point[0] - 28.220334047159405) + " " + (point[1] - 164.38434979659942 ) + " "
						+ (point[2] -  37.75929834994512 ) + " " + (point[3] - 180.07960823571474 ) + " "
						+ (length - 18.36662673003513));
			writer.write("\r\n");
			if (label == 4)
				writer.write((point[0] - 323.19957370820015 ) + " " + (point[1] - 279.68244187607655) + " "
						+ (point[2] - 333.76306541825585) + " " + (point[3] - 276.2791536595286) + " "
						+ (length - 11.098185788375933));
			writer.write("\r\n");
			if (label == 5)
				writer.write((point[0] - 121.79572577173442 ) + " " + (point[1] - 277.45720726092713) + " "
						+ (point[2] - 132.65834127497192) + " " + (point[3] - 287.29659315550555) + " "
						+ (length - 14.656395544389628));
			writer.write("\r\n");
		
			if (label == 6)
				writer.write((point[0] - 357.145870455239) + " " + (point[1] - 298.0757780155999 ) + " "
						+ (point[2] - 367.93878216586546) + " " + (point[3] - 315.5802625070824) + " "
						+ (length - 20.5643847587505));
			writer.write("\r\n");
			if (label == 7)
				writer.write((point[0] - 106.82706188254228 ) + " " + (point[1] - 317.70439690048266 ) + " "
						+ (point[2] - 117.1764204276652) + " " + (point[3] - 330.37225736794363) + " "
						+ (length - 16.357992270403155));
			writer.write("\r\n");
		
			if (label == 8)
				writer.write((point[0] - 36.1891360034064 ) + " " + (point[1] - 366.20177329461654) + " "
						+ (point[2] - 45.56449579796356) + " " + (point[3] - 382.0693576075372) + " "
						+ (length - 18.43034463066398));
			writer.write("\r\n");
			if (label == 9)
				writer.write((point[0] - 351.81427886479867) + " " + (point[1] - 371.74904382151175) + " "
						+ (point[2] - 363.0083416140858) + " " + (point[3] - 368.5496057444498 ) + " "
						+ (length - 11.642312693100601));
			writer.write("\r\n");
			if (label == 10)
				writer.write((point[0] - 341.2930905442583) + " " + (point[1] - 401.7183995097824 ) + " "
						+ (point[2] - 350.87179672046943) + " " + (point[3] - 392.6299117487737) + " "
						+ (length - 13.204250141230684));
			writer.write("\r\n");
			
			
			
			writer.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private void Testerrorseventh(final double[] point, int label, double length) {

		// Errorlist for Fake_big_file_nonoise.tif
		try {
			FileWriter writer = new FileWriter("error-Fake_seventh-offsettingLM.txt", true);
			
			if (label == 1)
				writer.write((point[0] - 271.4643571467613 ) + " " + (point[1] - 28.812548804883015) + " "
						+ (point[2] - 282.8951675969565  ) + " " + (point[3] - 51.710164559856814 ) + " "
						+ (length - 25.59226904380924));
			writer.write("\r\n");
			if (label == 2)
				writer.write((point[0] - 110.5306885004662 ) + " " + (point[1] - 65.44525892453427) + " "
						+ (point[2] - 121.62549571709461  ) + " " + (point[3] - 57.18401955845997 ) + " "
						+ (length - 13.83267230283887));
			writer.write("\r\n");
			if (label == 3)
				writer.write((point[0] - 68.05868325990714) + " " + (point[1] - 75.63455282822034 ) + " "
						+ (point[2] -  79.452538063318 ) + " " + (point[3] - 79.01775339391983 ) + " "
						+ (length - 11.88553630884859));
			writer.write("\r\n");
			if (label == 4)
				writer.write((point[0] - 206.6437991298365 ) + " " + (point[1] - 145.39791689745067) + " "
						+ (point[2] - 218.42214304937968) + " " + (point[3] - 144.92232783566837) + " "
						+ (length - 11.787941739028346));
			writer.write("\r\n");
			if (label == 5)
				writer.write((point[0] - 279.74510883639203 ) + " " + (point[1] - 185.48363393201532) + " "
						+ (point[2] - 288.9107332880293) + " " + (point[3] - 169.20120267841733) + " "
						+ (length - 18.68492009928315));
			writer.write("\r\n");
		
			if (label == 6)
				writer.write((point[0] - 88.40052110271773 ) + " " + (point[1] - 188.59400638940423 ) + " "
						+ (point[2] - 99.07470738911562) + " " + (point[3] - 207.37813796257564) + " "
						+ (length - 21.605134848802486));
			writer.write("\r\n");
			if (label == 7)
				writer.write((point[0] - 65.47844457154001 ) + " " + (point[1] - 252.69063416414554 ) + " "
						+ (point[2] - 76.51766577799131) + " " + (point[3] - 261.1045084130139) + " "
						+ (length - 13.880118325170509));
			writer.write("\r\n");
		
			if (label == 8)
				writer.write((point[0] - 0.48962832499659115 ) + " " + (point[1] - 321.2109153438377) + " "
						+ (point[2] - 11.288741744009124) + " " + (point[3] - 305.6177616973579) + " "
						+ (length - 18.96753255643506));
			writer.write("\r\n");
			if (label == 9)
				writer.write((point[0] - 292.58579113897326) + " " + (point[1] -378.7141000644622) + " "
						+ (point[2] - 302.9865792236662) + " " + (point[3] - 398.4700948971722 ) + " "
						+ (length - 22.326569924929256));
			writer.write("\r\n");
			if (label == 10)
				writer.write((point[0] - 172.44371203632338) + " " + (point[1] - 394.49432687477173 ) + " "
						+ (point[2] - 182.66983302824718) + " " + (point[3] - 399.9071680381849) + " "
						+ (length - 11.570324109626533));
			writer.write("\r\n");
			
			
			
			writer.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	private void Testerrorfirst(final double[] point, int label, double length) {

		// Errorlist for Fake_big_file_nonoise.tif
		try {
			FileWriter writer = new FileWriter("error-Fake_databigsnp-offsettingLM.txt", true);

			if (label == 1)
				writer.write((point[0] - 92.8676372377067) + " " + (point[1] - 5.499718994655133) + " "
						+ (point[2] - 112.8676372377067) + " " + (point[3] - 21.131534922414605) + " "
						+ (length - 25.384122383871276));
			writer.write("\r\n");
			if (label == 2)
				writer.write((point[0] - 147.99397954559853) + " " + (point[1] - 32.88119775551587) + " "
						+ (point[2] - 167.99397954559853) + " " + (point[3] - 45.984869141998246) + " "
						+ (length - 23.91037857928889));
			writer.write("\r\n");
			if (label == 3)
				writer.write((point[0] - 7.794350587783296) + " " + (point[1] - 33.10481325801554) + " "
						+ (point[2] - 27.794350587783295) + " " + (point[3] - 52.91969866168526) + " "
						+ (length - 28.153679751687225));
			writer.write("\r\n");
			if (label == 4)
				writer.write((point[0] - 55.346451860313366) + " " + (point[1] - 40.826958488144356) + " "
						+ (point[2] - 75.34645186031337) + " " + (point[3] - 56.4309289796462) + " "
						+ (length - 25.366984351705273));
			writer.write("\r\n");
			if (label == 5)
				writer.write((point[0] - 100.59945184657138) + " " + (point[1] - 47.24326969725762) + " "
						+ (point[2] - 120.59945184657138) + " " + (point[3] - 50.00073903143623) + " "
						+ (length - 20.18919604959384));
			writer.write("\r\n");
			if (label == 6)
				writer.write((point[0] - 46.95812894108671) + " " + (point[1] - 73.9193397629198) + " "
						+ (point[2] - 66.95812894108671) + " " + (point[3] - 78.3477854435783) + " "
						+ (length - 20.484411906289694));
			writer.write("\r\n");
			if (label == 7)
				writer.write((point[0] - 41.379541842397856) + " " + (point[1] - 80.073242656372) + " "
						+ (point[2] - 61.379541842397856 ) + " " + (point[3] - 94.4024072043764) + " "
						+ (length - 24.6033525488659));
			writer.write("\r\n");
			if (label == 8)
				writer.write((point[0] - 45.30959973137538) + " " + (point[1] - 158.02503112977763) + " "
						+ (point[2] - 65.30959973137539) + " " + (point[3] - 176.27699530569882) + " "
						+ (length - 27.07645095427225));
			writer.write("\r\n");
			if (label == 9)
				writer.write((point[0] - 106.06073273206955) + " " + (point[1] - 183.01279696011633 ) + " "
						+ (point[2] - 126.06073273206955) + " " + (point[3] - 184.70741244170682 ) + " "
						+ (length - 20.071664645226768));
			writer.write("\r\n");
			
			
			writer.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void Testerrorfifth(final double[] point, int label, double length) {

		// Errorlist for Fake_big_file_nonoise.tif
		try {
			FileWriter writer = new FileWriter("testerror-Fake_fifth-offsettingLM.txt", true);

			if (label == 1)
				writer.write((point[0] - 334.685446700782) + " " + (point[1] - 27.232991759972254) + " "
						+ (point[2] - 355.90213277915547) + " " + (point[3] - 71.43149034561719) + " "
						+ (length - 49.0270848141464));
			writer.write("\r\n");
			if (label == 2)
				writer.write((point[0] - 195.99940634144068) + " " + (point[1] - 60.110733648894104) + " "
						+ (point[2] - 216.0588472289147) + " " + (point[3] - 93.51548949286881) + " "
						+ (length - 38.96484161027254));
			writer.write("\r\n");
			if (label == 3)
				writer.write((point[0] - 292.7284350310164) + " " + (point[1] - 63.37806024717833) + " "
						+ (point[2] - 313.7854245963145) + " " + (point[3] - 86.83291705025314) + " "
						+ (length - 31.520265182986943));
			writer.write("\r\n");
			if (label == 4)
				writer.write((point[0] - 334.81180118812387) + " " + (point[1] - 175.71391696367218) + " "
						+ (point[2] - 356.0815811923552) + " " + (point[3] - 231.95156505253038) + " "
						+ (length - 60.12550709968821));
			writer.write("\r\n");
			if (label == 5)
				writer.write((point[0] - 0.5695104580630475) + " " + (point[1] - 180.11255226252794) + " "
						+ (point[2] - 21.90197486784737) + " " + (point[3] - 204.89056777814363) + " "
						+ (length - 32.69593385555468));
			writer.write("\r\n");
			if (label == 6)
				writer.write((point[0] - 288.45313521245055) + " " + (point[1] - 216.133940619811) + " "
						+ (point[2] - 309.3929746383774) + " " + (point[3] - 268.92640029469726) + " "
						+ (length - 56.79366754936761));
			writer.write("\r\n");
			if (label == 7)
				writer.write((point[0] - 337.326204913675) + " " + (point[1] - 250.106819309965) + " "
						+ (point[2] - 357.34629942387886 ) + " " + (point[3] - 307.2297384542443) + " "
						+ (length - 60.52959669253847));
			writer.write("\r\n");
			if (label == 8)
				writer.write((point[0] - 176.52956328690408) + " " + (point[1] - 278.30082154667275) + " "
						+ (point[2] - 197.30680558215627) + " " + (point[3] - 334.56959080632765) + " "
						+ (length - 59.982232297505405));
			writer.write("\r\n");
			if (label == 9)
				writer.write((point[0] - 6.653366259279758) + " " + (point[1] - 306.16197030167274 ) + " "
						+ (point[2] - 28.40503221617591) + " " + (point[3] - 349.77300714455924 ) + " "
						+ (length - 48.73456172381157));
			writer.write("\r\n");
			if (label == 10)
				writer.write((point[0] - 344.7227374264369) + " " + (point[1] - 389.3751531282371 ) + " "
						+ (point[2] - 366.6881936690487) + " " + (point[3] - 410.85382417972704 ) + " "
						+ (length - 30.721565358623998));
			writer.write("\r\n");
			
			writer.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	
	private void Testerrorsec(final double[] point, int label, double length) {

		// Errorlist for Fake_big_file_nonoise.tif
		try {
			FileWriter writer = new FileWriter("error-Fake_bigfile_noisy_sec-offsettingLM.txt", true);
			
			if (label == 1)
				writer.write((point[0] - 305.9925929042284) + " " + (point[1] - 6.386481793148133) + " "
						+ (point[2] - 322.0288782573045) + " " + (point[3] - 6.026247216588337) + " "
						+ (length - 16.040330946568158));
			writer.write("\r\n");
			
			if (label == 3)
				writer.write((point[0] - 353.15723234722964) + " " + (point[1] - 49.99143903881524) + " "
						+ (point[2] - 359.4669599689136) + " " + (point[3] - 43.039872569086285) + " "
						+ (length - 9.388127557873378));
			writer.write("\r\n");
			
			if (label == 5)
				writer.write((point[0] - 282.92829236420187) + " " + (point[1] - 51.861712067783714) + " "
						+ (point[2] - 294.6268069947843) + " " + (point[3] - 47.98838229640566) + " "
						+ (length - 12.323064881748989));
			writer.write("\r\n");
			if (label == 6)
				writer.write((point[0] - 150.17254613234996) + " " + (point[1] - 56.09908721936941) + " "
						+ (point[2] - 165.71900210752256) + " " + (point[3] - 52.70967491949685) + " "
						+ (length - 15.911643822261324));
			writer.write("\r\n");
			if (label == 7)
				writer.write((point[0] - 126.39749395604314) + " " + (point[1] - 91.10615287483455) + " "
						+ (point[2] - 144.3712600366361) + " " + (point[3] - 82.1183621577551 ) + " "
						+ (length - 20.095687325739664));
			writer.write("\r\n");
			if (label == 8)
				writer.write((point[0] - 185.06081941183263) + " " + (point[1] - 100.42682277520953) + " "
						+ (point[2] - 190.71523125712804) + " " + (point[3] - 90.4882958233705) + " "
						+ (length - 11.4344519452682));
			writer.write("\r\n");
			if (label == 9)
				writer.write((point[0] - 298.2248630132484) + " " + (point[1] - 124.5154128690997) + " "
						+ (point[2] - 311.81988528482725) + " " + (point[3] - 119.4678472985486) + " "
						+ (length - 14.50181191278309));
			writer.write("\r\n");
			if (label == 10)
				writer.write((point[0] - 207.15405474164575) + " " + (point[1] - 150.6899404519233) + " "
						+ (point[2] - 208.44789551107206) + " " + (point[3] - 139.5953917951248) + " "
						+ (length - 11.169737411179309));
			writer.write("\r\n");
			if (label == 11)
				writer.write((point[0] - 112.7552499916626) + " " + (point[1] - 166.99368119375615) + " "
						+ (point[2] - 119.33291011030018) + " " + (point[3] - 158.25754114246612) + " "
						+ (length - 10.935527222409945));
			writer.write("\r\n");
			if (label == 12)
				writer.write((point[0] - 236.22859339630475) + " " + (point[1] - 184.74586622212533) + " "
						+ (point[2] - 253.52499503650682) + " " + (point[3] - 174.24640098382434) + " "
						+ (length - 20.233741126877053));
			writer.write("\r\n");
			if (label == 13)
				writer.write((point[0] -118.41559571666849) + " " + (point[1] - 211.7010231662632) + " "
						+ (point[2] - 129.451058044193) + " " + (point[3] - 207.22262862709798) + " "
						+ (length - 11.909552738480041));
			writer.write("\r\n");
			if (label == 14)
				writer.write((point[0] - 239.43550593841613) + " " + (point[1] - 217.8967196461199) + " "
						+ (point[2] - 251.85701790196183) + " " + (point[3] - 209.36985573134737) + " "
						+ (length - 15.066564561357664));
			writer.write("\r\n");
			if (label == 15)
				writer.write((point[0] - 372.5730494703797 ) + " " + (point[1] -  219.2225283527957) + " "
						+ (point[2] - 387.07871048759705 ) + " " + (point[3] - 213.57221893243945) + " "
						+ (length - 15.567279726791872));
			writer.write("\r\n");
			if (label == 16)
				writer.write((point[0] - 144.8838518365699) + " " + (point[1] - 241.74727194739705) + " "
						+ (point[2] - 149.82141674236192) + " " + (point[3] - 232.2334431268912) + " "
						+ (length - 10.718791257637072));
			writer.write("\r\n");
			if (label == 17)
				writer.write((point[0] - 62.933723402780544) + " " + (point[1] - 257.46201405936245) + " "
						+ (point[2] - 73.24523552296085) + " " + (point[3] - 248.03191210113425) + " "
						+ (length - 13.973335505426206));
			writer.write("\r\n");
			if (label == 18)
				writer.write((point[0] - 51.588632388998214) + " " + (point[1] - 267.74624335786905) + " "
						+ (point[2] - 69.53297671887874) + " " + (point[3] - 270.33575591613373) + " "
						+ (length - 18.1302252804185));
			writer.write("\r\n");
			if (label == 19)
				writer.write((point[0] - 371.5598843638897) + " " + (point[1] - 302.3070640302173 ) + " "
						+ (point[2] - 382.3437183288822) + " " + (point[3] - 291.4978977236864) + " "
						+ (length - 15.268567425491131));
			writer.write("\r\n");
			if (label == 20)
				writer.write((point[0] - 300.6076952914607) + " " + (point[1] - 360.193761994309 ) + " "
						+ (point[2] - 312.80245594393347) + " " + (point[3] - 350.7875453664556) + " "
						+ (length - 15.400944731450855));
			writer.write("\r\n");
			if (label == 21)
				writer.write((point[0] - 325.8279550275642) + " " + (point[1] - 396.94385867009623 ) + " "
						+ (point[2] - 330.78480705510646 ) + " " + (point[3] - 385.9961289446052) + " "
						+ (length - 12.01761907223516));
			if (label == 22)
				writer.write((point[0] - 193.05380388784238 ) + " " + (point[1] - 394.24805388470986 ) + " "
						+ (point[2] - 209.00212892787482 ) + " " + (point[3] - 387.0908873699319) + " "
						+ (length - 17.48067802182697));
			writer.write("\r\n");
			if (label == 23)
				writer.write((point[0] - 363.18051478279176 ) + " " + (point[1] - 395.6581626679954 ) + " "
						+ (point[2] - 371.2717339920776  ) + " " + (point[3] - 387.40066285914924 ) + " "
						+ (length - 11.560888001611742));
			writer.write("\r\n");
			if (label == 24)
				writer.write((point[0] - 257.56853558741227 ) + " " + (point[1] - 392.1317053778325 ) + " "
						+ (point[2] - 273.208024381375 ) + " " + (point[3] - 395.6612398320712  ) + " "
						+ (length - 16.03281707000192));
			writer.write("\r\n");
			if (label == 25)
				writer.write((point[0] - 242.17648532439642 ) + " " + (point[1] - 401.6586191271966 ) + " "
						+ (point[2] - 255.71765961466426   ) + " " + (point[3] - 399.89316172057016 ) + " "
						+ (length - 13.655776836709897));
			writer.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void Testerrorthird(final double[] point, int label, double length) {

		// Errorlist for Fake_big_file_nonoise.tif
		try {
			FileWriter writer = new FileWriter("error-Fake_bigfile_noisy_third-offsettingLM.txt", true);
			
			if (label == 1)
				writer.write((point[0] - 159.70876069055475) + " " + (point[1] - 14.145735424295653) + " "
						+ (point[2] - 167.61212303542496) + " " + (point[3] - 3.8033283013758235) + " "
						+ (length - 13.016471159670546));
			writer.write("\r\n");
			if (label == 2)
				writer.write((point[0] - 136.20057982375772) + " " + (point[1] - 5.207904433739971) + " "
						+ (point[2] - 152.86561656707167) + " " + (point[3] - 7.9806486301158355) + " "
						+ (length - 16.89412797496633));
			writer.write("\r\n");
			if (label == 3)
				writer.write((point[0] - 47.45429996121627) + " " + (point[1] - 25.63722823136536) + " "
						+ (point[2] - 57.27725022588619) + " " + (point[3] - 14.891294327282901) + " "
						+ (length - 14.55903318813058));
			writer.write("\r\n");
			if (label == 4)
				writer.write((point[0] - 335.8381671659929) + " " + (point[1] - 24.31304234800028) + " "
						+ (point[2] - 352.92772515391044 ) + " " + (point[3] - 27.088979208768457) + " "
						+ (length - 17.31354434185467));
			writer.write("\r\n");
			if (label == 5)
				writer.write((point[0] - 177.3827474288163) + " " + (point[1] - 36.48388865262402) + " "
						+ (point[2] - 186.42256997573605 ) + " " + (point[3] - 25.703042029276453) + " "
						+ (length - 14.069294424239743));
			writer.write("\r\n");
			if (label == 6)
				writer.write((point[0] - 301.2171132081411 ) + " " + (point[1] - 30.315498396371794) + " "
						+ (point[2] - 317.15918882436847 ) + " " + (point[3] - 33.18008496056557) + " "
						+ (length - 16.197395813440853));
			writer.write("\r\n");
			if (label == 7)
				writer.write((point[0] - 335.03929214268237) + " " + (point[1] - 165.69191076004964) + " "
						+ (point[2] - 342.5256140270755 ) + " " + (point[3] - 157.48931895989702 ) + " "
						+ (length - 11.105292774018817));
			writer.write("\r\n");
			if (label == 8)
				writer.write((point[0] - 67.96287883132905) + " " + (point[1] - 173.86588761653982) + " "
						+ (point[2] - 77.93444046791349) + " " + (point[3] - 164.38574542277735) + " "
						+ (length - 13.758820352274325));
			writer.write("\r\n");
			if (label == 9)
				writer.write((point[0] - 132.08500177775954) + " " + (point[1] - 191.2328208796786) + " "
						+ (point[2] - 133.99927669679693) + " " + (point[3] - 179.71223488708554) + " "
						+ (length - 11.67854229252893));
			writer.write("\r\n");
			if (label == 10)
				writer.write((point[0] - 240.43561617353788 ) + " " + (point[1] - 215.53133293573083) + " "
						+ (point[2] - 255.6876330404176) + " " + (point[3] -205.4850084015699) + " "
						+ (length - 18.2634239712401));
			writer.write("\r\n");
			if (label == 11)
				writer.write((point[0] - 282.1835439050488 ) + " " + (point[1] - 253.03694030233416) + " "
						+ (point[2] - 292.61517066888075) + " " + (point[3] - 246.2773467971297) + " "
						+ (length - 12.43024301031551));
			writer.write("\r\n");
			if (label == 12)
				writer.write((point[0] - 294.153062640094) + " " + (point[1] - 261.93721368868546) + " "
						+ (point[2] - 302.13372226040036) + " " + (point[3] - 252.38483925060393) + " "
						+ (length - 12.44744091693156));
			writer.write("\r\n");
			if (label == 13)
				writer.write((point[0] -143.94147746939367) + " " + (point[1] - 297.3714118164239) + " "
						+ (point[2] - 157.24715327889555) + " " + (point[3] - 293.39652929013005) + " "
						+ (length - 13.886709467883303));
			writer.write("\r\n");
			if (label == 14)
				writer.write((point[0] - 305.2899421104744) + " " + (point[1] - 311.5437229719812) + " "
						+ (point[2] - 312.97810052264117) + " " + (point[3] - 302.27436866967133 ) + " "
						+ (length - 12.04278659415342));
			writer.write("\r\n");
			if (label == 15)
				writer.write((point[0] - 76.28332921728651 ) + " " + (point[1] -  326.1859505945892) + " "
						+ (point[2] - 84.60177818692252 ) + " " + (point[3] - 316.258685571946) + " "
						+ (length - 12.951725139541594));
			writer.write("\r\n");
			if (label == 16)
				writer.write((point[0] - 135.61245447425634) + " " + (point[1] - 316.38039979427714) + " "
						+ (point[2] - 137.49991384966583) + " " + (point[3] - 328.74893924805554) + " "
						+ (length - 12.511725345190966));
			writer.write("\r\n");
			if (label == 17)
				writer.write((point[0] - 38.84748653022617) + " " + (point[1] - 335.8949791635391) + " "
						+ (point[2] - 41.53933477646034) + " " + (point[3] - 324.41375518435126) + " "
						+ (length - 11.79256337871594));
			writer.write("\r\n");
			if (label == 18)
				writer.write((point[0] - 99.35229992942286) + " " + (point[1] - 339.2657583667901) + " "
						+ (point[2] - 109.0361829873792) + " " + (point[3] - 327.86290858949945) + " "
						+ (length - 14.96003255757259));
			writer.write("\r\n");
			if (label == 19)
				writer.write((point[0] - 232.91604534730675) + " " + (point[1] - 333.41185260962834 ) + " "
						+ (point[2] - 243.45367569566125) + " " + (point[3] - 328.4354683761082) + " "
						+ (length - 11.653585431024636));
			writer.write("\r\n");
			if (label == 20)
				writer.write((point[0] - 388.7911761490889) + " " + (point[1] - 339.7061239652844 ) + " "
						+ (point[2] - 395.4960477552932) + " " + (point[3] - 329.50112703979767) + " "
						+ (length - 12.210539116061916));
			writer.write("\r\n");
			if (label == 21)
				writer.write((point[0] - 315.65300280159875) + " " + (point[1] - 362.6874993570169 ) + " "
						+ (point[2] - 333.3737841658348  ) + " " + (point[3] - 354.04717540951754 ) + " "
						+ (length - 19.715001650438374));
			if (label == 22)
				writer.write((point[0] - 196.19090817342362 ) + " " + (point[1] - 371.8379923415543 ) + " "
						+ (point[2] - 198.92909691454048 ) + " " + (point[3] - 359.8836219091283) + " "
						+ (length - 12.263957355504797));
			writer.write("\r\n");
			
			if (label == 24)
				writer.write((point[0] - 335.4158263232457 ) + " " + (point[1] - 384.3346481250782 ) + " "
						+ (point[2] - 340.80011140941934 ) + " " + (point[3] - 376.6039312635851  ) + " "
						+ (length - 9.420961154880434));
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
}

package peakFitter;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import houghandWatershed.Boundingboxes;
import labeledObjects.Simpleobject;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.Point;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.BenchmarkAlgorithm;
import net.imglib2.algorithm.OutputAlgorithm;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import peakFitter.GaussianMaskFit.Endfit;
import preProcessing.GetLocalmaxmin;

public class SubpixelBentLength extends BenchmarkAlgorithm
implements OutputAlgorithm<ArrayList<double[]>> {
	
	private static final String BASE_ERROR_MSG = "[SubpixelLine] ";
	private final RandomAccessibleInterval<FloatType> source;
	private final RandomAccessibleInterval<IntType> intimg;
	private final ArrayList<Simpleobject> simpleobject;
	private final int ndims;
	private ArrayList<double[]> final_paramlist;
	private final double[] psf;
	private final double minlength;
	// LM solver iteration params
	final int maxiter = 500;
	final double lambda = 1e-3;
	final double termepsilon = 1e-1;
	//Mask fits iteration param
	final int iterations = 500;
	final double cutoffdistance = 20;
	final boolean halfgaussian = false;
	final double Intensityratio = 0.5;
	
	
	public SubpixelBentLength( final RandomAccessibleInterval<FloatType> source, 
			             final RandomAccessibleInterval<IntType> intimg,
			             final ArrayList<Simpleobject> simpleobject,
			             final double[] psf,
			             final double minlength){
		
		this.source = source;
		this.intimg = intimg;
		this.simpleobject = simpleobject;
		this.psf = psf;
		this.minlength = minlength;
		this.ndims = source.numDimensions();
		assert (source.numDimensions() == intimg.numDimensions());
		
	}

	@Override
	public boolean checkInput() {
		if (source.numDimensions() > 2) {
			errorMessage = BASE_ERROR_MSG + " Can only operate on 2D, make slices of your stack . Got "
					+ source.numDimensions() + "D.";
			return false;
		}
		return true;
	}

	@Override
	public boolean process() {
		
		final_paramlist = new ArrayList<double[]>();
		for (int index = 0; index < simpleobject.size(); ++index) {
			
			final int Label = simpleobject.get(index).Label;
			final double slope = simpleobject.get(index).slope;
			final double intercept = simpleobject.get(index).intercept;
			final double [] final_param = Getfinallineparam(Label, slope, intercept, psf, minlength);
			if (final_param!= null)
			final_paramlist.add(final_param);
		}
		
		

		return true;
	}

	@Override
	public ArrayList<double[]> getResult() {
		
		return final_paramlist;
	}

	
	private final double[] MakeimprovedLineguess(double slope, double intercept, double[] psf, int label)  {
		long[] newposition = new long[ndims];
		double[] minVal = { Double.MAX_VALUE, Double.MAX_VALUE };
		double[] maxVal = { -Double.MIN_VALUE, -Double.MIN_VALUE };

		RandomAccessibleInterval<FloatType> currentimg = Boundingboxes.CurrentLabelImage(intimg, source,
				label);
		

		double newintercept = intercept;


		final Cursor<FloatType> outcursor = Views.iterable(currentimg).localizingCursor();
		final RandomAccess<IntType> intcursor = intimg.randomAccess();
		final double maxintensityline = GetLocalmaxmin.computeMaxIntensity(currentimg);

		while (outcursor.hasNext()) {

			outcursor.fwd();
            intcursor.setPosition(outcursor);
            if (outcursor.get().get()/maxintensityline > Intensityratio){
            if (intcursor.get().get() == label){
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
		}

		final double[] MinandMax = new double[2 * ndims + 4];

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
		MinandMax[2 * ndims] =   Math.min(psf[0], psf[1]);
		MinandMax[2 * ndims + 1] = maxintensityline; 
		// This parameter guess estimates the background noise level
		MinandMax[2 * ndims + 2] = 1; 
		MinandMax[2 * ndims + 3] = 0; 
		
		System.out.println("Label: " + label + " " + "Hough Detection: " + " StartX: " + MinandMax[0] + " StartY: "
				+ MinandMax[1] + " EndX: " + MinandMax[2] + " EndY: " + MinandMax[3] + " " + "Bent param" + MinandMax[2 * ndims + 3]  );

		
		
			for (int d = 0; d < ndims; ++d) {

				if (MinandMax[d] == Double.MAX_VALUE || MinandMax[d + ndims] == -Double.MIN_VALUE)
					return null;
				if (MinandMax[d] >= source.dimension(d) || MinandMax[d + ndims] >= source.dimension(d))
					return null;
				if (MinandMax[d] <= 0 || MinandMax[d + ndims] <= 0)
					return null;

			}
		

		

		

		return MinandMax;

	}
	
	// Get line parameters for fitting line to a line in a label

		public double[] Getfinallineparam(final int label, final double slope, final double intercept, final double[] psf,
				final double minlength)  {

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

			final double[] start_param = MakeimprovedLineguess(slope, intercept, psf, label);
			if (start_param == null)
				return null;

			else {

				final double[] finalparamstart = start_param.clone();
				// LM solver part

				RandomAccessibleInterval<FloatType> currentimg = Boundingboxes.CurrentLabelImage(intimg, source,
						label);
				

				final double[] fixed_param = new double[ndims];

				for (int d = 0; d < ndims; ++d) {

					fixed_param[d] = 1.0 / Math.pow(psf[d], 2);
				}

			
				double inicutoffdistance = 0;
				final double[] inistartpos = { start_param[0], start_param[1] };
				final double[] iniendpos = { start_param[2], start_param[3] };
				inicutoffdistance = Distance(inistartpos, iniendpos);
				
				
				
				
				
				if (inicutoffdistance > minlength) {
					try {
						LevenbergMarquardtSolverLine.solve(X, finalparamstart, fixed_param, I, new GaussianBentLineds(), lambda,
								termepsilon, maxiter);
					} catch (Exception e) {
						e.printStackTrace();
					}

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
					try {
						startfit = peakFitter.GaussianMaskFit.sumofgaussianMaskFit(currentimg, startpos.clone(), psf,
								iterations, dxvector, newslope, newintercept, maxintensityline, halfgaussian, Endfit.Start,
								label);
					} catch (Exception e) {
						e.printStackTrace();
					}

					try {
						endfit = peakFitter.GaussianMaskFit.sumofgaussianMaskFit(currentimg, endpos.clone(), psf,
								iterations, dxvector, newslope, newintercept, maxintensityline,  halfgaussian, Endfit.End,
								label);
					} catch (Exception e) {
						e.printStackTrace();
					}

					final double Maskdist = sqDistance(startfit, endfit);
					// If mask fits fail, return LM solver results, very crucial for
					// noisy data

					double[] returnparam = new double[2 * ndims + 4];
					double[] LMparam = new double[2 * ndims];

					

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
					

					System.out.println("Fits :" + "StartX:" + returnparam[0] + " StartY:" + returnparam[1] + " " + "EndX:"
							+ returnparam[2] + "EndY: " + returnparam[3] + " " + "ds: " + finalparamstart[4] );

					System.out.println("Length: " + Distance(new double[]{returnparam[0],  returnparam[1]},new double[]{returnparam[2],  returnparam[3]} ));
					Testerrortwo(returnparam, label, Distance(new double[] { returnparam[0], returnparam[1] },
							new double[] { returnparam[2], returnparam[3] }));
					returnparam[2 * ndims] = finalparamstart[4];
					returnparam[2 * ndims + 1] = finalparamstart[5];
					returnparam[2 * ndims + 2] = finalparamstart[6];
					returnparam[2 * ndims + 3] = finalparamstart[7];
					return returnparam;

				}

				else
					return null;

			}
		}
		public int Getlabel(final Point linepoint) {

			RandomAccess<IntType> intranac = intimg.randomAccess();

			intranac.setPosition(linepoint);
			int currentlabel = intranac.get().get();

			return currentlabel;
		}


		

		private PointSampleList<FloatType> gatherfullData(final int label) {
			final PointSampleList<FloatType> datalist = new PointSampleList<FloatType>(ndims);

			RandomAccessibleInterval<FloatType> currentimg = Boundingboxes.CurrentLabelImage(intimg, source,
					label);
			
			boolean outofbounds = false;

			Cursor<FloatType> localcursor = Views.iterable(currentimg).localizingCursor();

			while (localcursor.hasNext()) {
				localcursor.fwd();

				for (int d = 0; d < ndims; d++) {

					if (localcursor.getDoublePosition(d) < 0 || localcursor.getDoublePosition(d) >= source.dimension(d)) {
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
		private void Testerrorone(final double[] point, int label, double length) {

			try {
				FileWriter writer = new FileWriter("../res/error-Pnoise1snr15Hough.txt", true);

				if (label == 1)
					writer.write((point[0] - 331.25407587479117) + " " + (point[1] - 20.852867366064714) + " "
							+ (point[2] - 336.0544738270745) + " " + (point[3] - 4.544704394726786 ) + " "
							+ (length - 17));
				writer.write("\r\n");
				if (label == 2)
					writer.write((point[0] - 357.7772829569185) + " " + (point[1] - 135.0526770048757) + " "
							+ (point[2] - 376.66554377520197) + " " + (point[3] - 120.24613779881884) + " "
							+ (length - 24));
				writer.write("\r\n");
				if (label == 3)
					writer.write((point[0] - 29.54178093981105) + " " + (point[1] - 134.020370559775) + " "
							+ (point[2] - 44.831239399479294) + " " + (point[3] - 138.7355015240304) + " "
							+ (length - 16));
				writer.write("\r\n");
				if (label == 4)
					writer.write((point[0] - 243.98675435953902) + " " + (point[1] - 137.3367331747739) + " "
							+ (point[2] - 268.95343229596267) + " " + (point[3] - 138.62707918220516) + " "
							+ (length - 24.999));
				writer.write("\r\n");
				if (label == 5)
					writer.write((point[0] - 350.117215979556) + " " + (point[1] - 149.37645829375725) + " "
							+ (point[2] - 365.03886102752364) + " " + (point[3] - 141.23124441918836) + " "
							+ (length - 17));
				writer.write("\r\n");
				if (label == 6)
					writer.write((point[0] - 398.4641506844532) + " " + (point[1] - 216.48151022847412 ) + " "
							+ (point[2] - 408.899852747816) + " " + (point[3] - 203.06152056958365) + " "
							+ (length - 17));
				writer.write("\r\n");

				if (label == 7)
					writer.write((point[0] - 107.82091544283591) + " " + (point[1] - 216.96614685962092) + " "
							+ (point[2] - 128.92832030710323) + " " + (point[3] - 210.7631266358985) + " "
							+ (length - 22));
				writer.write("\r\n");
				if (label == 8)
					writer.write((point[0] - 149.17307807786753) + " " + (point[1] - 259.7048195659357) + " "
							+ (point[2] - 164.23612622612032) + " " + (point[3] - 238.51273188011726) + " "
							+ (length - 25.999));
				writer.write("\r\n");
				
				if (label == 9)
					writer.write((point[0] - 4.129798064076837) + " " + (point[1] - 241.90174176947804) + " "
							+ (point[2] - 19.72450661678264) + " " + (point[3] - 260.1446873966559) + " "
							+ (length - 24));
				writer.write("\r\n");
				if (label == 10)
					writer.write((point[0] - 406.56643077926657) + " " + (point[1] - 285.42517366618745) + " "
							+ (point[2] - 423.02540981055705) + " " + (point[3] - 272.3828575663706) + " "
							+ (length - 21));
				writer.write("\r\n");
				if (label == 11)
					writer.write((point[0] - 188.33331401418496) + " " + (point[1] - 361.0687898756705) + " "
							+ (point[2] - 201.19360936039828) + " " + (point[3] - 343.21906148535084) + " "
							+ (length - 22));
				writer.write("\r\n");
				
				if (label == 13)
					writer.write((point[0] - 242.27073750693182 ) + " " + (point[1] - 401.0215776223561) + " "
							+ (point[2] - 257.94695655884016) + " " + (point[3] - 416.4571260815654 ) + " "
							+ (length - 22));
				writer.write("\r\n");

				if (label == 14)
					writer.write((point[0] - 128.1794901922759) + " " + (point[1] - 423.08894998368635) + " "
							+ (point[2] - 149.15996670775016) + " " + (point[3] - 438.44536900582943) + " "
							+ (length - 26));
				writer.write("\r\n");

				if (label == 15)
					writer.write((point[0] -383.5999637388662) + " " + (point[1] - 428.80166412226083) + " "
							+ (point[2] - 395.76151966251473) + " " + (point[3] - 447.13460053506606) + " "
							+ (length - 21.999));
				writer.write("\r\n");
				if (label == 17)
					writer.write((point[0] - 451.7368111880984) + " " + (point[1] - 431.26071471752243) + " "
							+ (point[2] - 480.50980808810493) + " " + (point[3] - 439.7527491847318) + " "
							+ (length - 30));
				writer.write("\r\n");
				if (label == 16)
					writer.write((point[0] - 211.8768803373426) + " " + (point[1] - 430.4545540295625) + " "
							+ (point[2] - 218.0397374553543) + " " + (point[3] - 447.36665597293324 ) + " "
							+ (length - 17.999));
				writer.write("\r\n");
				
				if (label == 18)
					writer.write((point[0] - 132.33597530343263) + " " + (point[1] -481.09798518247544) + " "
							+ (point[2] - 143.55251539520583) + " " + (point[3] - 503.4405281973079) + " "
							+ (length - 24.999));
				writer.write("\r\n");
				if (label == 19)
					writer.write((point[0] - 15.635640299061473) + " " + (point[1] -506.5683812906543) + " "
							+ (point[2] - 32.41110113306011) + " " + (point[3] - 503.81449062343876) + " "
							+ (length - 16.999));
				writer.write("\r\n");
				if (label == 20)
					writer.write((point[0] - 107.40037876264329 ) + " " + (point[1] - 525.4010542967413) + " "
							+ (point[2] - 129.18791411562972) + " " + (point[3] - 515.3361000870906 ) + " "
							+ (length - 24));
				writer.write("\r\n");
				
				

				writer.close();

			} catch (IOException e) {
				e.printStackTrace();
			}
		
		}

		private void Testerrortwo(final double[] point, int label, double length) {

			// Errorlist for Fake_big_file.tif
			try {
				FileWriter writer = new FileWriter("../res/error-Pnoise2snr15Hough.txt", true);

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
				
				if (label == 13)
					writer.write((point[0] - 221.05784248127878) + " " + (point[1] - 287.98083384565257) + " "
							+ (point[2] - 246.92398516835158) + " " + (point[3] - 301.09352474004817) + " "
							+ (length - 29));
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
		private void Testerrorthree(final double[] point, int label, double length) {

			// Errorlist for Fake_big_file.tif
			try {
				FileWriter writer = new FileWriter("../res/error-Pnoise3snr45Hough.txt", true);

		if (label == 1)
			writer.write((point[0] -120.8286527722272) + " " + (point[1] - 35.90331959617559 ) + " "
					+ (point[2] - 128.3358379749917) + " " + (point[3] - 11.010706829907168) + " "
					+ (length - 25.99999));
		writer.write("\r\n");
		if (label == 2)
			writer.write((point[0] - 383.1094865028246) + " " + (point[1] - 24.211622828136434) + " "
					+ (point[2] - 397.58964339022407) + " " + (point[3] - 15.305161324469086 ) + " "
					+ (length - 17));
		writer.write("\r\n");
		if (label == 3)
			writer.write((point[0] - 461.6064166590004) + " " + (point[1] - 32.56142533075916) + " "
					+ (point[2] - 469.6329721344923) + " " + (point[3] - 19.889629593814597) + " "
					+ (length - 14.99999999));
		writer.write("\r\n");
		if (label == 4)
			writer.write((point[0] - 130.67281403216373) + " " + (point[1] - 67.94708537739378) + " "
					+ (point[2] - 145.10413223838316) + " " + (point[3] - 81.79400037066326) + " "
					+ (length - 19.99999999));
		writer.write("\r\n");
		if (label == 5)
			writer.write((point[0] - 65.55650288878368) + " " + (point[1] -70.37778387737362) + " "
					+ (point[2] - 76.18895202493668) + " " + (point[3] - 93.0041186527437) + " "
					+ (length - 25));
		writer.write("\r\n");
		if (label == 6)
			writer.write((point[0] - 365.0211262912541 ) + " " + (point[1] - 122.57525237020812) + " "
					+ (point[2] - 365.5585088819804 ) + " " + (point[3] - 103.58285336557503) + " "
					+ (length - 19));
		writer.write("\r\n");
		if (label == 7)
			writer.write((point[0] - 456.13362283344867) + " " + (point[1] - 139.8294237980426) + " "
					+ (point[2] - 457.9747823668729 ) + " " + (point[3] - 158.74000632668202) + " "
					+ (length - 19));
		writer.write("\r\n");
		
		if (label == 8)
			writer.write((point[0] - 187.7644848368468) + " " + (point[1] - 162.4083026246098) + " "
					+ (point[2] - 205.57157794640466) + " " + (point[3] - 159.78010267263954) + " "
					+ (length - 17.999999));
		writer.write("\r\n");
		if (label == 9)
			writer.write((point[0] - 227.09840196562163) + " " + (point[1] - 211.67634722341552) + " "
					+ (point[2] - 237.03830634648833) + " " + (point[3] - 193.17774454918006) + " "
					+ (length - 21));
		writer.write("\r\n");
		if (label == 10)
			writer.write((point[0] - 364.53298941799824) + " " + (point[1] - 226.10175493456273) + " "
					+ (point[2] - 385.708753416097) + " " + (point[3] - 239.39035979349333) + " "
					+ (length - 24.99999));
		writer.write("\r\n");
		if (label == 11)
			writer.write((point[0] - 302.5742317731157) + " " + (point[1] - 326.989833342397) + " "
					+ (point[2] - 315.1970359599219) + " " + (point[3] - 314.157654475909) + " "
					+ (length -17.99999999));
		writer.write("\r\n");
		if (label == 12)
			writer.write((point[0] - 367.8641291577637) + " " + (point[1] - 324.33109676872925) + " "
					+ (point[2] - 382.4586555012279) + " " + (point[3] - 317.77367343254616) + " "
					+ (length - 15.9999999));
		writer.write("\r\n");
		if (label == 13)
			writer.write((point[0] - 2.1268212268463365) + " " + (point[1] - 389.72581922926076) + " "
					+ (point[2] - 23.40689682082701) + " " + (point[3] - 368.57969942645775) + " "
					+ (length - 30));
		writer.write("\r\n");
		if (label == 14)
			writer.write((point[0] - 86.01548609621949) + " " + (point[1] - 433.1702050680524) + " "
					+ (point[2] - 99.25674261766153) + " " + (point[3] - 418.1812382679907) + " "
					+ (length - 19.99999999));
		writer.write("\r\n");

		if (label == 15)
			writer.write((point[0] - 255.7207889127908 ) + " " + (point[1] - 427.2351699441273 ) + " "
					+ (point[2] - 270.3455272306394) + " " + (point[3] - 448.73207737599887) + " "
					+ (length - 26));
		writer.write("\r\n");

		

		if (label == 16)
			writer.write((point[0] - 107.60082443095591 ) + " " + (point[1] - 470.35947742086205) + " "
					+ (point[2] - 117.59361709892363) + " " + (point[3] - 443.1355156464434 ) + " "
					+ (length - 28.999));
		writer.write("\r\n");
		if (label == 17)
			writer.write((point[0] - 54.80889136248092) + " " + (point[1] - 520.5932127634953) + " "
					+ (point[2] - 63.19008609304531) + " " + (point[3] - 501.33819848439543) + " "
					+ (length - 20.999999));
		writer.write("\r\n");
		if (label == 18)
			writer.write((point[0] - 328.08457330752304) + " " + (point[1] - 530.5438183274165) + " "
					+ (point[2] - 333.4959718939446) + " " + (point[3] - 504.0916583902992 ) + " "
					+ (length - 26.9999999));
		writer.write("\r\n");
		if (label == 19)
			writer.write((point[0] - 235.08296631632936) + " " + (point[1] - 530.472022828654) + " "
					+ (point[2] - 238.0281000059525) + " " + (point[3] - 513.7290784176422) + " "
					+ (length - 16.999999));
		writer.write("\r\n");
		if (label == 20)
			writer.write((point[0] - 93.28850153024646) + " " + (point[1] - 528.8331031875234) + " "
					+ (point[2] - 96.61733426106944) + " " + (point[3] - 545.504002162989) + " "
					+ (length - 16.99999999));
		writer.write("\r\n");
		

		writer.close();

	} catch (IOException e) {
		e.printStackTrace();
	}
		
		
		}	
		
}
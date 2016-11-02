package peakFitter;

import java.util.ArrayList;

import edu.mines.jtk.sgl.BoundingBox;
import houghandWatershed.Boundingboxes;
import ij.gui.EllipseRoi;
import labeledObjects.LabelledImg;
import labeledObjects.Simpleobject;
import net.imglib2.Cursor;
import net.imglib2.Point;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.BenchmarkAlgorithm;
import net.imglib2.algorithm.OutputAlgorithm;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import peakFitter.GaussianMaskFitMSER.EndfitMSER;
import preProcessing.GetLocalmaxmin;

public class SubpixelLengthMSER extends BenchmarkAlgorithm
implements OutputAlgorithm<ArrayList<double[]>> {
	
	private static final String BASE_ERROR_MSG = "[SubpixelLineMSER] ";
	private final RandomAccessibleInterval<FloatType> source;
	private final ArrayList<LabelledImg> imgs;
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
	
	
	public SubpixelLengthMSER( final RandomAccessibleInterval<FloatType> source, 
			             final ArrayList<LabelledImg> imgs,
			             final ArrayList<Simpleobject> simpleobject,
			             final double[] psf,
			             final double minlength){
		
		this.source = source;
		this.imgs = imgs;
		this.simpleobject = simpleobject;
		this.psf = psf;
		this.minlength = minlength;
		this.ndims = source.numDimensions();
		
	}

	@Override
	public boolean checkInput() {
		if (source.numDimensions() > 2) {
			errorMessage = BASE_ERROR_MSG + " Can only operate on 1D, 2D, make slices of your stack . Got "
					+ source.numDimensions() + "D.";
			return false;
		}
		return true;
	}

	@Override
	public boolean process() {
		
		final_paramlist = new ArrayList<double[]>();
		for (int index = 0; index < simpleobject.size() ; ++index) {
			
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

		RandomAccessibleInterval<FloatType> currentimg = Boundingboxes.CurrentLabelImage(imgs, label);
				//imgs.get(label).Actualroiimg;
		EllipseRoi roi = imgs.get(label).roi;


		final Cursor<FloatType> inputcursor = Views.iterable(currentimg).localizingCursor();

		final double maxintensityline = GetLocalmaxmin.computeMaxIntensity(currentimg);

           while(inputcursor.hasNext()){
			
			inputcursor.fwd();
			
			final int x = inputcursor.getIntPosition(0);
			final int y = inputcursor.getIntPosition(1);
			
			if (roi.contains(x, y)){
				if (inputcursor.get().get() / maxintensityline > Intensityratio){
				inputcursor.localize(newposition);
				long pointonline = (long) (newposition[1] - slope * newposition[0] - intercept);

				// To get the min and max co-rodinates along the line so we have
				// starting points to
				// move on the line smoothly

				if (pointonline == 0) {

					for (int d = 0; d < ndims; ++d) {
						if (inputcursor.getDoublePosition(d) <= minVal[d])
							minVal[d] = inputcursor.getDoublePosition(d);

						if (inputcursor.getDoublePosition(d) >= maxVal[d])
							maxVal[d] = inputcursor.getDoublePosition(d);

					}

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
		
		
		System.out.println("Label: " + label + " " + "MSER Detection: " + " StartX: " + MinandMax[0] + " StartY: "
				+ MinandMax[1] + " EndX: " + MinandMax[2] + " EndY: " + MinandMax[3]);

		
		
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

				RandomAccessibleInterval<FloatType> currentimg = Boundingboxes.CurrentLabelImage(imgs, label);
						//imgs.get(label).Actualroiimg;

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
						LevenbergMarquardtSolverLine.solve(X, finalparamstart, fixed_param, I, new GaussianLineds(), lambda,
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
					final EllipseRoi roi = imgs.get(label).roi;


					System.out.println("Doing Mask Fits: ");
					try {
								
							startfit =	peakFitter.GaussianMaskFitMSER.sumofgaussianMaskFit(currentimg, roi, startpos.clone(), psf,
								iterations, dxvector, newslope, newintercept, maxintensityline, halfgaussian, EndfitMSER.Start,
								label);
					} catch (Exception e) {
						e.printStackTrace();
					}

					try {
						endfit = peakFitter.GaussianMaskFitMSER.sumofgaussianMaskFit(currentimg,roi, endpos.clone(), psf,
								iterations, dxvector, newslope, newintercept, maxintensityline,  halfgaussian, EndfitMSER.End,
								label);
					} catch (Exception e) {
						e.printStackTrace();
					}

					final double Maskdist = sqDistance(startfit, endfit);
					// If mask fits fail, return LM solver results, very crucial for
					// noisy data

					double[] returnparam = new double[2 * ndims + 3];
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

					returnparam[2 * ndims] = finalparamstart[4];
					returnparam[2 * ndims + 1] = finalparamstart[5];
					returnparam[2 * ndims + 2] = finalparamstart[6];
					return returnparam;

				}

				else
					return null;

			}
		}
		public int Getlabel(final Point linepoint) {

			final int x = linepoint.getIntPosition(0);
			final int y = linepoint.getIntPosition(1);
			int currentlabel = Integer.MIN_VALUE;
			for (int index = 0; index < imgs.size(); ++index){
				
				EllipseRoi ellipse = imgs.get(index).roi;
				
				if (ellipse.contains(x, y)){
					
					currentlabel = index;
					break;
				}
				
			}
			
		

			return currentlabel;
		}


		

		private PointSampleList<FloatType> gatherfullData(final int label) {
			final PointSampleList<FloatType> datalist = new PointSampleList<FloatType>(ndims);

			RandomAccessibleInterval<FloatType> currentimg = imgs.get(label).Actualroiimg;

			EllipseRoi roi = imgs.get(label).roi;
			
			Cursor<FloatType> localcursor = Views.iterable(currentimg).localizingCursor();

			while (localcursor.hasNext()) {
				localcursor.fwd();
				int x = localcursor.getIntPosition(0);
				int y = localcursor.getIntPosition(1);
				if (roi.contains(x, y)){
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

		public double sqDistance(final double[] cordone, final double[] cordtwo) {

			double distance = 0;

			for (int d = 0; d < ndims; ++d) {

				distance += Math.pow((cordone[d] - cordtwo[d]), 2);

			}
			return (distance);
		}

		
}

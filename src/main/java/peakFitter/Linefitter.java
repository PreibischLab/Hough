package peakFitter;

import houghandWatershed.Finalfunction;
import net.imglib2.Cursor;
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
	private final  double[] MakeLineguess(
			double slope, 
			double intercept,
			int label) throws Exception {

		
		final double[] realpos = new double[ndims];
		double sigmasq, sigma = 1.0;
		sigmasq = sigma * sigma;
		RandomAccessibleInterval<FloatType> imgout = new ArrayImgFactory<FloatType>().create(inputimg,
				new FloatType());;
		final Cursor<FloatType> outcursor = Views.iterable(imgout).localizingCursor();
		long[] newposition = new long[ndims];
		RandomAccess<IntType> ranac = intimg.randomAccess();
		final RandomAccess<FloatType> ranacinput = inputimg.randomAccess();
		double[] minVal = { Double.MAX_VALUE, Double.MAX_VALUE };
		double[] maxVal = { Double.MIN_VALUE, Double.MIN_VALUE };
		final double maxintensity =	 GetLocalmaxmin.computeMaxIntensityinlabel(inputimg,intimg,label );
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
				if (ranacinput.get().get()/maxintensity >= 0.5){
				outcursor.localize(newposition);

				long pointonline = (long) (newposition[1] - slope * newposition[0] - intercept);

				// To get the min and max co-rodinates along the line so we have starting points to
				// move on the line smoothly
				
				if (pointonline == 0 ) {
					
						
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
		
		final double[] MinandMax = new double[2*ndims + 1];
		
		if (slope >= 0){
		for (int d = 0; d < ndims; ++d) {
			
			MinandMax[d] = minVal[d];
			MinandMax[d + ndims] = maxVal[d];
		}
		
		}
		
		if (slope < 0){
			
			MinandMax[0] = minVal[0];
			MinandMax[1] = maxVal[1];
			MinandMax[2] = maxVal[0];
			MinandMax[3] = minVal[1];
			
		}
		
		MinandMax[2*ndims] = maxintensity;
		
		for (int d = 0; d < ndims; ++d) {
			
			if (MinandMax[d] == Double.MAX_VALUE || MinandMax[d + ndims] == Double.MIN_VALUE)
				 return null;
			
		}
		
		
		return MinandMax;
		
		
	
	}

	// Get line parameters for fitting line to a line in a label
	
	public double[] Getfinallineparam(final int label, final double slope, final double intercept, final double[] psf,
			final double [] sigma, double minlength) throws Exception{
		
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

		
			final double[] start_param = MakeLineguess(slope, intercept, label);
			if(start_param==null)
				return null;
			else{
		final double[] fixed_param = new double[2*ndims];
		
		for (int d = 0; d < ndims; ++d){
			
			fixed_param[d] = 1.0 / Math.pow(sigma[d] , 2);
		}
			fixed_param[ndims] = slope;
			fixed_param[ndims + 1] = intercept;
		final double[] finalparam = start_param.clone();
				
		
		// LM solver part
		int maxiter = 300;
		double lambda = 1e-3;
		double termepsilon = 1e-2;
		
			
		
		LevenbergMarquardtSolverLine.solve(X, finalparam, fixed_param, I, new GaussianLine(), lambda, termepsilon, maxiter);
		
		
		// NaN protection: we prefer returning the crude estimate than NaN
		for (int j = 0; j < finalparam.length; j++) {
			if (Double.isNaN(finalparam[j]))
				finalparam[j] = start_param[j];
		}
		
		final double[] startpos = {finalparam[0], finalparam[1]};
		final double[] endpos = {finalparam[2], finalparam[3]};
		System.out.println((endpos[1] -startpos[1])/(endpos[0]-startpos[0]) + " "
		+ (startpos[1] - ((endpos[1] -startpos[1])/(endpos[0]-startpos[0]))*startpos[0])+ " " + label);
        System.out.println(startpos[0] + " " +  startpos[1] + " " +  endpos[0] + " " + endpos[1]);
		double maxintensity = finalparam[4];
		double cutoffdistance = Distance(startpos, endpos);
		if (cutoffdistance > minlength ){
		
			int iterations = 100;
		
		double newslope = (finalparam[3] - finalparam[1]) / (finalparam[2] - finalparam[0]);
		double newintercept = finalparam[1] - slope * finalparam[0];
		final double radius = 0.25 * (psf[0] + psf[1]);
		final double[] bigsigma = {sigma[0]/2, sigma[1]/2};
		final double[] startfit = peakFitter.GaussianMaskFit.gaussianMaskFit(inputimg, intimg, startpos, sigma,
				iterations, maxintensity,radius,  newslope, newintercept, Endfit.Start);
		final double[] endfit =	peakFitter.GaussianMaskFit.gaussianMaskFit(inputimg, intimg, endpos, sigma,
				iterations, maxintensity,radius,  newslope, newintercept, Endfit.End);
		
		final double[] refindedparam = {startfit[0],startfit[1], 
				endfit[0], endfit[1], finalparam[4]};
		
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

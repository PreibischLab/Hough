package peakFitter;

import java.util.ArrayList;

import com.sun.tools.javac.util.Pair;

import drawandOverlay.AddGaussian;
import drawandOverlay.PushCurves;
import houghandWatershed.Finalfunction;
import ij.plugin.HyperStackReducer;
import labeledObjects.Finalobject;
import labeledObjects.Indexedlength;
import labeledObjects.LabelMax;
import labeledObjects.PreFinalobject;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.IterableInterval;
import net.imglib2.Localizable;
import net.imglib2.Point;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealCursor;
import net.imglib2.RealLocalizable;
import net.imglib2.RealPoint;
import net.imglib2.RealPointSampleList;
import net.imglib2.RealRandomAccess;
import net.imglib2.RealRandomAccessibleRealInterval;
import net.imglib2.algorithm.neighborhood.HyperSphereNeighborhood;
import net.imglib2.algorithm.neighborhood.RectangleShape;
import net.imglib2.algorithm.neighborhood.SquareStrelTest;
import net.imglib2.algorithm.region.hypersphere.HyperSphere;
import net.imglib2.algorithm.region.hypersphere.HyperSphereCursor;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.roi.RectangleRegionOfInterest;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Intervals;
import net.imglib2.view.Views;
import peakFitter.GaussianMastFit.Endfit;
import preProcessing.GetLocalmaxmin;
import preProcessing.GetLocalmaxmin.IntensityType;

public class LengthDetection {

	private final RandomAccessibleInterval<FloatType> inputimg;
	private final RandomAccessibleInterval<IntType> intimg;
	private final int ndims;

	public LengthDetection(RandomAccessibleInterval<FloatType> inputimg, RandomAccessibleInterval<IntType> intimg) {

		this.inputimg = inputimg;
		this.intimg = intimg;
		this.ndims = inputimg.numDimensions();
		assert inputimg.numDimensions() == intimg.numDimensions();

	}
	
	public static enum HalforFull{
		
		Half, Full;
		
	}

	
	
	
	private final double[] makeBestGuess(final Localizable point, final double[][] X, final double[] I) {

		double[] start_param = new double[2 * ndims + 2];

		double I_sum = 0;
		double[] X_sum = new double[ndims];
		for (int j = 0; j < ndims; j++) {
			X_sum[j] = 0;
			for (int i = 0; i < X.length; i++) {
				X_sum[j] += X[i][j] * I[i];
			}
		}
		double max_I = Double.NEGATIVE_INFINITY;
		for (int i = 0; i < X.length; i++) {
			I_sum += I[i];
			if (I[i] > max_I) {
				max_I = I[i];
			}

		}

		start_param[0] = max_I;

		for (int j = 0; j < ndims; j++) {
			start_param[j + 1] = X_sum[j] / I_sum;
		}

		for (int j = 0; j < ndims; j++) {
			double C = 0;
			double dx;
			for (int i = 0; i < X.length; i++) {
				dx = X[i][j] - start_param[j + 1];
				C += I[i] * dx * dx;
			}
			C /= I_sum;
			start_param[ndims + j + 1] = 1 / C;
		}
		start_param[2 * ndims + 1] = 0;

		return start_param;
	}

	private final double[] makeBestfixedsigmaGuess(final Localizable point, final double[][] X, final double[] I) {

		double[] start_param = new double[2 * ndims + 2];

		double I_sum = 0;
		double[] X_sum = new double[ndims];
		for (int j = 0; j < ndims; j++) {
			X_sum[j] = 0;
			for (int i = 0; i < X.length; i++) {
				X_sum[j] += X[i][j] * I[i];
			}
		}
		double max_I = Double.NEGATIVE_INFINITY;
		for (int i = 0; i < X.length; i++) {
			I_sum += I[i];
			if (I[i] > max_I) {
				max_I = I[i];
			}

		}

		start_param[0] = max_I;

		for (int j = 0; j < ndims; j++) {
			start_param[j + 1] = X_sum[j] / I_sum;
		}

		for (int j = 0; j < ndims; j++) {
			double C = 0;
			double dx;
			for (int i = 0; i < X.length; i++) {
				dx = X[i][j] - start_param[j + 1];
				C += I[i] * dx * dx;
			}
			C /= I_sum;
			start_param[ndims + j + 1] = 1 / C;
		}
		start_param[2 * ndims + 1] = 0;

		return start_param;
	}

	private final double[] makeNoiseGuess() {
		double[] start_param = new double[ndims - 1];

		start_param[ndims - 2] = 0.5;

		return start_param;
	}

	private final double[] makeBestpointsGuess(final Localizable point, final double[][] X, final double[] I) {

		double[] start_param = new double[2 * ndims + 2];

		double I_sum = 0;

		double[] X_sum = new double[ndims];
		for (int j = 0; j < ndims; j++) {
			X_sum[j] = 0;
			for (int i = 0; i < X.length; i++) {
				X_sum[j] += X[i][j] * I[i];
			}
		}
		double max_I = Double.NEGATIVE_INFINITY;
		for (int i = 0; i < X.length; i++) {
			I_sum += I[i];
			if (I[i] > max_I) {
				max_I = I[i];
			}

		}

		start_param[0] = max_I;

		for (int j = 0; j < ndims; j++) {
			start_param[j + 1] = X_sum[j] / I_sum;
		}

		for (int j = 0; j < ndims; j++) {
			double C = 0;
			double dx;
			for (int i = 0; i < X.length; i++) {
				dx = X[i][j] - start_param[j + 1];
				C += I[i] * dx * dx;
			}
			C /= I_sum;
			start_param[ndims + j + 1] = 1 / C;
		}

		start_param[2 * ndims + 1] = 0;

		return start_param;
	}
	
	private final  double[] MakeLineguess(
			double slope, 
			double intercept, final double[][] X, final double[] I,
			int label) {

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

			if (distance < 5 * sigma)
				intensity = (1 / (sigma * Math.sqrt(2 * Math.PI))) * Math.exp(-distance * distance / (2 * sigmasq));
			else
				intensity = 0;
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
		for (int d = 0; d < ndims; ++d) {
			
			MinandMax[d] = minVal[d];
			MinandMax[d + ndims] = maxVal[d];
		}
		
		MinandMax[2*ndims] = maxintensity;
		
		return MinandMax;
		
		}
	
	

	// Get line parameters for fitting line to a line in a label
	
	public double[] Getfinallineparam(final int label, final double slope, final double intercept, final double [] sigma) throws Exception{
		
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

		final double[] start_param = MakeLineguess(slope, intercept,X,I, label);
		
		System.out.println(start_param[0] + " " +  start_param[1] + " " +  start_param[2] + " " + start_param[3]);
		final double[] fixed_param = new double[2*ndims];
		
		for (int d = 0; d < ndims; ++d){
			
			fixed_param[d] = 1.0 / Math.pow(sigma[d] , 2);
		}
			fixed_param[ndims] = slope;
			fixed_param[ndims + 1] = intercept;
		final double[] finalparam = start_param.clone();
		
		// LM solver part
		int maxiter = 1000;
		double lambda = 1e-2;
		double termepsilon = 1e-3;

		LevenbergMarquardtSolverLine.solve(X, finalparam, fixed_param, I, new GaussianLine(), lambda, termepsilon, maxiter);

		// NaN protection: we prefer returning the crude estimate than NaN
		for (int j = 0; j < finalparam.length; j++) {
			if (Double.isNaN(finalparam[j]))
				finalparam[j] = start_param[j];
		}
		
		
		return finalparam;
		
	}
	
	
	// Get final parameters for the Gaussian fit along the line, input the
	// centroid position and value where the Gaussian has to be fitted
	// in the image
	public double[] Getfinalparam(final Localizable point, final long radius) throws Exception {

		PointSampleList<FloatType> datalist = gatherData(point, radius);

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

		final double[] start_param = makeBestGuess(point, X, I);

		final double[] finalparam = start_param.clone();
		int maxiter = 1000;
		double lambda = 1e-3;
		double termepsilon = 1e-3;

		LevenbergMarquardtSolverLocal.solve(X, finalparam, I, new GaussianMultiDLM(), lambda, termepsilon, maxiter);

		// NaN protection: we prefer returning the crude estimate than NaN
		for (int j = 0; j < finalparam.length; j++) {
			if (Double.isNaN(finalparam[j]))
				finalparam[j] = start_param[j];
		}
		
		

		return finalparam;

	}

	

	public double[] Getfinalpointsparam(final Localizable point, final long radius) throws Exception {

		PointSampleList<FloatType> datalist = gatherPointsData(point, radius);

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

		final double[] start_param = makeBestpointsGuess(point, X, I);

		final double[] finalparam = start_param.clone();
		int maxiter = 1000;
		double lambda = 1e-3;
		double termepsilon = 1e-3;

		LevenbergMarquardtSolverLocal.solve(X, finalparam, I, new GaussianandConstNoise(), lambda, termepsilon,
				maxiter);

		// NaN protection: we prefer returning the crude estimate than NaN
		for (int j = 0; j < finalparam.length; j++) {
			if (Double.isNaN(finalparam[j]))
				finalparam[j] = start_param[j];
		}

		return finalparam;

	}

	public double[] Getnoiseparam(final Localizable point, long radius) throws Exception {

		PointSampleList<FloatType> datalist = gatherData(point, radius);

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

		final double[] start_param = makeNoiseGuess();

		final double[] finalparam = start_param.clone();
		int maxiter = 1000;
		double lambda = 1e-3;
		double termepsilon = 1e-3;

		LevenbergMarquardtSolverLocal.solve(X, finalparam, I, new Nonoisepoiss(), lambda, termepsilon, maxiter);

		// NaN protection: we prefer returning the crude estimate than NaN
		for (int j = 0; j < finalparam.length; j++) {
			if (Double.isNaN(finalparam[j]))
				finalparam[j] = start_param[j];
		}

		return finalparam;

	}

	private PointSampleList<FloatType> gatherData(final Localizable point, final long radius) {
		final PointSampleList<FloatType> datalist = new PointSampleList<FloatType>(ndims);

		for (int d = 0; d < ndims; ++d)
			assert inputimg.dimension(d) == intimg.dimension(d);

		
		RandomAccess<FloatType> ranac = inputimg.randomAccess();

		ranac.setPosition(point);

		RandomAccess<IntType> intranac = intimg.randomAccess();

		
		intranac.setPosition(point);

		final double[] position = new double[ndims];
		point.localize(position);

		// Gather data around the point
		boolean outofbounds = false;
		for (int d = 0; d < ndims; d++) {
			
			if (point.getDoublePosition(d) <= 0 || point.getDoublePosition(d)>= inputimg.dimension(d)){
				
				outofbounds = true;
				break;
			}
		}
		
		HyperSphere<FloatType> region = new HyperSphere<FloatType>(inputimg, point, radius);

		HyperSphereCursor<FloatType> localcursor = region.localizingCursor();

		final int label = intranac.get().get();

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

	
	private PointSampleList<FloatType> gatherPointsData(final Localizable point, final long radius) {
		final PointSampleList<FloatType> datalist = new PointSampleList<FloatType>(ndims);

		RandomAccess<IntType> intranac = intimg.randomAccess();
		final double[] position = new double[ndims];
		point.localize(position);
		intranac.setPosition(point);
		final int label = intranac.get().get();
		// Gather data around the point
		HyperSphere<FloatType> region = new HyperSphere<FloatType>(inputimg, point, radius);

		HyperSphereCursor<FloatType> localcursor = region.localizingCursor();

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

	public ArrayList<Finalobject> Updateslopeandintercept(ArrayList<Finalobject> finalparam) {

		int Maxlabel = houghandWatershed.PerformWatershedding.GetMaxlabelsseeded(intimg);
		final ArrayList<Finalobject> updateparamlist = new ArrayList<Finalobject>();
		for (int label = 1; label < Maxlabel - 1; ++label) {

			ArrayList<RealPoint> centroidlist = new ArrayList<RealPoint>();

			int labelindex = 0;

			for (int index = 0; index < finalparam.size(); ++index) {

				if (finalparam.get(index).Label == label) {

					labelindex = index;
					centroidlist.add(finalparam.get(index).centroid);

				}

			}

			double newslope = 0;
			double newintercept = 0;

			final double[] pointone = new double[ndims];
			final double[] pointtwo = new double[ndims];

			if (centroidlist.size() > 0) {
				centroidlist.get(0).localize(pointone);
				centroidlist.get(centroidlist.size() - 1).localize(pointtwo);

				if (pointtwo[0] != pointone[0] ){
				newslope = (pointtwo[1] - pointone[1]) / (pointtwo[0] - pointone[0]);
				newintercept = pointtwo[1] - newslope * pointtwo[0];
				}
				
				else{
					centroidlist.get(centroidlist.size()/2).localize(pointone);
					centroidlist.get(3*centroidlist.size()/4).localize(pointtwo);
					newslope = (pointtwo[1] - pointone[1]) / (pointtwo[0] - pointone[0]);
					newintercept = pointtwo[1] - newslope * pointtwo[0];
				}
				

				if (Double.isNaN(newslope) ){
					newslope = finalparam.get(labelindex).slope ;
					newintercept = finalparam.get(labelindex).intercept;
				}
				
					
				
			}

			System.out.println("old: " + finalparam.get(labelindex).slope + " " + finalparam.get(labelindex).intercept);

			System.out.println("new: " + newslope + " " + newintercept);

			for (int index = 0; index < finalparam.size(); ++index) {

				if (finalparam.get(index).Label == label) {
					Finalobject update = new Finalobject(label, finalparam.get(index).centroid,
							finalparam.get(index).Intensity, finalparam.get(index).sigmaX, finalparam.get(index).sigmaY,
							newslope, newintercept);

					updateparamlist.add(update);

				}
			}

		}

		return updateparamlist;
	}

	public ArrayList<Finalobject> Removepoints(ArrayList<Finalobject> finalparam, ArrayList<LabelMax> labelmaxlist) {

		int Maxlabel = houghandWatershed.PerformWatershedding.GetMaxlabelsseeded(intimg);

		final ArrayList<Finalobject> correctparamlist = new ArrayList<Finalobject>();
		for (int label = 1; label < Maxlabel - 1; ++label) {
			
			float[] pos = new float[ndims];
			double maxintensity = Double.MIN_VALUE;

			for (int index = 0; index < finalparam.size(); ++index) {

				if (finalparam.get(index).Label == label) {

					
					if (finalparam.get(index).Intensity > maxintensity) {

						maxintensity = finalparam.get(index).Intensity;

					}

				}

			}

			final LabelMax labelmax = new LabelMax(label, maxintensity);
			labelmaxlist.add(labelmax);
			System.out.println("Label :" + label + " " + maxintensity);

			for (int listindex = 0; listindex < finalparam.size(); ++listindex) {

				if (finalparam.get(listindex).Label == label) {

					if (finalparam.get(listindex).Intensity / maxintensity > 0.5   )

						correctparamlist.add(finalparam.get(listindex));
				}

			}
		}

		return correctparamlist;

	}

	public void Returnlengths(ArrayList<Finalobject> finalparam, ArrayList<Indexedlength> finallength,
			ArrayList<LabelMax> labelmaxlist, double[] sigma) throws Exception {

		int Maxlabel = houghandWatershed.PerformWatershedding.GetMaxlabelsseeded(intimg);

		for (int label = 1; label < Maxlabel - 1; ++label) {
			final double[] minVal = { Double.MAX_VALUE, Double.MAX_VALUE };
			final double[] maxVal = { Double.MIN_VALUE, Double.MIN_VALUE };
			double length = 0;
			double lengthpre = 0;

			double[] minposition = new double[ndims];
			double[] maxposition = new double[ndims];

			double[] startpos = new double[ndims];
			double[] endpos = new double[ndims];
			double slope = 0;
			double intercept = 0;
			int labelindex = 0;
			for (int index = 0; index < finalparam.size(); ++index) {

				if (finalparam.get(index).Label == label) {
					labelindex = index;

					for (int d = 0; d < ndims; ++d) {
						minposition[d] = finalparam.get(index).centroid.getDoublePosition(d);
						maxposition[d] = finalparam.get(index).centroid.getDoublePosition(d);

						
						if (minposition[d] <= minVal[d]) {
							minVal[d] = minposition[d];
						}
						if (maxposition[d] >= maxVal[d]) {
							maxVal[d] = maxposition[d];
						}

					}

					if (finalparam.get(index).slope >= 0) {

						for (int d = 0; d < ndims; ++d) {

							startpos[d] = minVal[d];
							endpos[d] = maxVal[d];

						}
					}

					if (finalparam.get(index).slope < 0) {

						startpos[0] = minVal[0];
						startpos[1] = maxVal[1];

						endpos[0] = maxVal[0];
						endpos[1] = minVal[1];

					}

				}

			}

			slope = finalparam.get(labelindex).slope;
			intercept = finalparam.get(labelindex).intercept;
			lengthpre = Distance(startpos, endpos);

			System.out.println("Label :" + label + " " + " StartX: " + startpos[0] + " " + " StartY: " + startpos[1]
					+ "EndX :" + endpos[0] + " " + "EndY :" + endpos[1] + " " + "Length: " + lengthpre);

			final int iterations = 5000;
			double maxintensity = 0;

			for (int listlabelindex = 0; listlabelindex < labelmaxlist.size(); ++listlabelindex) {

				if (labelmaxlist.get(listlabelindex).Label == label)

					maxintensity = labelmaxlist.get(listlabelindex).maxIntensity;

			}

			final double[] newsigma = { sigma[0], sigma[1] };
			final double[] startfit = peakFitter.GaussianMastFit.gaussianMaskFit(inputimg, intimg, startpos, newsigma,
					iterations, maxintensity, 1.0, slope, intercept, Endfit.Start);
			final double[] endfit = peakFitter.GaussianMastFit.gaussianMaskFit(inputimg, intimg, endpos, newsigma,
					iterations, maxintensity, 1.0, slope, intercept, Endfit.End);

			
			

			length = Distance(startfit, endfit);
			final Indexedlength currentlength = new Indexedlength(label, length, startfit, endfit, slope, intercept);

			finallength.add(currentlength);

			System.out.println(
					"New:" + "Label :" + label + " " + " StartX: " + startfit[0] + " " + " StartY: " + startfit[1] + " "
							+ "EndX :" + endfit[0] + " " + "EndY :" + endfit[1] + "  " + "Length: " + length);

		}

	}

}

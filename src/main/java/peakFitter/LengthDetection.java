package peakFitter;

import drawandOverlay.AddGaussian;
import ij.plugin.HyperStackReducer;
import net.imglib2.Cursor;
import net.imglib2.IterableInterval;
import net.imglib2.Localizable;
import net.imglib2.Point;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealLocalizable;
import net.imglib2.algorithm.neighborhood.HyperSphereNeighborhood;
import net.imglib2.algorithm.neighborhood.RectangleShape;
import net.imglib2.algorithm.region.hypersphere.HyperSphere;
import net.imglib2.algorithm.region.hypersphere.HyperSphereCursor;
import net.imglib2.roi.RectangleRegionOfInterest;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;

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

	
	private final double[] makeBestGuess(final Localizable point, final double[][] X, final double[] I, final double [] psf) {

		double[] start_param = new double[2*ndims+2];

		
		double I_sum = 0;
		double max_I = Double.NEGATIVE_INFINITY;
		for (int i = 0; i < X.length; i++) {
			I_sum += I[i];
			if (I[i] > max_I) {
				max_I = I[i];
			}

		}

		start_param[0] = max_I;

		
		for (int d = 0; d < ndims; ++d) {

			start_param[d + 1] = point.getDoublePosition(d);
		}

		for (int j = 0; j < ndims; j++) {
			
			start_param[ndims + j + 1] = 1.0/Math.pow(psf[j],2);
		}
		
		
		
		start_param[2*ndims+1] = 0; 
		
		return start_param;		
	}

	private final double[] makeNoiseGuess(){
		double[] start_param = new double[ndims];
		
		for (int d = 0; d < ndims; ++d)
			start_param[d] = 0.5;
		
		return start_param;
	}
	private final double[] makeBestpointsGuess(final Localizable point, final double[][] X, final double[] I) {

		double[] start_param = new double[2*ndims+2];

		
		double I_sum = 0;
		double max_I = Double.NEGATIVE_INFINITY;
		for (int i = 0; i < X.length; i++) {
			I_sum += I[i];
			if (I[i] > max_I) {
				max_I = I[i];
			}

		}

		start_param[0] = max_I;

		
		for (int d = 0; d < ndims; ++d) {

			start_param[d + 1] = point.getDoublePosition(d);
		}

		for (int j = 0; j < ndims; j++) {
			double C = 0;
			double dx;
			for (int i = 0; i < X.length; i++) {
				dx = X[i][j] - start_param[j+1];
				C += I[i] * dx * dx;
			}
			C /= I_sum;
			start_param[ndims + j + 1] = 1 / C;
		}
		
		start_param[2*ndims+1] = 0; 
		
		return start_param;		
	}
	
	
	
	// Get final parameters for the Gaussian fit along the line, input the
	// centroid position and value where the Gaussian has to be fitted
	// in the image
	public double[] Getfinalparam(final Localizable point, final long radius, final double[] psf)
			throws Exception {



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

		final double[] start_param = makeBestGuess(point, X, I, psf);
		
		final double[] finalparam = start_param.clone();
		int maxiter = 1000;
		double lambda = 1e-3;
		double termepsilon = 1e-1;

		LevenbergMarquardtSolverLocal.solve(X, finalparam, I, new GaussianMultiDLM(), lambda, termepsilon, maxiter);
		
		// NaN protection: we prefer returning the crude estimate than NaN
		for (int j = 0; j < finalparam.length; j++) {
			if (Double.isNaN(finalparam[j])  )
				finalparam[j] = start_param[j];
		}
	
	
		
		return finalparam;

	}
	
	
	public double[] Getfinalpointsparam(final Localizable point, final long radius)
			throws Exception {



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
		double termepsilon = 1e-1;

		LevenbergMarquardtSolverLocal.solve(X, finalparam, I, new GaussianandConstNoise(), lambda, termepsilon, maxiter);
	
		// NaN protection: we prefer returning the crude estimate than NaN
		for (int j = 0; j < finalparam.length; j++) {
			if (Double.isNaN(finalparam[j])  )
				finalparam[j] = start_param[j];
		}
		
	
		
		return finalparam;

	}

    public double[] Getnoiseparam(final Localizable point, long radius, double[] psf) 
	throws Exception {



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
		double termepsilon = 1e-1;

		LevenbergMarquardtSolverLocal.solve(X, finalparam, I, new Nonoisepoiss(), lambda, termepsilon, maxiter);
	
		// NaN protection: we prefer returning the crude estimate than NaN
		for (int j = 0; j < finalparam.length; j++) {
			if (Double.isNaN(finalparam[j])  )
				finalparam[j] = start_param[j];
		}
		
	
		
		return finalparam;

	}
	
	
	private PointSampleList<FloatType> gatherData(final Localizable point, final long radius){
		final PointSampleList<FloatType> datalist = new PointSampleList<FloatType>(ndims);
		
		RandomAccess<FloatType> ranac = inputimg.randomAccess();
		
		ranac.setPosition(point);
		final double[] position = new double[ndims];
		point.localize(position);
		
		
		// Gather data around the point
		
		
		
        HyperSphere<FloatType> region = new HyperSphere<FloatType>(inputimg, point, radius);
		
		HyperSphereCursor<FloatType> localcursor = region.localizingCursor();
		  RandomAccess<IntType> intranac = intimg.randomAccess();
		  intranac.setPosition(point);
		 final int label = intranac.get().get();
			boolean outofbounds = false;
			  while(localcursor.hasNext()){
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
				  if (i == label){
				  
				    Point newpoint = new Point(localcursor);
					datalist.add(newpoint, localcursor.get().copy());
					
					
				  }
			  }
			
				
		return datalist;
	}
	
	private PointSampleList<FloatType> gatherPointsData(final Localizable point, final long radius){
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
	  while(localcursor.hasNext()){
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
		  if (i == label){
		  
		    Point newpoint = new Point(localcursor);
			datalist.add(newpoint, localcursor.get().copy());
			
			
		  }
	  }
		
		return datalist;
	}


	
	
	
	
}

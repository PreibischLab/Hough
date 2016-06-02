package varun;

import java.util.ArrayList;

import com.sun.tools.javac.util.Pair;

import net.imglib2.Cursor;
import net.imglib2.Interval;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.localization.LevenbergMarquardtSolver;
import net.imglib2.algorithm.neighborhood.Neighborhood;
import net.imglib2.algorithm.neighborhood.RectangleShape;
import net.imglib2.algorithm.region.hypersphere.HyperSphere;
import net.imglib2.algorithm.region.hypersphere.HyperSphereCursor;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Intervals;
import net.imglib2.view.Views;
import varun.OverlayLines.Simulatedline;
import varun.PerformWatershedding.Lineobjects;

public class LengthDetection {
	public static final class Labelparam {
		final int Label;
		final double[] point;
		final double slope;

		protected Labelparam(final int Label, final double[] point, final double slope) {
			this.Label = Label;
			this.point = point;
			this.slope = slope;

		}
	}
	
	/*** Make the Error of Match function b/w the model and the data ***/
	
	
	
	public static void ErrorofMatch(
			RandomAccessibleInterval<FloatType> ModelImage,
			RandomAccessibleInterval<FloatType> DataImage,
			RandomAccessibleInterval<IntType> IntImage,
			int label){
		int n = ModelImage.numDimensions();
		Cursor<IntType> intCursor = Views.iterable(IntImage).localizingCursor();
		RandomAccess<FloatType> modelranac = ModelImage.randomAccess();
		RandomAccess<FloatType> dataranac = DataImage.randomAccess();
		ArrayList<Simulatedline> listdata = new ArrayList<Simulatedline>();
		ArrayList<Simulatedline> listmodel = new ArrayList<Simulatedline>();
		while(intCursor.hasNext()){
			intCursor.fwd();
			
			int i = intCursor.get().get();
			if (i == label) {
				
				final double[] modelposition = new double[n];
				final double[] dataposition = new double[n];
				modelranac.setPosition(intCursor);
				dataranac.setPosition(intCursor);
				modelranac.localize(modelposition);
				dataranac.localize(dataposition);
				final Simulatedline dataline = new Simulatedline(label, modelposition, modelranac.get());
				listdata.add(dataline);
				
				final Simulatedline modeline = new Simulatedline(label, dataposition, dataranac.get());
				listmodel.add(modeline);
				
			}

			
		}
		
		
		
	}
	
	

	
	/**
	Input a guess list containing the list of centroids determined by the HT, the psf of the microscope and 
	@returns a doube array startparam[] containing guess for the Gaussian at one centroid point
	 
	 *  start_param[0] = A
	 * start_param[1 → ndims] = x_i, y_i (by HT)
	 * start_param[ndims+1 → 2 × ndims]  = 1 / σ_i2
	 *        
	 */
	
	
	public static double[] makeBestGuess(
			final ArrayList<Simulatedline> guessline, 
			final double[] psf,
			int listindex,
			int ndims, 
			int label) {
		double[][] X = new double[guessline.size()][ndims];
		double[] I = new double[guessline.size()];
	
			for (int d = 0; d < ndims; ++d) {
                 if (guessline.get(listindex).Label == label)
				X[listindex][d] = guessline.get(listindex).point[d];

			}

		
			 if (guessline.get(listindex).Label == label)
			I[listindex] = guessline.get(listindex).Value.get();
		

		

		double[] start_param = new double[2 * ndims + 1];

		start_param[0] = I[listindex];

		for (int j = 0; j < ndims; j++) {
			start_param[j + 1] = X[listindex][j];
		}

		for (int j = 0; j < ndims; j++) {
			
			start_param[ndims + j + 1] = 1 / psf[j];
		}
		
		return start_param;
	}
	
	
	public static void makeBeads(
			 RandomAccessibleInterval<FloatType> inputimg,
			RandomAccessibleInterval<IntType> intimg,
			int label){
		
		int n = intimg.numDimensions();
		 Interval interval = Intervals.expand( inputimg, -1 );
		 
	        inputimg = Views.interval( inputimg, interval );
	        
	        Interval intervalint = Intervals.expand( intimg, -1 );
			 
	        intimg = Views.interval( intimg, intervalint );

	        final Cursor< IntType > center = Views.iterable( intimg ).cursor();
	        
	       ArrayList<Simulatedline> alllines = new ArrayList<Simulatedline>();
	        final RectangleShape shape = new RectangleShape( 1, true );
	        for ( final Neighborhood< IntType > localNeighborhood : shape.neighborhoods( intimg ) )
	        {
	            
	            final IntType centerValue = center.next();
	 
	            // keep this boolean true as long as no other value in the local neighborhood
	            // is larger or equal
	            boolean isinLabel = true;
	 
	            // check if all pixels in the local neighborhood that are smaller
	            for ( final IntType value : localNeighborhood )
	            {
	                // test if the center is smaller than the current pixel value
	                if ( centerValue.get()!= value.get() )
	                {
	                    isinLabel = false;
	                    break;
	                }
	            }
	 
	            if ( isinLabel )
	            {
	            	
	            	int i = centerValue.get();
	            	if (i == label){
	                // draw a sphere of radius one in the new image
	                HyperSphere< IntType > hyperSphere = new HyperSphere< IntType >( intimg, center, 1 );
	                HyperSphereCursor<IntType> hypercursor = hyperSphere.localizingCursor();
	                RandomAccess<FloatType> ranac = inputimg.randomAccess();

	                while(hypercursor.hasNext()){
	                	hypercursor.fwd();
	                	double[] position = new double[n];
	                	hypercursor.localize(position);
	                	ranac.setPosition(hypercursor);
	                Simulatedline simline = new Simulatedline(label, position, ranac.get());

	                alllines.add(simline);
	                
	                }
	                
	                
	            	}
	            }
	        }

		
		
	}
	
	
	public static ArrayList<double[]> makeBeads(
			final RandomAccessibleInterval<FloatType> inputimg,
			Img<IntType> intimg,
			final ArrayList<Simulatedline> linelist, 
			double[] typical_sigma,
			int ndims, 
			int label, 
			int setlength) {
		
		
		ArrayList<double[]> beadlist = new ArrayList<double[]>();
		
	
			
			double[][] X = new double[linelist.size()][ndims];
			double[] I = new double[linelist.size()];
		
			Pair<double[][], double[]> MLEset = makeMLEsubsets(linelist, ndims, label, 0, linelist.size()-1);
		
		X = MLEset.fst;

		I = MLEset.snd;
		
		
		double[] start_param = new double[2 * ndims + 1];
		double[] finalparam = new double[2 * ndims + 1];
		double[] X_sum = new double[ndims];
		for (int j = 0; j < ndims; j++) {
			X_sum[j] = 0;
			for (int i = 0; i < X.length; i++) {
				X_sum[j] += X[i][j] * I[i];
			}
		}

		double I_sum = 0;
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
		
		finalparam = getExactparam(inputimg, intimg,start_param,typical_sigma,label);
		
		beadlist.add(finalparam);
		
	
		return beadlist;
	}
	

	
	public static double[] getExactparam(
			final RandomAccessibleInterval<FloatType> inputimg, 
			Img<IntType> intimg,
			double[] start_param,
			double[] typical_sigma,
			int label){
		    int n = inputimg.numDimensions();
		
		Cursor<IntType> intCursor = intimg.localizingCursor();
		RandomAccess<FloatType> ranac = inputimg.randomAccess();
		ArrayList<Simulatedline> listdata = new ArrayList<Simulatedline>();
		
		while(intCursor.hasNext()){
			intCursor.fwd();
			
			int i = intCursor.get().get();
			if (i == label) {
				
				final double[] position = new double[n];
				ranac.setPosition(intCursor);
				ranac.localize(position);
				final Simulatedline line = new Simulatedline(label, position, ranac.get());
				listdata.add(line);
				
			}

			
		}
		
		double[][] X = new double[listdata.size()][n];
		double[] I = new double[listdata.size()];
		
		for (int index = 0; index < listdata.size(); ++index) {
			for (int d = 0; d < n; ++d) {
                 if (listdata.get(index).Label == label)
				X[index][d] = listdata.get(index).point[d];

			}

		}
		for (int index = 0; index < listdata.size(); ++index) {
			 if (listdata.get(index).Label == label)
			I[index] = listdata.get(index).Value.get();
		}
		
		
		
		int maxiter = 500;
		double lambda = 1e-2;
		double termepsilon = 1e-1;
		
		final double[] finalparam = start_param.clone();
		
		try {
			LevenbergMarquardtSolverLocal.solve(X, finalparam, I, new GaussianMultiDLM(), lambda , termepsilon, maxiter);
		} catch (Exception e) {
			e.printStackTrace();
		} 
		
		
		
		return finalparam;
	}
	public static Pair<double[][], double[]> makeMLEsets(final ArrayList<Simulatedline> linelist, int ndims, int label){
		double[][] X = new double[linelist.size()][ndims];
		double[] I = new double[linelist.size()];
		
		for (int index = 0; index < linelist.size(); ++index) {
			for (int d = 0; d < ndims; ++d) {
                 if (linelist.get(index).Label == label && linelist.get(index).Value.get()>0)
				X[index][d] = linelist.get(index).point[d];

			}

		}
		for (int index = 0; index < linelist.size(); ++index) {
			 if (linelist.get(index).Label == label && linelist.get(index).Value.get()>0)
			   I[index] = linelist.get(index).Value.get();
		}
		
		Pair<double[][], double[]> MLEset = new Pair<double[][], double[]>(X,I);
		
		return MLEset;
		
	}
	
	public static Pair<double[][], double[]> makeMLEsubsets(
			final ArrayList<Simulatedline> linelist, 
			int ndims, 
			int label, 
			int startindex, 
			int endindex){
		double[][] X = new double[linelist.size()][ndims];
		double[] I = new double[linelist.size()];
		
		for (int index = startindex; index < endindex; ++index) {
			for (int d = 0; d < ndims; ++d) {
                 if (linelist.get(index).Label == label && linelist.get(index).Value.get()>0)
				X[index][d] = linelist.get(index).point[d];

			}

		}
		for (int index = startindex; index < endindex; ++index) {
			 if (linelist.get(index).Label == label && linelist.get(index).Value.get()>0)
			   I[index] = linelist.get(index).Value.get();
			
		}
		
		
        Pair<double[][], double[]> MLEset = new Pair<double[][], double[]>(X,I);
		
		return MLEset;
		
		
	}
}

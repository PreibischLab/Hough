package peakFitter;

import net.imglib2.Cursor;
import net.imglib2.IterableInterval;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
public class GaussianMaskFit {
	public static enum Endfit {
		Start, End
	}

	protected Endfit Endfit;

	public static double[] gaussianMaskFit(final RandomAccessibleInterval<FloatType> signalInterval,
			final RandomAccessibleInterval<IntType> intimg, final double[] location, final double[] sigma,
			final int iterations, final double maxintensity, final double deltas, final double slope, final double intercept,
			final Endfit startorend) {
		final int n = signalInterval.numDimensions();

		// pre-compute sigma^2
		final double[] sq_sigma = new double[n];
		for (int d = 0; d < n; ++d)
			sq_sigma[d] =  sigma[d] * sigma[d];

		// make the interval we fit on iterable
		final IterableInterval<FloatType> signalIterable = Views.iterable(signalInterval);

		// create the mask image
		final Img<FloatType> gaussianMask = new ArrayImgFactory<FloatType>().create(signalInterval,
				signalIterable.firstElement());

		// set the mask image to the same location as the interval we fit on and
		// make it iterable
		final long[] translation = new long[n];
		for (int d = 0; d < n; ++d)
			translation[d] = signalInterval.min(d);

		final RandomAccessibleInterval<FloatType> translatedMask = Views.translate(gaussianMask, translation);
		final IterableInterval<FloatType> translatedIterableMask = Views.iterable(translatedMask);

		// remove background in the input
		final double bg = removeBackground(signalIterable);

		double N = 0;
		int i = 0;
		do {
			
			switch (startorend) {

			case Start:
				setstartGaussian(translatedIterableMask, location, sq_sigma, maxintensity,deltas, slope, intercept);
				break;

			case End:
				setendGaussian(translatedIterableMask, location, sq_sigma, maxintensity,deltas,  slope, intercept);
				break;

			}
			/*
			final long[] longintlocation = new long[n];
			for (int d = 0; d <n; ++d){
				longintlocation[d] = (long) location[d];
				
			}
			
            final RandomAccess<IntType> intranac = intimg.randomAccess();
			
             boolean outofbounds = false;
            for (int d = 0; d <n ; ++d ){
			
            	if (longintlocation[d] <= 0 || longintlocation[d] >= intimg.dimension(d)){
            	outofbounds = true;
            		break;
            	}
            	else
            	intranac.setPosition(longintlocation);
            }
            if (outofbounds == false){
            final int label = intranac.get().get();
            
		
            */
			//ImageJFunctions.show(gaussianMask);
			// compute the sums
			final Cursor<FloatType> cMask = gaussianMask.cursor();
			final Cursor<FloatType> cImg = signalIterable.localizingCursor();

			double sumLocSN[] = new double[n]; // int_{all_px} d * S[ d ] * N[ d
												// ]
			double sumSN = 0; // int_{all_px} S[ d ] * N[ d ]
			double sumSS = 0; // int_{all_px} S[ d ] * S[ d ]

			while (cMask.hasNext()) {
				cMask.fwd();
				cImg.fwd();
			//	intranac.setPosition(cImg);
			//	if (intranac.get().get() == label){
					final double signal = cImg.get().getRealDouble();
					final double mask = cMask.get().getRealDouble();
					final double weight = 8;

					final double signalmask = signal * mask * weight;

					sumSN += signalmask;
					sumSS += signal * signal * weight;

					for (int d = 0; d < n; ++d) {
						final double l = cImg.getDoublePosition(d);
						sumLocSN[d] += l * signalmask;
					}

				}

				for (int d = 0; d < n; ++d)
					location[d] = sumLocSN[d] / sumSN;

				N = sumSN / sumSS;

				++i;
		//	}
		//	}
		//	else
         //   	break;
		
		} while (i < iterations);
		restoreBackground(signalIterable, bg);
		//ImageJFunctions.show(gaussianMask);
		return location;
	}

	public static double removeBackground(final IterableInterval<FloatType> iterable) {
		double i = 0;

		for (final FloatType t : iterable)
			i += t.getRealDouble();

		i /= (double) iterable.size();

		for (final FloatType t : iterable)
			t.setReal(t.get() - i);

		return i;
	}

	public static void restoreBackground(final IterableInterval<FloatType> iterable, final double value) {
		for (final FloatType t : iterable)
			t.setReal(t.get() + value);
	}

	final public static void setstartGaussian(final IterableInterval<FloatType> image, final double[] location,
			final double[] sq_sigma, final double maxintensity,final double deltas,  final double slope, final double intercept) {
		final int numDimensions = image.numDimensions();

		final Cursor<FloatType> cursor = image.localizingCursor();

		while (cursor.hasNext()) {
			cursor.fwd();

			double value = maxintensity;
			double constbox, secondconstbox, thirdconstbox;
		
			constbox = Math.exp(-deltas*deltas/((1+slope*slope)*sq_sigma[0]))*
					Math.exp(-slope*slope*deltas*deltas/((1+slope*slope)*sq_sigma[1]));
			secondconstbox = Math.exp(-4*deltas*deltas/((1+slope*slope)*sq_sigma[0]))*
					Math.exp(-slope*slope*4*deltas*deltas/((1+slope*slope)*sq_sigma[1]));
			thirdconstbox = Math.exp(-9*deltas*deltas/((1+slope*slope)*sq_sigma[0]))*
					Math.exp(-slope*slope*9*deltas*deltas/((1+slope*slope)*sq_sigma[1]));
		
			double totalbox = constbox;
			double secondtotalbox = secondconstbox;
			double thirdtotalbox = thirdconstbox;
			
			for (int d = 0; d < numDimensions; ++d) {
				final double x = location[d] - cursor.getDoublePosition(d);
				
				// Full Gaussian fit
			
				
			if (d == 0){
				if (slope >= 0){
				totalbox *= Math.exp(2*x*deltas/(sq_sigma[d]*Math.sqrt(1+slope*slope)));
				secondtotalbox *= Math.exp(2*2*x*deltas/(sq_sigma[d]*Math.sqrt(1+slope*slope)));
				thirdtotalbox *= Math.exp(2*3*x*deltas/(sq_sigma[d]*Math.sqrt(1+slope*slope)));
				}
				else{
				totalbox *= Math.exp(-2*x*deltas/(sq_sigma[d]*Math.sqrt(1+slope*slope)));
				secondtotalbox *= Math.exp(-2*2*x*deltas/(sq_sigma[d]*Math.sqrt(1+slope*slope)));
				thirdtotalbox *= Math.exp(-2*3*x*deltas/(sq_sigma[d]*Math.sqrt(1+slope*slope)));
				}
				
			}
			if (d == 1){
				if (slope >=0){
				totalbox *=Math.exp(2*slope*x*deltas/(sq_sigma[d]*Math.sqrt(1+slope*slope)));
				secondtotalbox *=Math.exp(2*2*slope*x*deltas/(sq_sigma[d]*Math.sqrt(1+slope*slope)));
				thirdtotalbox *=Math.exp(2*3*slope*x*deltas/(sq_sigma[d]*Math.sqrt(1+slope*slope)));
				}
				else{
				totalbox *=Math.exp(2*slope*x*deltas/(sq_sigma[d]*Math.sqrt(1+slope*slope)));
				secondtotalbox *=Math.exp(2*2*slope*x*deltas/(sq_sigma[d]*Math.sqrt(1+slope*slope)));
				thirdtotalbox *=Math.exp(2*3*slope*x*deltas/(sq_sigma[d]*Math.sqrt(1+slope*slope)));
				}
			}
				  
				value *= Math.exp(-(x * x) / sq_sigma[d]) * (1 + totalbox + secondtotalbox + thirdtotalbox )  ;
				
			}

			cursor.get().setReal(value);
		}
	}

	final public static void setendGaussian(final IterableInterval<FloatType> image, final double[] location,
			final double[] sq_sigma, final double maxintensity, final double deltas,  final double slope, final double intercept) {
		final int numDimensions = image.numDimensions();

		final Cursor<FloatType> cursor = image.localizingCursor();

		while (cursor.hasNext()) {
			cursor.fwd();

			double value = maxintensity;
			
			double constbox, secondconstbox, thirdconstbox;
			
			constbox = Math.exp(-deltas*deltas/((1+slope*slope)*sq_sigma[0]))*
					Math.exp(-slope*slope*deltas*deltas/((1+slope*slope)*sq_sigma[1]));
			secondconstbox = Math.exp(-4*deltas*deltas/((1+slope*slope)*sq_sigma[0]))*
					Math.exp(-slope*slope*4*deltas*deltas/((1+slope*slope)*sq_sigma[1]));
			thirdconstbox = Math.exp(-9*deltas*deltas/((1+slope*slope)*sq_sigma[0]))*
					Math.exp(-slope*slope*9*deltas*deltas/((1+slope*slope)*sq_sigma[1]));
			
			double totalbox = constbox;
			double secondtotalbox = secondconstbox;
			double thirdtotalbox = thirdconstbox;
			for (int d = 0; d < numDimensions; ++d) {
				final double x = location[d] - cursor.getDoublePosition(d);
				
				// Full Gaussian fit
				
		
				if (d == 0){
					if (slope >= 0){
					totalbox *=Math.exp(-2*x*deltas/(sq_sigma[d]*Math.sqrt(1+slope*slope)));
					secondtotalbox *=Math.exp(-2*2*x*deltas/(sq_sigma[d]*Math.sqrt(1+slope*slope)));
					thirdtotalbox *=Math.exp(-2*3*x*deltas/(sq_sigma[d]*Math.sqrt(1+slope*slope)));
					}
					else{
					totalbox *=Math.exp(2*x*deltas/(sq_sigma[d]*Math.sqrt(1+slope*slope)));	
					secondtotalbox *=Math.exp(2*2*x*deltas/(sq_sigma[d]*Math.sqrt(1+slope*slope)));	
					thirdtotalbox *=Math.exp(2*3*x*deltas/(sq_sigma[d]*Math.sqrt(1+slope*slope)));	
					}
				}
				if (d == 1){
					if (slope >= 0){
					totalbox *=Math.exp(-2*slope*x*deltas/(sq_sigma[d]*Math.sqrt(1+slope*slope)));
					secondtotalbox *=Math.exp(-2*2*slope*x*deltas/(sq_sigma[d]*Math.sqrt(1+slope*slope)));
					thirdtotalbox *=Math.exp(-2*3*slope*x*deltas/(sq_sigma[d]*Math.sqrt(1+slope*slope)));
					}
					else{
					totalbox *=Math.exp(-2*slope*x*deltas/(sq_sigma[d]*Math.sqrt(1+slope*slope)));
					secondtotalbox *=Math.exp(-2*2*slope*x*deltas/(sq_sigma[d]*Math.sqrt(1+slope*slope)));
					thirdtotalbox *=Math.exp(-2*3*slope*x*deltas/(sq_sigma[d]*Math.sqrt(1+slope*slope)));
					}
					
				}
				
				value *= Math.exp(-(x * x) / sq_sigma[d]) * (1 + totalbox + secondtotalbox + thirdtotalbox) ;
				
			
				
			
				
			}

			cursor.get().setReal(value);
		}
	
	
	}
	

}
package peakFitter;

import com.sun.tools.javac.util.Pair;

import houghandWatershed.PerformWatershedding;
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

	public static double[] sumofgaussianMaskFit(final RandomAccessibleInterval<FloatType> signalInterval,
			final RandomAccessibleInterval<IntType> intimg, final double[] location, final double[] sigma,
			final int iterations,  final double[] dxvector, final double slope,
			final double intercept,final double maxintensityline, final int numberofgaussians, final Endfit startorend, int label) throws Exception {
		final int n = signalInterval.numDimensions();

		// pre-compute sigma^2
		final double[] sq_sigma = new double[n];
		for (int d = 0; d < n; ++d)
			sq_sigma[d] = sigma[d] * sigma[d];

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
				beststartfitsumofGaussian(translatedIterableMask, location, sq_sigma, dxvector, slope,
						intercept, numberofgaussians);
				break;

			case End:
				bestendfitsumofGaussian(translatedIterableMask, location, sq_sigma, dxvector, slope,
						intercept, numberofgaussians);
				break;

			}

			final long[] longintlocation = new long[n];
			for (int d = 0; d < n; ++d) {
				longintlocation[d] = (long) location[d];

			}

			// ImageJFunctions.show(gaussianMask);
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
				final double signal = cImg.get().getRealDouble();
				final double mask = cMask.get().getRealDouble();
				final double weight = maxintensityline;

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

		} while (i < iterations);
		restoreBackground(signalIterable, bg);

		// ImageJFunctions.show(gaussianMask);

		return location;

	}

	public static double[] singlegaussianMaskFit(final RandomAccessibleInterval<FloatType> signalInterval,
			final RandomAccessibleInterval<IntType> intimg, final double[] location, final double[] sigma,
			final int iterations, final double maxintensity) {
		final int n = signalInterval.numDimensions();

		// pre-compute sigma^2
		final double[] sq_sigma = new double[n];
		for (int d = 0; d < n; ++d)
			sq_sigma[d] = sigma[d] * sigma[d];

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

			setsingleGaussian(translatedIterableMask, location, sq_sigma);

			// ImageJFunctions.show(gaussianMask);
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

				final double signal = cImg.get().getRealDouble();
				final double mask = cMask.get().getRealDouble();
				final double weight = maxintensity;

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

		} while (i < iterations);
		restoreBackground(signalIterable, bg);
		// ImageJFunctions.show(gaussianMask);
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

	final public static void beststartfitsumofGaussian(final IterableInterval<FloatType> image, final double[] location,
			final double[] sq_sigma,  final double[] dxvector, final double slope,
			final double intercept, int numberofgaussians) {
		final int ndims = image.numDimensions();
		final Cursor<FloatType> cursor = image.localizingCursor();

		double sumofgaussians = 0;
		while (cursor.hasNext()) {
			cursor.fwd();

			double value = 1.0;

			for (int d = 0; d < ndims; ++d) {
				final double x =  cursor.getDoublePosition(d) - location[d] ;
				sumofgaussians = Math.exp(-(x * x) / sq_sigma[d]);
				double y = 0;
				for (int i = 1; i < numberofgaussians; ++i) {

					y = x + i * dxvector[d];

					sumofgaussians += Math.exp(-(y * y) / sq_sigma[d]);

				}

				value *= sumofgaussians;

			}
			
			 
			cursor.get().setReal(value);

		}
	}

	final public static void bestendfitsumofGaussian(final IterableInterval<FloatType> image, final double[] location,
			final double[] sq_sigma,  final double[] dxvector, final double slope,
			final double intercept, int numberofgaussians) {
		final int ndims = image.numDimensions();
		final Cursor<FloatType> cursor = image.localizingCursor();
		double sumofgaussians = 0;
		while (cursor.hasNext()) {
			cursor.fwd();

			double value = 1.0;

			for (int d = 0; d < ndims; ++d) {
				final double x = cursor.getDoublePosition(d) - location[d] ;
				sumofgaussians = Math.exp(-(x * x) / sq_sigma[d]);
				double y = 0;
				for (int i = 1; i < numberofgaussians; ++i) {

					y = x - i * dxvector[d];

					sumofgaussians += Math.exp(-(y * y) / sq_sigma[d]);

				}

				value *= sumofgaussians;

			}
			
			  
			 
			cursor.get().setReal(value);

		}

	}

	
	final public static void setsingleGaussian(final IterableInterval<FloatType> image, final double[] location,
			final double[] sq_sigma) {
		final int numDimensions = image.numDimensions();

		final Cursor<FloatType> cursor = image.localizingCursor();

		while (cursor.hasNext()) {
			cursor.fwd();

			double value = 1.0;

			for (int d = 0; d < numDimensions; ++d) {
				final double x = location[d] - cursor.getIntPosition(d);

				value *= Math.exp(-(x * x) / sq_sigma[d]);

			}

			cursor.get().setReal(value);
		}
	}

}
package peakFitter;

import net.imglib2.Cursor;
import net.imglib2.IterableInterval;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;

public class GaussianMastFit {
	public static enum Endfit {
		Start, End
	}
	protected Endfit Endfit;
	public static double[] gaussianMaskFit(
			final RandomAccessibleInterval<FloatType> signalInterval,
			final RandomAccessibleInterval<IntType> intimg,
			final double[] location,
			final double[] sigma,
			final int iterations,
			final Endfit startorend)
	{
		final int n = signalInterval.numDimensions();
		
		// pre-compute 2*sigma^2
		final double[] two_sq_sigma = new double[ n ];
		for ( int d = 0; d < n; ++d )
			two_sq_sigma[ d ] = 2 * sigma[ d ] * sigma[ d ];

		// make the interval we fit on iterable
		final IterableInterval< FloatType > signalIterable = Views.iterable( signalInterval );
		
		// create the mask image
		final Img< FloatType > gaussianMask = new ArrayImgFactory< FloatType >().create( signalInterval, signalIterable.firstElement() );
		
		// set the mask image to the same location as the interval we fit on and make it iterable
		final long[] translation = new long[ n ];
		for ( int d = 0; d < n; ++d )
			translation[ d ] = signalInterval.min( d );
		
		final RandomAccessibleInterval<FloatType> translatedMask = Views.translate( gaussianMask, translation );
		final IterableInterval<FloatType> translatedIterableMask = Views.iterable( translatedMask );
		
		// remove background in the input
		final double bg = removeBackground( signalIterable );
		
		double N = 0;
		int i = 0;
		
		do
		{
			
			int [] intlocation = new int[n];
			
			for (int d = 0; d < n; ++d){
				
				intlocation[d] = (int) location[d];
				
			}
			
			RandomAccess<IntType> ranac = intimg.randomAccess();
			ranac.setPosition(intlocation);
			
			int label = ranac.get().get();
			
			switch (startorend){
			
			case Start:
				setstartGaussian( translatedIterableMask, location, two_sq_sigma );
			break;
			
			case End:
				setendGaussian( translatedIterableMask, location, two_sq_sigma );
			break;
				
			
			
			}
			
			
			
			// compute the sums
			final Cursor< FloatType > cMask = gaussianMask.cursor();
			final Cursor< FloatType > cImg = signalIterable.localizingCursor();

			double sumLocSN[] = new double[ n ]; // int_{all_px} d * S[ d ] * N[ d ]
			double sumSN = 0; // int_{all_px} S[ d ] * N[ d ]
			double sumSS = 0; // int_{all_px} S[ d ] * S[ d ]
			
			while ( cMask.hasNext() )
			{
				cMask.fwd();
				cImg.fwd();
				
				ranac.setPosition(cImg);
				if (ranac.get().get() == label){
				final double signal = cImg.get().getRealDouble();
				final double mask = cMask.get().getRealDouble();
				final double weight = 10;
				
				final double signalmask = signal * mask * weight;
				
				sumSN += signalmask;
				sumSS += signal * signal * weight;
				
				for ( int d = 0; d < n; ++d )
				{
					final double l = cImg.getDoublePosition( d );
					sumLocSN[ d ] += l * signalmask;
				}
				
			}
			
			for ( int d = 0; d < n; ++d )
				location[ d ] = sumLocSN[ d ] / sumSN;
			
			N = sumSN / sumSS;
			
			
			
			++i;
		
		}
		}
		while ( i < iterations );
		restoreBackground( signalIterable, bg );
		
		
		return location;
	}
	
	public static double removeBackground( final IterableInterval< FloatType > iterable )
	{
		double i = 0;
		
		for ( final FloatType t : iterable )
			i += t.getRealDouble();
		
		i /= (double)iterable.size();
		
		for ( final FloatType t : iterable )
			t.setReal( t.get() - i );
		
		return i;
	}

	public static void restoreBackground( final IterableInterval< FloatType > iterable, final double value )
	{
		for ( final FloatType t : iterable )
			t.setReal( t.get() + value );
	}

	final public static void setstartGaussian( final IterableInterval< FloatType > image, final double[] location, final double[] two_sq_sigma )
	{
		final int numDimensions = image.numDimensions();
		
		final Cursor< FloatType > cursor = image.localizingCursor();
		
		while ( cursor.hasNext() )
		{
			cursor.fwd();
			
			double value = 1;
			
			for ( int d = 0; d < numDimensions; ++d )
			{
				final double x = location[ d ] - cursor.getDoublePosition( d );
				
				if (x <=  0)
					
				value *= Math.exp( -(x * x) / two_sq_sigma[ d ] );
				
				else value = 0;
				
			}
			
			cursor.get().setReal( value );
		}
	}
	final public static void setendGaussian( final IterableInterval< FloatType > image, final double[] location, final double[] two_sq_sigma )
	{
		final int numDimensions = image.numDimensions();
		
		final Cursor< FloatType > cursor = image.localizingCursor();
		
		while ( cursor.hasNext() )
		{
			cursor.fwd();
			
			double value = 1;
			
			for ( int d = 0; d < numDimensions; ++d )
			{
				final double x = location[ d ] - cursor.getDoublePosition( d );
				
				if ( x >= 0)
					
				value *= Math.exp( -(x * x) / two_sq_sigma[ d ] );
				
				else value = 0;
			}
			
			cursor.get().setReal( value );
		}
	}
	
	
}

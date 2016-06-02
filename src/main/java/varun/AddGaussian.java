package varun;

import net.imglib2.Cursor;
import net.imglib2.IterableInterval;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;

public class AddGaussian {
	
	
	final public static void addGaussian( final RandomAccessibleInterval< FloatType > image, final double[] location, final double[] sigma,
			boolean Normalize )
	{
	final int numDimensions = image.numDimensions();
	final int[] size = new int[ numDimensions ];

	final long[] min = new long[ numDimensions ];
	final long[] max = new long[ numDimensions ];

	final double[] two_sq_sigma = new double[ numDimensions ];

	for ( int d = 0; d < numDimensions; ++d )
	{
	size[ d ] = getSuggestedKernelDiameter( sigma[ d ] ) * 2;
	min[ d ] = (int)Math.round( location[ d ] ) - size[ d ]/2;
	max[ d ] = min[ d ] + size[ d ] - 1;
	two_sq_sigma[ d ] = 2 * sigma[ d ] * sigma[ d ];
	}

	final RandomAccessible< FloatType > infinite = Views.extendZero( image );
	final RandomAccessibleInterval< FloatType > interval = Views.interval( infinite, min, max );
	final IterableInterval< FloatType > iterable = Views.iterable( interval );
	final Cursor< FloatType > cursor = iterable.localizingCursor();
	double sum = 0;
	while ( cursor.hasNext() )
	{
	cursor.fwd();

	double value = 1;

	for ( int d = 0; d < numDimensions; ++d )
	{
	final double x = location[ d ] - cursor.getIntPosition( d );
	value *= Math.exp( -(x * x) / two_sq_sigma[ d ] );
	}
	
	if (Normalize){
      
           sum += value;
       
           value /= sum;
	}
	
	cursor.get().set( cursor.get().get() + (float)value );
	
	
	}
	
	
	
	
	}

	public static int getSuggestedKernelDiameter( final double sigma )
	{
	int size = 3;
    int cutoff = 3; // This number means cutoff is chosen to be cutoff times sigma. 
    if ( sigma > 0 )
	size = Math.max( cutoff, ( 2 * ( int ) ( cutoff * sigma + 0.5 ) + 1 ) );

	return size;
	}
	
	
	final public static void add2DaxisGaussian( final RandomAccessibleInterval< FloatType > image, final double[] location, final double[] sigma, double slope )
	{
	final int numDimensions = image.numDimensions();

	final long[] min = new long[ numDimensions ];
	final long[] max = new long[ numDimensions ];

	final double[] two_sq_sigma = new double[ numDimensions ];

	for ( int d = 0; d < numDimensions; ++d )
	{
	
	two_sq_sigma[ d ] = 2 * sigma[ d ] * sigma[ d ];
	}

	final RandomAccessible< FloatType > infinite = Views.extendZero( image );
	final RandomAccessibleInterval< FloatType > interval = Views.interval( infinite, min, max );
	final IterableInterval< FloatType > iterable = Views.iterable( interval );
	final Cursor< FloatType > cursor = iterable.localizingCursor();
	while ( cursor.hasNext() )
	{
	cursor.fwd();

	double value = 1;

	final double sintheta = slope / Math.sqrt(1+slope*slope);
	final double costheta = 1.0 / Math.sqrt(1+slope*slope);
	final double x = location[ 0 ] - cursor.getIntPosition( 0 );
	final double y = location[1] - cursor.getIntPosition( 1 );
	final double xprime = x*costheta  - y*sintheta;
	final double yprime = x*sintheta +y*costheta;
	value *= Math.exp( -(xprime * xprime) / two_sq_sigma[ 0 ] );
	
	value *= Math.exp( -(yprime * yprime) / two_sq_sigma[ 1 ] );
	
	cursor.get().set( cursor.get().get() + (float)value );
	
	
	}
	
	
	
	
	}
	
	
	
}

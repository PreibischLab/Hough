package varun;

import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.Interval;
import net.imglib2.IterableInterval;
import net.imglib2.Localizable;
import net.imglib2.Point;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.neighborhood.Neighborhood;
import net.imglib2.algorithm.neighborhood.RectangleShape;
import net.imglib2.algorithm.region.hypersphere.HyperSphere;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.Cursor;
import net.imglib2.Point;
import net.imglib2.PointSampleList;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Intervals;
import net.imglib2.view.Views;

public class GetLocalmaxima {

	public static void Thresholding(RandomAccessibleInterval<FloatType> img, RandomAccessibleInterval<FloatType> imgout,
			FloatType ThresholdValue) {

		final Cursor<FloatType> bound = Views.iterable(img).localizingCursor();

		final RandomAccess<FloatType> outbound = imgout.randomAccess();

		while (bound.hasNext()) {

			bound.fwd();

			outbound.setPosition(bound);

			if (bound.get().compareTo(ThresholdValue) > 0) {

				outbound.get().set(bound.get());

			}

			else {

				outbound.get().setZero();

			}

		}
	}

	 public static  Img<FloatType> FindandDisplayLocalMaxima(RandomAccessibleInterval<FloatType> img,
         ImgFactory<FloatType> imageFactory )
 {
     // Create a new image for the output
     Img<FloatType> output = imageFactory.create( img, new FloatType() );

     // define an interval that is one pixel smaller on each side in each dimension
     
     Interval interval = Intervals.expand( img, -1 );

     // create a view on the source with this interval
     img = Views.interval( img, interval );

     // create a Cursor that iterates over the source and checks in a 8-neighborhood
     // if it is a maxima
     final Cursor< FloatType > center = Views.iterable(img).cursor();

     // instantiate a RectangleShape to access rectangular local neighborhoods
    
     final RectangleShape shape = new RectangleShape( 1, true );

     // iterate over the set of neighborhoods in the image
     for ( final Neighborhood<FloatType> localNeighborhood : shape.neighborhoods(img) )
     {
         // what is the value that we investigate?
         // (the center cursor runs over the image in the same iteration order as neighborhood)
         final FloatType centerValue = center.next();

         // keep this boolean true as long as no other value in the local neighborhood
         // is larger or equal
         boolean isMaximum = true;

         // check if all pixels in the local neighborhood that are smaller
         for ( final FloatType value : localNeighborhood )
         {
             // test if the center is smaller than the current pixel value
             if ( centerValue.compareTo( value ) < 0 )
             {
                 isMaximum = false;
                 break;
             }
         }

         if ( isMaximum )
         {
        	 final RandomAccess<FloatType> outbound = output.randomAccess();
        	 outbound.setPosition(center);
        	 outbound.get().set(center.get());
          
         }
     }

     return output;
 }
}

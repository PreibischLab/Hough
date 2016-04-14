package varun;

import java.util.ArrayList;

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
import net.imglib2.RealPoint;
import net.imglib2.RealPointSampleList;
import net.imglib2.algorithm.localextrema.RefinedPeak;
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

public class GetLocalmaxmin {

	public static ArrayList<RefinedPeak<Point>> Removesimilar(ArrayList<RefinedPeak<Point>> SubpixelMinlist, 
			double thetatolerance, double rhotolerance){
		/********
		 * The part below removes the close values in theta and rho coordinate 
		 * (keeps only single of the multiple values)
		 ********/

		int j = 0;

		for (int i = 0; i < SubpixelMinlist.size(); ++i) {

			j = i + 1;
			while (j < SubpixelMinlist.size()) {

				if (Math.abs(SubpixelMinlist.get(i).getDoublePosition(0) - SubpixelMinlist.get(j)
						.getDoublePosition(0))<thetatolerance && Math.abs(SubpixelMinlist.get(i).getDoublePosition(1) - SubpixelMinlist.get(j)
								.getDoublePosition(1))<rhotolerance ) {

					SubpixelMinlist.remove(j);

				}

				else {
					++j;
				}

			}

		}
		
		
		return SubpixelMinlist;
	}
	
	
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

	public static Img<FloatType> FindandDisplayLocalMaxima(RandomAccessibleInterval<FloatType> img,
			ImgFactory<FloatType> imageFactory) {
		// Create a new image for the output
		Img<FloatType> output = imageFactory.create(img, new FloatType());

		// define an interval that is span number of pixel smaller on each side
		// in each dimension
		int span = 1;

		Interval interval = Intervals.expand(img, -span);

		// create a view on the source with this interval
		img = Views.interval(img, interval);

		// create a Cursor that iterates over the source and checks in a
		// 8-neighborhood
		// if it is a maxima
		final Cursor<FloatType> center = Views.iterable(img).cursor();

		// instantiate a RectangleShape to access rectangular local
		// neighborhoods

		final RectangleShape shape = new RectangleShape(span, true);

		// iterate over the set of neighborhoods in the image
		for (final Neighborhood<FloatType> localNeighborhood : shape.neighborhoods(img)) {
			final FloatType centerValue = center.next();

			// keep this boolean true as long as no other value in the local
			// neighborhood
			// is smaller
			boolean isMaximum = true;

			// check if all pixels in the local neighborhood that are smaller
			for (final FloatType value : localNeighborhood) {
				// test if the center is smaller than the current pixel value
				if (centerValue.compareTo(value) < 0) {
					isMaximum = false;
					break;
				}
			}

			if (isMaximum) {
				final RandomAccess<FloatType> outbound = output.randomAccess();
				outbound.setPosition(center);
				outbound.get().set(center.get());

			}
		}

		return output;
	}
	 public static  Img<FloatType> FindandDisplayLocalMinima(RandomAccessibleInterval<FloatType> img,
	         ImgFactory<FloatType> imageFactory )
	 {
	     // Create a new image for the output
	     Img<FloatType> output = imageFactory.create( img, new FloatType() );

	     // define an interval that is span number of pixel smaller on each side in each dimension
	     int span = 1;
	     
	     Interval interval = Intervals.expand( img, -span );

	     // create a view on the source with this interval
	     img = Views.interval( img, interval );

	     // create a Cursor that iterates over the source and checks in a 8-neighborhood
	     // if it is a maxima
	     final Cursor< FloatType > center = Views.iterable(img).cursor();

	     // instantiate a RectangleShape to access rectangular local neighborhoods
	    
	     final RectangleShape shape = new RectangleShape( span, true );

	     // iterate over the set of neighborhoods in the image
	     for ( final Neighborhood<FloatType> localNeighborhood : shape.neighborhoods(img) )
	     {
	         final FloatType centerValue = center.next();

	         // keep this boolean true as long as no other value in the local neighborhood
	         // is smaller
	         boolean isMinimum = true;

	         // check if all pixels in the local neighborhood that are smaller
	         for ( final FloatType value : localNeighborhood )
	         {
	             // test if the center is smaller than the current pixel value
	             if ( centerValue.compareTo( value ) >= 0 )
	             {
	                 isMinimum = false;
	                 break;
	             }
	         }

	         if ( isMinimum )
	         {
	        	 final RandomAccess<FloatType> outbound = output.randomAccess();
	        	 outbound.setPosition(center);
	        	 outbound.get().set(center.get());
	        	 
	          
	         }
	     }

	     return output;
	 }	
	
	 
	 public static  ArrayList<RealPoint> FindLocalMinima(RandomAccessibleInterval<FloatType> img)
	 {
		 
		 int n = img.numDimensions();
		 
	     ArrayList<RealPoint> Minlist = new ArrayList<RealPoint>(n);
		 
	     // define an interval that is span number of pixel smaller on each side in each dimension
	     int span = 1;
	     
	     Interval interval = Intervals.expand( img, -span );

	     // create a view on the source with this interval
	     img = Views.interval( img, interval );

	     // create a Cursor that iterates over the source and checks in a 8-neighborhood
	     // if it is a maxima
	     final Cursor< FloatType > center = Views.iterable(img).cursor();

	     // instantiate a RectangleShape to access rectangular local neighborhoods
	    
	     final RectangleShape shape = new RectangleShape( span, true );

	     // iterate over the set of neighborhoods in the image
	     for ( final Neighborhood<FloatType> localNeighborhood : shape.neighborhoods(img) )
	     {
	         final FloatType centerValue = center.next();

	         // keep this boolean true as long as no other value in the local neighborhood
	         // is smaller
	         boolean isMinimum = true;

	         // check if all pixels in the local neighborhood that are smaller
	         for ( final FloatType value : localNeighborhood )
	         {
	             // test if the center is smaller than the current pixel value
	             if ( centerValue.compareTo( value ) >= 0 )
	             {
	                 isMinimum = false;
	                 break;
	             }
	         }

	         if ( isMinimum )
	         {
	        		 RealPoint Minpoints = new RealPoint(center);
	        	     Minlist.add(Minpoints);
	         }
	     }

	     return Minlist;
	 }	
	
	
}

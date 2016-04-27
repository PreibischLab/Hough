package varun;

import net.imglib2.Cursor;
import net.imglib2.KDTree;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealPoint;
import net.imglib2.RealPointSampleList;
import net.imglib2.algorithm.labeling.Watershed;
import net.imglib2.algorithm.stats.Normalize;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.labeling.Labeling;
import net.imglib2.neighborsearch.NearestNeighborSearchOnKDTree;
import net.imglib2.roi.labeling.ImgLabeling;
import net.imglib2.roi.labeling.LabelingType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.integer.ShortType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import varun.GetLocalmaxmin.IntensityType;

public class PerformWatershedding {
	
	public static enum InverseType {
		Straight, Inverse
	}
	public static void DistanceTransformImage(RandomAccessibleInterval<FloatType> inputimg,
			RandomAccessibleInterval<FloatType> outimg, final InverseType invtype){
int n = inputimg.numDimensions();
		
		final Img<BitType> bitimg = new ArrayImgFactory<BitType>().create(inputimg, new BitType());
		// make an empty list
				final RealPointSampleList< BitType > list = new RealPointSampleList< BitType >( n );

				FloatType threshold = new FloatType(100);
				GetLocalmaxmin.ThresholdingBit(inputimg, bitimg,
						threshold);
				
				// cursor on the binary image
				final Cursor< BitType > cursor = bitimg.localizingCursor();

				// for every pixel that is 1, make a new RealPoint at that location
				while ( cursor.hasNext() )
					if ( cursor.next().getInteger() == 1 )
						list.add( new RealPoint( cursor ), cursor.get() );

				// build the KD-Tree from the list of points that == 1
				final KDTree< BitType > tree = new KDTree< BitType >( list );

				// Instantiate a nearest neighbor search on the tree (does not modifiy the tree, just uses it)
				final NearestNeighborSearchOnKDTree< BitType > search = new NearestNeighborSearchOnKDTree< BitType >( tree );

				// randomaccess on the output
				final RandomAccess< FloatType > ranac = outimg.randomAccess();

				// reset cursor for the input (or make a new one)
				cursor.reset();

				// for every pixel of the binary image
				while ( cursor.hasNext() )
				{
					cursor.fwd();

					// set the randomaccess to the same location
					ranac.setPosition( cursor );

					// if value == 0, look for the nearest 1-valued pixel
					if ( cursor.get().getInteger() == 0 )
					{
						// search the nearest 1 to the location of the cursor (the current 0)
						search.search( cursor );

						// get the distance (the previous call could return that, this for generality that it is two calls)
						switch(invtype){
						
						case Straight:
							ranac.get().setReal( search.getDistance() );
							break;
						case Inverse:
							ranac.get().setReal( -search.getDistance() );
							break;
							
						}
						
						
					}
					else
					{
						// if value == 1, no need to search 
						ranac.get().setZero();
					}
				}
		
	}
	
	
	
	public static void WatershedImage(RandomAccessibleInterval<FloatType> inputimg,RandomAccessibleInterval<FloatType> seedimg){
	
		int n = inputimg.numDimensions();
		long[] dimensions = new long[] {inputimg.dimension(0),inputimg.dimension(1)};
		
		RandomAccessibleInterval<FloatType> outputimg = new ArrayImgFactory<FloatType>().create(inputimg, new FloatType());
		ImgLabeling<Double, ShortType> seeds =
	            new ImgLabeling<Double, ShortType>(new ArrayImgFactory<ShortType>().create(seedimg, new ShortType()));
		ImgLabeling< Double, ShortType > seedLabeling = 
				new ImgLabeling< Double, ShortType >( new ArrayImgFactory< ShortType >().create( dimensions, new ShortType() ) );
		RandomAccessibleInterval<LabelingType<Float>> watershedResult =
	            new ImgLabeling<Float, ShortType>(new ArrayImgFactory<ShortType>().create(inputimg, new ShortType()));
		
		final Cursor<LabelingType<Double>> seedcursor = seeds.localizingCursor();
		
		
		
	
	}
}

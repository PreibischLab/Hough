package preProcessing;

import ij.ImagePlus;
import ij.ImageStack;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.IterableInterval;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.BenchmarkAlgorithm;
import net.imglib2.algorithm.OutputAlgorithm;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Util;
import net.imglib2.view.IntervalView;
import net.imglib2.view.Views;

public class Addstacks extends BenchmarkAlgorithm implements OutputAlgorithm< RandomAccessibleInterval< FloatType >> {

	
	private static final String BASE_ERROR_MSG = "[CouldnotBiggify] ";
	private final RandomAccessibleInterval<FloatType> sourceA;
	private final RandomAccessibleInterval<FloatType> sourceB;
	private RandomAccessibleInterval< FloatType > output;
	public Addstacks(RandomAccessibleInterval<FloatType> sourceA, RandomAccessibleInterval<FloatType> sourceB){
		
		this.sourceA = sourceA;
		this.sourceB = sourceB;
	}
	
	
	
	@Override
	public boolean process() {
		 final int ndims = sourceA.numDimensions();
		final FloatType type = sourceA.randomAccess().get().createVariable();
			final ImgFactory< FloatType > factory = Util.getArrayOrCellImgFactory( sourceA, type );
			 this.output = factory.create( sourceA, type );
		
			
			
			
		   
		 
			

				
				
			  final long[] dimensions = {sourceA.dimension(0) , sourceA.dimension(1) , sourceA.dimension(2)};
			
			  
			   
			   
				for (int frame = 0; frame < sourceA.dimension(ndims - 1); ++frame) {

					IntervalView<FloatType> currentframeA = Views.hyperSlice(sourceA, ndims - 1, frame);
					IntervalView<FloatType> currentframeB = Views.hyperSlice(sourceB, ndims - 1, frame);
					IntervalView<FloatType> currentframeout = Views.hyperSlice(output, ndims - 1, frame);
					
					processSlice( currentframeA, currentframeB, currentframeout);
					
					
				}
				
				  
				   
		   
		
		return true;
	}
	
	
	private static void processSlice(RandomAccessibleInterval<FloatType> currentframeA,RandomAccessibleInterval<FloatType> currentframeB, RandomAccessibleInterval<FloatType> currentframeout) {
		int ndims = currentframeA.numDimensions();
		
		
		
			final Cursor< FloatType > cursorImg1 = Views.iterable(currentframeA).localizingCursor();
			final RandomAccess< FloatType > ra1 = currentframeB.randomAccess();
			final RandomAccess< FloatType > ra2 = currentframeout.randomAccess();
        while ( cursorImg1.hasNext())
        {
            cursorImg1.fwd();
            ra2.setPosition( cursorImg1 );
            ra1.setPosition(cursorImg1);
            ra2.get().set( cursorImg1.get().get() + ra1.get().get());
        }
        
		
		
	}
	

	


	@Override
	public boolean checkInput() {
		if ( sourceA.numDimensions() > 3 )
		{
			errorMessage = BASE_ERROR_MSG + " Can only operate on 1D, 2D or 3D images. Got " + sourceA.numDimensions() + "D.";
			return false;
		}
		return true;
	}


	


	@Override
	public RandomAccessibleInterval<FloatType> getResult() {
		return output;
	}
	
	
	
	
}

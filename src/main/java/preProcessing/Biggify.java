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

public class Biggify extends BenchmarkAlgorithm implements OutputAlgorithm< RandomAccessibleInterval< FloatType >> {

	
	private static final String BASE_ERROR_MSG = "[CouldnotBiggify] ";
	private final RandomAccessibleInterval<FloatType> source;
	private RandomAccessibleInterval< FloatType > output;
	private RandomAccessibleInterval< FloatType > bigoutput;
	private final int numberofPixels;
	public Biggify(RandomAccessibleInterval<FloatType> source, final int numberofPixels){
		
		this.source = source;
		this.numberofPixels = numberofPixels;
	}
	
	
	
	@Override
	public boolean process() {
		 final int ndims = source.numDimensions();
		final FloatType type = source.randomAccess().get().createVariable();
			final ImgFactory< FloatType > factory = Util.getArrayOrCellImgFactory( source, type );
			 this.output = factory.create( source, type );
		
			
			
			
		   
		   if (ndims > 2 ){
			

				
				
			  final long[] dimensions = {source.dimension(0) + numberofPixels, source.dimension(1) +numberofPixels, source.dimension(2)};
			 this.bigoutput = factory.create(dimensions, type);
			
			  
			   
			   
				for (int frame = 0; frame < source.dimension(ndims - 1); ++frame) {

					IntervalView<FloatType> currentframe = Views.hyperSlice(source, ndims - 1, frame);
					IntervalView<FloatType> currentframeout = Views.hyperSlice(bigoutput, ndims - 1, frame);
					
					processSlice( currentframe, currentframeout, numberofPixels );
					
					
				}
				
				  
				   
		   }
		   
		   else{
			   final long[] dimensions = {source.dimension(0) + numberofPixels, source.dimension(1) +numberofPixels};
			   this.bigoutput = factory.create( dimensions, type );
			   processSlice(source, output, numberofPixels);
			   
		   }
		
		return true;
	}
	
	
	private static void processSlice(RandomAccessibleInterval<FloatType> currentframe, RandomAccessibleInterval<FloatType> currentframeout, int numberofPixels) {
		int ndims = currentframe.numDimensions();
		
		
		
			final Cursor< FloatType > cursorImg1 = Views.iterable(currentframe).localizingCursor();
			final RandomAccess< FloatType > ra2 = currentframeout.randomAccess();
        while ( cursorImg1.hasNext())
        {
            cursorImg1.fwd();
            ra2.setPosition( cursorImg1 );
            ra2.get().set( cursorImg1.get());
        }
        
		
		
	}
	

	public static IntervalView<FloatType> biggifyimage(IntervalView<FloatType> inputimg, final int numberofPixels){
		  
		
		   final int ndims = inputimg.numDimensions();
		   
		   if (ndims == 2){
		   long[] min = new long[ inputimg.numDimensions() ];
	       long[] max = new long[ inputimg.numDimensions() ];
		   for (int d = 0; d < ndims; ++d){
			min[d] = inputimg.min(d) - numberofPixels;
			max[d] = inputimg.max(d) + numberofPixels;
  		}
		
		FinalInterval interval = new FinalInterval(min, max);
		 
		IntervalView<FloatType> outimg  = Views.offsetInterval(Views.extendBorder(inputimg), interval);
		
		return outimg;
		   }
		  if (ndims > 2 ){
			  
			   IntervalView<FloatType> bigimg = inputimg;
			   
				for (int frame = 1; frame < inputimg.dimension(ndims - 1); ++frame) {

					IntervalView<FloatType> currentframe = Views.hyperSlice(inputimg, ndims - 1, frame);
					IntervalView<FloatType> currentframeout = Views.hyperSlice(bigimg, ndims - 1, frame);
					processSlice( currentframe, currentframeout, numberofPixels );
				
			  
			  
		  }
				
				return bigimg;
		  }
		  
		  else
			  return null;
		   
	}


	@Override
	public boolean checkInput() {
		if ( source.numDimensions() > 3 )
		{
			errorMessage = BASE_ERROR_MSG + " Can only operate on 1D, 2D or 3D images. Got " + source.numDimensions() + "D.";
			return false;
		}
		return true;
	}


	


	@Override
	public RandomAccessibleInterval<FloatType> getResult() {
		return bigoutput;
	}
	
	
	
	
}

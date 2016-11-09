package preProcessing;

import net.imglib2.FinalInterval;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.IntervalView;
import net.imglib2.view.Views;

public class Biggify {

	
	public static RandomAccessibleInterval<FloatType> biggifyimage(RandomAccessibleInterval<FloatType> inputimg, final int numberofPixels){
		  
		
		   final int ndims = inputimg.numDimensions();
		   long[] min = new long[ inputimg.numDimensions() ];
	       long[] max = new long[ inputimg.numDimensions() ];
		   for (int d = 0; d < ndims; ++d){
			min[d] = inputimg.min(d) - numberofPixels;
			max[d] = inputimg.max(d) + numberofPixels;
     		}
		
		FinalInterval interval = new FinalInterval(min, max);
		 
		RandomAccessibleInterval<FloatType> outimg  = Views.offsetInterval(Views.extendBorder(inputimg), interval);
		
		return outimg;
	}
	
	
	public static IntervalView<FloatType> biggifyimage(IntervalView<FloatType> inputimg, final int numberofPixels){
		  
		
		   final int ndims = inputimg.numDimensions();
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
	
	
	
	
}

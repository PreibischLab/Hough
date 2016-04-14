package varun;

import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;

public class SplitImages {

	public static RandomAccessibleInterval<FloatType> Splitk( RandomAccessibleInterval<FloatType> inputimg,
			int[] index, int[] size){
		
		final long[] dimensions = new long[] { inputimg.dimension(0)/size[0], inputimg.dimension(1)/size[1] };
		
		RandomAccessibleInterval< FloatType > outimg =
                Views.interval( inputimg, new long[] { index[0], index[1] }, new long[]{ index[0]+dimensions[0], index[1]+dimensions[1] } );

		return outimg;
	}
	
}

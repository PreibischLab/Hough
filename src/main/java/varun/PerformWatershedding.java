package varun;

import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.roi.labeling.ImgLabeling;
import net.imglib2.roi.labeling.LabelingType;
import net.imglib2.type.numeric.integer.ShortType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;

public class PerformWatershedding {
	
	
	public static void WatershedImage(RandomAccessibleInterval<FloatType> inputimg){
	
		RandomAccessibleInterval<FloatType> outputimg = new ArrayImgFactory<FloatType>().create(inputimg, new FloatType());
		ImgLabeling<Double, ShortType> seeds =
	            new ImgLabeling<Double, ShortType>(new ArrayImgFactory<ShortType>().create(inputimg, new ShortType()));
		
		RandomAccessibleInterval<LabelingType<Float>> watershedResult =
	            new ImgLabeling<Float, ShortType>(new ArrayImgFactory<ShortType>().create(inputimg, new ShortType()));
	
		
		
	
	}
}

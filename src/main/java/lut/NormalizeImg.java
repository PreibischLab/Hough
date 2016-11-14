package lut;

import java.io.File;

import ij.ImageJ;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.algorithm.stats.Normalize;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import preProcessing.Addstacks;

public class NormalizeImg {

	
	public static void main(String[] args) throws Exception {
		
		new ImageJ();

		// Load the stack of images
		RandomAccessibleInterval<FloatType> img = util.ImgLib2Util.openAs32Bit(
				
				new File("../res/C1-23-09-16-laevis-6uM-25C-dup.tif"), 
				new ArrayImgFactory<FloatType>());
		
        RandomAccessibleInterval<FloatType> secimg = util.ImgLib2Util.openAs32Bit(
				
				new File("../res/C2-23-09-16-laevis-6uM-25C-dup.tif"), 
				new ArrayImgFactory<FloatType>());
		
		
        
        FloatType minval = new FloatType(0);
		FloatType maxval = new FloatType(1);
		Normalize.normalize(Views.iterable(img), minval, maxval);
		Normalize.normalize(Views.iterable(secimg), minval, maxval);
		Addstacks add = new Addstacks(img, secimg);
		add.process();
		RandomAccessibleInterval<FloatType> imgout = add.getResult();
		Normalize.normalize(Views.iterable(imgout), minval, maxval);
		ImageJFunctions.show(imgout);
        
	}
	
}

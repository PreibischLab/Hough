package lut;

import java.io.File;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
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
				
				new File("../res/C1-Bovine_12uM_37C-1.tif"), 
				new ArrayImgFactory<FloatType>());
		
        RandomAccessibleInterval<FloatType> secimg = util.ImgLib2Util.openAs32Bit(
				
				new File("../res/C2-Bovine_12uM_37C-1.tif"), 
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

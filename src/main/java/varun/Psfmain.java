package varun;

import java.io.File;
import java.util.ArrayList;

import ij.ImageJ;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.gauss.Gauss;
import net.imglib2.algorithm.stats.Normalize;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import util.ImgLib2Util;
import varun.Kernels.ProcessingType;

public class Psfmain {
	
	
	public static void main(String[] args) throws Exception {
		RandomAccessibleInterval<FloatType> biginputimg = ImgLib2Util
				.openAs32Bit(new File("src/main/resources/Fresh_data/psf_488_02.tif"));
		// small_mt.tif image to be used for testing
		// 2015-01-14_Seeds-1.tiff for actual
		// mt_experiment.tif for big testing
		new ImageJ();
      
		new Normalize();
		FloatType minval = new FloatType(0);
		FloatType maxval = new FloatType(1);
		Normalize.normalize(Views.iterable(biginputimg), minval, maxval);
		ImageJFunctions.show(biginputimg);
		// Initialize empty images to be used later
				RandomAccessibleInterval<FloatType> inputimg = new ArrayImgFactory<FloatType>().create(biginputimg,
						new FloatType());
				
				 
				 inputimg = Kernels.Preprocess(biginputimg, ProcessingType.SupressThresh);
				 ImageJFunctions.show(inputimg);
		Extractpsfinfo getpsf = new Extractpsfinfo(inputimg);
		
		
		ArrayList<double[]> totalgausslist = new ArrayList<double[]>();
		getpsf.Extractparams(totalgausslist);
		
	}
}

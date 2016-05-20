package varun;

import java.io.File;

import ij.ImageJ;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.gauss3.Gauss3;
import net.imglib2.algorithm.stats.Normalize;
import net.imglib2.exception.IncompatibleTypeException;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import util.ImgLib2Util;
import varun.Kernels.ProcessingType;

public class HoughpostWater {
	
	public static void main(String[] args) throws IncompatibleTypeException {

		RandomAccessibleInterval<FloatType> biginputimg = ImgLib2Util
				.openAs32Bit(new File("src/main/resources/2015-01-14_Seeds-1.tiff"));

		new ImageJ();
		
		new Normalize();
		FloatType minval = new FloatType(0);
		FloatType maxval = new FloatType(1);
		Normalize.normalize(Views.iterable(biginputimg), minval, maxval);
		ImageJFunctions.show(biginputimg);
		double[] sigma = new double[ biginputimg.numDimensions()];

		for ( int d = 0; d < sigma.length; ++d )
			sigma[ d ] = 1;
		
		RandomAccessibleInterval<FloatType> inputimg = new ArrayImgFactory<FloatType>().create(biginputimg, new FloatType());
		RandomAccessibleInterval<FloatType> tmpinputimg = new ArrayImgFactory<FloatType>().create(biginputimg, new FloatType());
		
		// Preprocess image
		tmpinputimg = Kernels.Preprocess(biginputimg, ProcessingType.Meanfilter);
		inputimg = Kernels.Preprocess(tmpinputimg, ProcessingType.CannyEdge);
		
		
		ImageJFunctions.show(inputimg).setTitle("Preprocessed image");
		
		// Do watershedding and Hough
		
		

		PerformWatershedding.DowatersheddingandHough(biginputimg,inputimg);
		
		
	
	}
}
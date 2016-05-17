package varun;

import java.io.File;

import ij.ImageJ;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.stats.Normalize;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import util.ImgLib2Util;
import varun.Kernels.ProcessingType;

public class HoughpostWater {
	
	public static void main(String[] args) {

		RandomAccessibleInterval<FloatType> biginputimg = ImgLib2Util
				.openAs32Bit(new File("src/main/resources/2015-01-14_Seeds-1.tiff"));

		new ImageJ();
		
		new Normalize();
		FloatType minval = new FloatType(0);
		FloatType maxval = new FloatType(1);
		Normalize.normalize(Views.iterable(biginputimg), minval, maxval);
		ImageJFunctions.show(biginputimg);
		
		RandomAccessibleInterval<FloatType> tmpinputimg = new ArrayImgFactory<FloatType>().create(biginputimg, new FloatType());
		RandomAccessibleInterval<FloatType> inputimg = new ArrayImgFactory<FloatType>().create(biginputimg, new FloatType());
		
		// Preprocess image
	//	tmpinputimg = Kernels.Preprocess(biginputimg, ProcessingType.Gradientmag);
	//	inputimg = Kernels.Preprocess(tmpinputimg, ProcessingType.SupressThresh);
		
		inputimg = Kernels.Preprocess(biginputimg, ProcessingType.SupressThresh);
		ImageJFunctions.show(inputimg).setTitle("Preprocessed image");
		
		// Do watershedding and Hough
		
		// Set size of pixels in Hough space
		int mintheta = 0;
		// Usually is 180 but to allow for detection of vertical lines,allowing
		// // a few more degrees
		int maxtheta = 200;
		double size = Math
				.sqrt((inputimg.dimension(0) * inputimg.dimension(0) + inputimg.dimension(1) * inputimg.dimension(1)));
		int minRho = (int) -Math.round(size);
		int maxRho = -minRho;
				double thetaPerPixel = 1;
				double rhoPerPixel = 1;
				double[] min = { mintheta, minRho };
				double[] max = { maxtheta, maxRho };
				int pixelsTheta = (int) Math.round((maxtheta - mintheta) / thetaPerPixel);
				int pixelsRho = (int) Math.round((maxRho - minRho) / rhoPerPixel);

		PerformWatershedding.DowatersheddingandHough(inputimg, min, max, pixelsTheta, pixelsRho, thetaPerPixel, rhoPerPixel);
		
		
	
	}
}
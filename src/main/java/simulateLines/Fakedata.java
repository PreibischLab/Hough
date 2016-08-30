package simulateLines;

import java.io.File;
import java.util.Random;

import ij.ImageJ;
import net.imglib2.FinalInterval;
import net.imglib2.Interval;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.stats.Normalize;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Intervals;
import net.imglib2.view.Views;
import poissonSimulator.Poissonprocess;
import preProcessing.Kernels;
import util.ImgLib2Util;

public class Fakedata {

	
	public static void main (String[] args){
		
		new ImageJ();
		
		final FinalInterval range = new FinalInterval(200, 200);
		
		
		RandomAccessibleInterval<FloatType> imgout = new ArrayImgFactory<FloatType>().create(range, new FloatType());
		RandomAccessibleInterval<FloatType> noisyimg = new ArrayImgFactory<FloatType>().create(imgout, new FloatType());
		final int ndims = imgout.numDimensions();
		final Random rnd = new Random(150);
		final double [] sigma = {1.7,1.8};
		final double [] Ci = new double[ndims];
		
		for (int d = 0; d < ndims; ++d)
			Ci[d] = 1.0 / Math.pow(sigma[d],2);
		
		Kernels.SaltandPepperNoise(imgout);
		final int numlines = 10;
		Gaussianlines.Drawsimulatedlines(imgout, range,rnd,sigma, numlines);
		
		ImageJFunctions.show(imgout);
		
		noisyimg = Poissonprocess.poissonProcess(imgout, 25f);
		//noisyimg = imgout;
		
		FloatType minval = new FloatType(0);
		FloatType maxval = new FloatType(1);
		Normalize.normalize(Views.iterable(noisyimg), minval, maxval);
		//Kernels.GaussianBlur(noisyimg);
		
		ImageJFunctions.show(noisyimg);
		
		
	}
	
}

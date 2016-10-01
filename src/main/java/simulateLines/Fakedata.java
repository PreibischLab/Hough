package simulateLines;

import java.io.File;
import java.util.Random;

import ij.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.Interval;
import net.imglib2.IterableInterval;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.stats.Normalize;
import net.imglib2.exception.IncompatibleTypeException;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.Type;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Intervals;
import net.imglib2.view.Views;
import poissonSimulator.Poissonprocess;
import preProcessing.Kernels;
import util.ImgLib2Util;

public class Fakedata {

	
	public static void main (String[] args) throws IncompatibleTypeException{
		
		new ImageJ();
		
		final FinalInterval range = new FinalInterval(600, 600);
		final FinalInterval smallrange = new FinalInterval(450, 450);
		RandomAccessibleInterval<FloatType> imgout = new ArrayImgFactory<FloatType>().create(range, new FloatType());
		RandomAccessibleInterval<FloatType> noisyimg = new ArrayImgFactory<FloatType>().create(imgout, new FloatType());
		RandomAccessibleInterval<FloatType> noisyimgsec = new ArrayImgFactory<FloatType>().create(imgout, new FloatType());
		final int ndims = imgout.numDimensions();
		
		FloatType minval = new FloatType(0);
		FloatType maxval = new FloatType(1);
		final double [] sigma = {1.4,1.5};
		final double [] Ci = new double[ndims];
		
		for (int d = 0; d < ndims; ++d)
			Ci[d] = 1.0 / Math.pow(sigma[d],2);
		
		
		final int numlines = 20;
		Gaussianlines.Drawsimulatedlines(imgout, smallrange,  sigma, numlines);
		Normalize.normalize(Views.iterable(imgout), minval, maxval);
		Kernels.addBackground(Views.iterable(imgout), 0.2);
		ImageJFunctions.show(imgout);
		
		noisyimg = Poissonprocess.poissonProcess(imgout, 45);
		noisyimgsec = Poissonprocess.poissonProcess(imgout, 35);
		//noisyimg = imgout;
		
		
		Normalize.normalize(Views.iterable(noisyimg), minval, maxval);
		Normalize.normalize(Views.iterable(noisyimgsec), minval, maxval);
		
		
		ImageJFunctions.show(noisyimg);
		ImageJFunctions.show(noisyimgsec);
		
	}
	
	
	


	
}

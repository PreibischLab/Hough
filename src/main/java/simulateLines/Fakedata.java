package simulateLines;

import java.util.Random;

import ij.ImageJ;
import net.imglib2.FinalInterval;
import net.imglib2.Interval;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Intervals;
import net.imglib2.view.Views;

public class Fakedata {

	
	public static void main (String[] args){
		
		new ImageJ();
		
		final FinalInterval range = new FinalInterval(200, 200);
		
		
		RandomAccessibleInterval<FloatType> imgout = new ArrayImgFactory<FloatType>().create(range, new FloatType());
		final int ndims = imgout.numDimensions();
		final Random rnd = new Random(150);
		final double [] sigma = {1.7,1.8};
		final double [] Ci = new double[ndims];
		
		for (int d = 0; d < ndims; ++d)
			Ci[d] = 1.0 / Math.pow(sigma[d],2);
		
		
		final int numlines = 4;
		Gaussianlines.Drawsimulatedlines(imgout, range,rnd,Ci, numlines);
		
		
		ImageJFunctions.show(imgout);
		
		
	}
	
}

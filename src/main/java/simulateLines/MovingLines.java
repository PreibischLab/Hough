package simulateLines;

import java.util.ArrayList;
import java.util.Random;

import ij.ImageJ;
import net.imglib2.FinalInterval;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.exception.IncompatibleTypeException;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import poissonSimulator.Poissonprocess;
import preProcessing.Kernels;

public class MovingLines {

	
	
	public static void main (String[] args) throws IncompatibleTypeException{
        new ImageJ();
		
		final FinalInterval range = new FinalInterval(512, 512);
		final FinalInterval smallrange = new FinalInterval(412, 412);
		
		
		
		final int ndims = range.numDimensions();
		final double [] sigma = {1.7,1.8};
		final double [] Ci = new double[ndims];
		
		for (int d = 0; d < ndims; ++d)
			Ci[d] = 1.0 / Math.pow(sigma[d],2);
		
		final int numframes = 15;
		final Random rnd = new Random(130);
		final int numlines = 10;

		ArrayList<Fakeline> linearrayini = Gaussianlines.Drawstartline(smallrange, rnd, sigma, numlines);
		for (int frame = 0; frame < numframes; ++frame){
			
		RandomAccessibleInterval<FloatType> lineimage = new ArrayImgFactory<FloatType>().create(range, new FloatType());
		RandomAccessibleInterval<FloatType> noisylines = new ArrayImgFactory<FloatType>().create(range, new FloatType());
		Kernels.SaltandPepperNoise(lineimage);
		
		Gaussianlines.Drawmovingsimulatedlines(lineimage, smallrange, frame, linearrayini, sigma);
	
		noisylines = Poissonprocess.poissonProcess(lineimage, 25);
		
		ImageJFunctions.show(noisylines);
		
		
		}
	}
}

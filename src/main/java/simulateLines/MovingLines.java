package simulateLines;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;

import com.sun.tools.javac.util.Pair;

import ij.ImageJ;
import mISC.Tree.Distance;
import net.imglib2.FinalInterval;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.stats.Normalize;
import net.imglib2.exception.IncompatibleTypeException;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import poissonSimulator.Poissonprocess;
import preProcessing.Kernels;

public class MovingLines {

	
	
	public static void main (String[] args) throws IncompatibleTypeException, IOException{
        new ImageJ();
		
		final FinalInterval range = new FinalInterval(1212, 1212);
		final FinalInterval smallrange = new FinalInterval(912, 912);
		
		
		
		final int ndims = range.numDimensions();
		final double [] sigma = {1.65,1.47};
		final double [] Ci = new double[ndims];
		
		for (int d = 0; d < ndims; ++d)
			Ci[d] = 1.0 / Math.pow(sigma[d],2);
		
		final int numframes = 20;
		final int numlines = 5;

		
		
			
		RandomAccessibleInterval<FloatType> lineimage = new ArrayImgFactory<FloatType>().create(range, new FloatType());
		RandomAccessibleInterval<FloatType> noisylines = new ArrayImgFactory<FloatType>().create(range, new FloatType());
		
		
		ArrayList<double[]> startseeds = new ArrayList<double[]>();
		ArrayList<double[]> endseeds = new ArrayList<double[]>();
		Gaussianlines.GetSeeds(lineimage, startseeds, endseeds, smallrange, numlines, sigma);
		
		
		

		FloatType minval = new FloatType(0);
		FloatType maxval = new FloatType(1);
		Normalize.normalize(Views.iterable(lineimage), minval, maxval);
		Kernels.addBackground(Views.iterable(lineimage), 0.2);
		noisylines = Poissonprocess.poissonProcess(lineimage, 20);
		ImageJFunctions.show(noisylines);
		ArrayList<double[]> startseedscopy =  new ArrayList<double[]>();
		ArrayList<double[]> endseedscopy = new ArrayList<double[]>();
		for (int index = 0; index < startseeds.size(); ++index){
			
			startseedscopy.add(index, startseeds.get(index));
			
		}
       for (int index = 0; index < endseeds.size(); ++index){
			
			endseedscopy.add(index, endseeds.get(index));
			
		}
		
		for (int frame = 1; frame < numframes; ++frame){
			
			RandomAccessibleInterval<FloatType> noisylinesframe = new ArrayImgFactory<FloatType>().create(range, new FloatType());
			RandomAccessibleInterval<FloatType> lineimageframe = new ArrayImgFactory<FloatType>().create(range, new FloatType());
			
			Pair<ArrayList<double[]>, ArrayList<double[]>> pair	 = Gaussianlines.Growseeds(lineimageframe, startseeds, endseeds, frame, sigma);
		
			double[] prevst = new double[ndims];
			double[] nextst = new double[ndims];
			double[] preven = new double[ndims];
			double[] nexten = new double[ndims];
			for (int d = 0; d < ndims; ++d){
				prevst[d] = pair.fst.get(0)[d];
				nextst[d] = pair.fst.get(1)[d];
				
				preven[d] = pair.snd.get(0)[d];
				nexten[d] = pair.snd.get(1)[d];
				
				
			}
			
	    	double startlength = Distance(prevst, nextst);
	    	double endlength = Distance(preven, nexten);
	
	   
		
		Normalize.normalize(Views.iterable(lineimageframe), minval, maxval);
		Kernels.addBackground(Views.iterable(lineimageframe), 0.2);
		noisylinesframe = Poissonprocess.poissonProcess(lineimageframe, 20);
		
	
		ImageJFunctions.show(noisylinesframe);
		
		
		
		
		
		FileWriter writerend = new FileWriter("../res/Actuallength-movingend.txt", true);
		   FileWriter writerstart = new FileWriter("../res/Actuallength-movingstart.txt", true);
		
			
			writerend.write( frame + " " + endlength);
			writerend.write("\r\n");
		
	
		
    
		
			
			writerstart.write( frame + " " + startlength );
			writerstart.write("\r\n");
			writerend.close();
			writerstart.close();
		
	
		}
	
		
	}
	public static double Distance(final double[] cordone, final double[] cordtwo) {

		double distance = 0;

		for (int d = 0; d < cordone.length; ++d) {

			distance += Math.pow((cordone[d] - cordtwo[d]), 2);

		}
		return Math.sqrt(distance);
	}
}

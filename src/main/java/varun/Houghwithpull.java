package varun;

import java.io.File;
import java.io.FileNotFoundException;

import ij.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import util.ImgLib2Util;

public class Houghwithpull {

	public static <T extends RealType<T>> void Houghspace(Img<T> inputimage, RandomAccessibleInterval<T> imgout,
			double[] min, double[] max) {

		int n = inputimage.numDimensions();
		final double[] inputposition = new double[n];
		final double[] position = new double[n];
		final long[] secondposition = new long[n];
		double sigmasq = 0, sigma, cutoff;
		double[] delta = new double[n];
		double[] actualposition = new double[n];
		final double[] realpos = new double[n];
		final double[] secondrealpos = new double[n];
		final Cursor<T> inputcursor = inputimage.localizingCursor();

		
		for (int d = 0; d < n; ++d) {
			delta[d] = (max[d] - min[d]) / (imgout.dimension(d));

			sigmasq += Math.pow(2*delta[d], 2);

		}
		sigma = Math.sqrt(sigmasq);
		cutoff = 5 * sigma;
		
		// for every function (as defined by an individual pixel)
		//while (inputcursor.hasNext()) {
		for (int i =0; i<10; ++i){
					inputcursor.fwd();
					inputcursor.localize(inputposition);

			final Cursor<T> outputcursor = Views.iterable(imgout).localizingCursor();

			final RandomAccess<T> outbound = imgout.randomAccess();
			
			
			while (outputcursor.hasNext()) {
				outputcursor.fwd();
				outputcursor.localize(position);
				final Cursor<T> second = Views.iterable(imgout).localizingCursor();
				
				// Forward transformation
				for (int d = 0; d < n; ++d) 
					realpos[d] = position[d] * delta[d] + min[d];
				
				double distance = 0;
				double intensity = 0;
				outbound.setPosition(outputcursor);
				double mindistance = Double.MAX_VALUE;
				double secmindistance = Double.MAX_VALUE;
				// To set the pixel intensity to the distance from the curve 
				while (second.hasNext()) {
					second.fwd();
					second.localize(secondposition);
					
					// Forward transformation
					for (int d = 0; d < n; ++d) 
						secondrealpos[d] = secondposition[d] * delta[d] + min[d];
		         
					/** Sin function, x=secondrealpos  A*sin(B*x+C) **/
					double A = Math.pow(inputposition[0], 2) + Math.pow(inputposition[1], 2);
					double C;
					
					C = Math.atan2(inputposition[0],inputposition[1]);
					
					
					Finalfunction sinfunction = new Finalfunction(secondrealpos, Math.sqrt(A), 1, Math.toDegrees(C));
					final double functionvalue = sinfunction.Sinfunction();
					final double functionderiv = sinfunction.DerivSinfunction();
					
				// This gives the distance from realpos to the Normal line constructed at point secondrealpos
				final double distanceline = Finaldistance.disttocurve(secondrealpos, realpos, functionvalue, functionderiv);
				
				if (distanceline <= mindistance) {
					mindistance = distanceline;
					actualposition[0] = secondrealpos[0];
					actualposition[1] = functionvalue;
					distance = Finaldistance.Generalfunctiondist(actualposition, realpos);
					if (distance <= secmindistance)
						secmindistance = distance;
				}
				
				intensity += (1 / (sigma * Math.sqrt(2 * Math.PI))) * Math.exp(-secmindistance * secmindistance / (2 * sigmasq));
				outbound.get().setReal(intensity);	
				
					
				}
			
					
		}
					
		
		}
	}
	
	public static void main(String[] args) throws FileNotFoundException {

		final Img<FloatType> inputimg = ImgLib2Util.openAs32Bit(new File("src/main/resources/small_line.tif"));
		// ImageJFunctions.show(inputimg);
		double size = Math.sqrt(Math.pow(inputimg.dimension(0), 2) + Math.pow(inputimg.dimension(1), 2));
		
		int maxRho = (int) Math.round( size ); 
		
		double[] min = { 0, -maxRho };
		double[] max = { 180, maxRho };

		

		final Img<FloatType> houghimage = new ArrayImgFactory<FloatType>().create(inputimg,
				new FloatType());

		Houghspace(inputimg, houghimage, min, max);

		new ImageJ();

		ImageJFunctions.show(houghimage);
	}

}

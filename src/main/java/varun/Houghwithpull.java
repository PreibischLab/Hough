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
		while (inputcursor.hasNext()) {
		
					inputcursor.fwd();
					inputcursor.localize(inputposition);

			final Cursor<T> outputcursor = Views.iterable(imgout).localizingCursor();

			final RandomAccess<T> outbound = imgout.randomAccess();
			
			final RandomAccess<T> poly = imgout.randomAccess();
			double[] gradient = new double[n];
			
			while (outputcursor.hasNext()) {
				outputcursor.fwd();
				outputcursor.localize(position);
				
				// Forward transformation
				for (int d = 0; d < n; ++d) 
					realpos[d] = position[d] * delta[d] + min[d];
				
				double distance = 0;
				double intensity = 0;
				outbound.setPosition(outputcursor);
				double mindistance = Double.MAX_VALUE;
				double secmindistance = Double.MAX_VALUE;
				double t = min[0];
				double step = 0.5;
				
				while (true) {
		
					
					/** Sin function, x=secondrealpos  A*sin(x+C) **/
					double A = Math.pow(inputposition[0], 2) + Math.pow(inputposition[1], 2);
					double C;
					
					C = Math.atan2(inputposition[0],inputposition[1]);
					
					  Finalfunction function = new Finalfunction(new double[] {t,0}, A, C); 
					  final double functionvalue = function.Sinfunction(); 
					  final double functionderiv = function.DerivSinfunction();
					 
					gradient[0] = 1;
					gradient[1] = functionderiv;

					

					poly.setPosition(Math.round(t + gradient[0]), 0);
					poly.setPosition(Math.round(functionvalue + gradient[1]), 1);

					double newx = poly.getDoublePosition(0);
					double newy = poly.getDoublePosition(1);
					// If slope of the tangent is greater than 10 degrees, go slow
					if (Math.abs(Math.toDegrees(Math.atan(functionvalue/t)))>5)
						step = 0.02;

					final double distanceline = Finaldistance.disttocurve(new double[] { newx, newy }, realpos, newy,
							functionderiv);
					if (distanceline <= mindistance) {
						mindistance = distanceline;
						actualposition[0] = newx;
						actualposition[1] = newy;
						distance = Finaldistance.Generalfunctiondist(actualposition, realpos);
						if (distance < secmindistance)
							secmindistance = distance;
					}

					intensity = (1 / (sigma * Math.sqrt(2 * Math.PI)))
							* Math.exp(-secmindistance * secmindistance / (2 * sigmasq));
					outbound.get().setReal(intensity);
					t += step;

					if (t >= max[0])
						break;
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

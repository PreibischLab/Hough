package varun;

import java.io.File;

import ij.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
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

public class PlotFunctionPush {
	

	public static <T extends RealType<T>> void push(RandomAccessibleInterval<T> imgout, double[] min,
			double[] max) {

		int n = imgout.numDimensions();
		final int[] position = new int[n];

		final int[] realposone = new int[n];
		final int[] realpostwo = new int[n];
		final int[] realposthree = new int[n];
		final int[] realposfour = new int[n];
		final int[] realposfive = new int[n];
		final double[] realpos = new double[n];
		double functionone, functiontwo, functionthree, functionfour, functionfive;
		int newposone, newpostwo, newposthree, newposfour, newposfive;
		final double sigma = 100.0;
		final double radius = 100;
		final double radiussecond = 200;

		double[] delta = new double[n];

		double[] centercircle = new double[n]; // Defines the center of the circle 1
		
		double[] centersecond = new double[n];
		
		
		for (int d = 0; d < n; ++d) {
			delta[d] = (max[d] - min[d]) / imgout.dimension(d);
			centercircle[d] = 0; // (imgout.dimension(d) / 2) * delta[d] + min[d];
			centersecond[d] = 0; // (3*imgout.dimension(d) / 4) * delta[d] + min[d];
		}

		
		// For centering Gaussian and Quadratic function at the center of the image
		
		double center = (imgout.dimension(0) / 2) * delta[0] + min[0]; 
		
		final Cursor<T> inputcursor = Views.iterable(imgout).localizingCursor();
		final RandomAccess<T> outbound = imgout.randomAccess();

		while (inputcursor.hasNext()) {
			inputcursor.fwd();
			inputcursor.localize(position);

			// Forward transformation
			
			for (int d = 0; d < n; ++d) {

				realpos[d] = inputcursor.getDoublePosition(d) * delta[d] + min[0];

			}
			// Gaussian function

			functionone = -100 * Math.exp(-Math.pow((realpos[0] - center), 2) / sigma); 

			// Sin function
			
			functiontwo = 10 * Math.sin(Math.toRadians(realpos[0] * 10)); 
			
			// Quadratic function

			functionthree = (Math.pow((realpos[0] - center), 2) ); 
			
			// Half circle function
			
			functionfour = -Math.sqrt(radius*radius-(realpos[0]-centercircle[0])*(realpos[0]-centercircle[0]))+centercircle[1];
			
			// Second Half circle
			
			functionfive = -Math.sqrt(radiussecond * radiussecond - (realpos[0] - centersecond[0]) * (realpos[0] - centersecond[0]))
					+ centersecond[1];
																			 

			// Back transformation
			
			//Transforming back the y coordinate for Gaussian
			
			newposone = (int) Math.round((functionone - min[1]) / delta[1]); 
			
			//Transforming back the y coordinate for Sin function
			
			newpostwo = (int) Math.round((functiontwo - min[1]) / delta[1]); 
			
			//Transforming back the y coordinate for Quadratic function
			
			newposthree = (int) Math.round((functionthree - min[1]) / delta[1]);
			
            //Transforming back the y coordinate for Circle function
			
			newposfour = (int) Math.round((functionfour - min[1]) / delta[1]);
			
           //Transforming back the y coordinate for second Circle function
			
			newposfive = (int) Math.round((functionfive - min[1]) / delta[1]);
			
			

			// To get the new pixels, first set them to original values and then
			// change the y value to the new one

			inputcursor.localize(realposone);
			inputcursor.localize(realpostwo);
			inputcursor.localize(realposthree);
			inputcursor.localize(realposfour);
			inputcursor.localize(realposfive);

			realposone[1] = newposone;
			realpostwo[1] = newpostwo;
			realposthree[1] = newposthree;
			realposfour[1] = newposfour;
			realposfive[1] = newposfive;
/*
			// To make sure that the new y values are in the proper range, plot
			// only when in range

			// Plotting for function one
			if (newposone < imgout.dimension(1)) {

				outbound.setPosition(realposone);

				outbound.get().setReal(1);

			}

			// Plotting for function two
			if (newpostwo < imgout.dimension(1)) {

				outbound.setPosition(realpostwo);

				outbound.get().setReal(1);

			}

			// Plotting for function three
			if (newposthree < imgout.dimension(1)) {

				outbound.setPosition(realposthree);

				outbound.get().setReal(1);

			}
		*/	
			// Plotting for function four
						if (newposfour < imgout.dimension(1)) {

							outbound.setPosition(realposfour);

							outbound.get().setReal(1);

						}
						
			// Plotting for function five
			if (newposfive < imgout.dimension(1)) {

				outbound.setPosition(realposfive);

				outbound.get().setReal(1);

			}
			

		}
	}

	public static void main(String[] args) {

		

	

		double[] min = { -400, -500 };
		double[] max = { 400, 0 };

		final double ratio = (max[1]-min[1]) / (max[0]-min[0]);
		final int sizeX = 1000;
		final int sizeY =  (int)Math.round( sizeX * ratio ); 


		final Img<FloatType> houghquadimage = new ArrayImgFactory<FloatType>().create(new long[]{sizeX, sizeY}, new FloatType());

		push(houghquadimage, min, max);
  //new ImageJ();
		ImageJFunctions.show(houghquadimage).setTitle("Push-circle function");

	}
}

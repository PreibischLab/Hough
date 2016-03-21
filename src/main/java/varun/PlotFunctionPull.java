package varun;

import java.io.File;
import java.io.FileNotFoundException;

import javax.annotation.Generated;

import ij.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.Localizable;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.neighborsearch.NearestNeighborSearch;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import util.ImgLib2Util;

public class PlotFunctionPull {
	// function: (x-x0)^2 +(y-y0)^2 = R^2

	// for every pixel in the output image, compute the distance to the
	// function
	// intensity = distance_function (circle)
	// %% in general intensity = sum_f gauss( distance_f ); f ... function

	public static <T extends RealType<T>> void pull(RandomAccessibleInterval<T> imgout, double[] min, double[] max) {

		int n = imgout.numDimensions();
		final int[] position = new int[n];

		final double[] realpos = new double[n];

		final double[] secondposition = new double[n];

		double[] actualposition = new double[n];

		final double[] secondrealpos = new double[n];

		double radius = 200, sigmasq = 0, sigma, cutoff;

		double[] delta = new double[n];

		double[] center = new double[n];

		for (int d = 0; d < n; ++d) {
			delta[d] = (max[d] - min[d]) / (imgout.dimension(d));
			center[d] = 0; // Center of the circle
			sigmasq += Math.pow(delta[d], 2);

		}

		// sigmasq = Math.pow(delta[0], 2) + Math.pow(delta[1], 2);
		sigma =  Math.sqrt(sigmasq); // Ensures resolution for small of big
											// size boxes
		cutoff = 5 * sigma;

		final Cursor<T> inputcursor = Views.iterable(imgout).localizingCursor();

		final RandomAccess<T> outbound = imgout.randomAccess();

		while (inputcursor.hasNext()) {
			inputcursor.fwd();

			inputcursor.localize(position);
			final Cursor<T> second = Views.iterable(imgout).localizingCursor();

			// Forward transformation

			for (int d = 0; d < n; ++d) {

				realpos[d] = position[d] * delta[d] + min[d];

			}

			double mindistance = Double.MAX_VALUE;
			double distance = 0;
			/*
			 * 
			 * For circle and line with known distance formula
			 * 
			 * 
			 * // for the circle function (center and radius describes a circle)
			 * 
			 * Finalfunction circlefunction = new Finalfunction(realpos, center,
			 * radius, 0);
			 * 
			 * distance = circlefunction.Circlefunctiondist();
			 * 
			 * // for a line the slope and the constant along y axis describe
			 * the curve
			 * 
			 * Finalfunction linefunction = new Finalfunction(realpos,10,4);
			 * distancesecond = linefunction.Linefunctiondist();
			 * 
			 */
			outbound.setPosition(inputcursor);

			while (second.hasNext()) {
				second.fwd();
				second.localize(secondposition);

				// Forward transformation

				for (int d = 0; d < n; ++d) {

					secondrealpos[d] = secondposition[d] * delta[d] + min[d];

				}

				/** Sin function, x=secondrealpos  A*sin(B*x+C) **/
				
				Finalfunction sinfunction = new Finalfunction(secondrealpos, 10, 10, 4);
				final double functionvalue = sinfunction.Sinfunction();
				final double functionderiv = sinfunction.DerivSinfunction();
            
				
				/** Quadratic function, x=secondrealpos, A*x^2+B*x+C **/
			/*	
				Finalfunction Quadfunction = new Finalfunction(secondrealpos, 0.01,0.02, 4);
				final double functionvalue = Quadfunction.Quadfunction();
				final double functionderiv = Quadfunction.DerivQuadfunction();
				
			*/
				/** Cubic function, x=secondrealpos, A*x^3+B*x^2+C*x+D **/
			/*	
				Finalfunction Cubicfunction = new Finalfunction(secondrealpos, 1.1,0.1,-4.1, 4);
				final double functionvalue = Cubicfunction.Cubicfunction();
				final double functionderiv = Cubicfunction.DerivCubicfunction();
				
				/** Biquad function, x=secondrealpos, A*x^4+B*x^3+C*x^2+D*x+E **/
			/*	
				Finalfunction Biquadfunction = new Finalfunction(secondrealpos,2.1,-2.4,-2.4,1.8,4.0);
				final double functionvalue = Biquadfunction.Biquadfunction();
				final double functionderiv = Biquadfunction.DerivBiquadfunction();
			*/	
				Finalfunction Normalline = new Finalfunction(realpos, -1.0 / functionderiv,
						secondrealpos[0] / functionderiv + functionvalue);
				final double distanceline = Normalline.Linefunctiondist();

				if (distanceline <= mindistance) {

					mindistance = distanceline;

					actualposition[0] = secondrealpos[0];
					actualposition[1] = functionvalue;

					distance = Finaldistance.Generalfunctiondist(actualposition, realpos);
				}

			}

			if (Math.abs(distance) < cutoff)

			//	outbound.get().setReal(distance);
				outbound.get().setReal(Math.exp(-distance * distance / sigmasq));

			else

				outbound.get().setReal(0);

		}

	}

	public static void main(String[] args) throws FileNotFoundException {

		double[] min = { -60, -20 };
		double[] max = { 60, 20 };

		final double ratio = (max[1] - min[1]) / (max[0] - min[0]);
		final int sizeX = 350;
		final int sizeY = (int) Math.round(sizeX * ratio);

		final Img<FloatType> houghquadimage = new ArrayImgFactory<FloatType>().create(new long[] { sizeX, sizeY },
				new FloatType());

		pull(houghquadimage, min, max);

		 new ImageJ();
		ImageJFunctions.show(houghquadimage).setTitle("Pull function");

	}
}

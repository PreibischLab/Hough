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

		sigma = 0.7 * Math.sqrt(sigmasq); // Ensures resolution for small of big
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

				Finalfunction sinfunction = new Finalfunction(secondrealpos, 10, 10, 4);

				final double functionvalue = sinfunction.Sinfunction();
				final double functionderiv = sinfunction.DerivSinfunction();

				Finalfunction Normalline = new Finalfunction(realpos, -1.0 / functionderiv,
						secondrealpos[0] / functionderiv + functionvalue);
				final double distanceline = Normalline.Linefunctiondist();

				if (distanceline <= mindistance) {

					mindistance = distanceline;

					actualposition[0] = secondrealpos[0];
					actualposition[1] = functionvalue;

					distance = Finalfunction.Generalfunctiondist(actualposition, realpos);
				}

			}

			if (Math.abs(distance) < cutoff)

				outbound.get().setReal(Math.exp(-distance * distance / sigmasq));

			else

				outbound.get().setReal(0);

		}

	}

	public static void main(String[] args) throws FileNotFoundException {

		double[] min = { -100, -100 };
		double[] max = { 100, 100 };

		final double ratio = (max[1] - min[1]) / (max[0] - min[0]);
		final int sizeX = 400;
		final int sizeY = (int) Math.round(sizeX * ratio);

		final Img<FloatType> houghquadimage = new ArrayImgFactory<FloatType>().create(new long[] { sizeX, sizeY },
				new FloatType());

		pull(houghquadimage, min, max);

		// new ImageJ();
		ImageJFunctions.show(houghquadimage).setTitle("Pull function");

	}
}

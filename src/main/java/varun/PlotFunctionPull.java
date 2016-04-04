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
import net.imglib2.roi.EllipseRegionOfInterest;
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

		final double[] position = new double[n];

		final double[] realpos = new double[n];

		final double[] secondposition = new double[n];
		final double[] secondrealpos = new double[n];
		double[] actualposition = new double[n];

		double sigmasq = 0, sigma;

		double[] delta = new double[n];

		for (int d = 0; d < n; ++d) {

			delta[d] = (max[d] - min[d]) / (imgout.dimension(d));

			sigmasq += Math.pow(delta[d], 2);

		}
		sigma =  Math.sqrt(sigmasq);
System.out.println(sigma);
		final Cursor<T> inputcursor = Views.iterable(imgout).localizingCursor();

		final RandomAccess<T> outbound = imgout.randomAccess();
		
		double[] center = { 0, 5 }; double radius = 10;
		while (inputcursor.hasNext()) {
			inputcursor.fwd();
			inputcursor.localize(position);
			final Cursor<T> second = Views.iterable(imgout).localizingCursor();

			// Forward transformation
			for (int d = 0; d < n; ++d)
				realpos[d] = position[d] * delta[d] + min[d];
			outbound.setPosition(inputcursor);
			double mindistance = Double.MAX_VALUE;
			double secmindistance = Double.MAX_VALUE;
			// To set the pixel intensity to the distance from the curve
			double distance = 0;
			double intensity = 0;

			Finalfunction circlefunction = new Finalfunction(realpos, center,
					  radius, 0);
					  
					  distance = circlefunction.Circlefunctiondist();
			
	/*		while (second.hasNext()) {
				second.fwd();
				second.localize(secondposition);

				// Forward transformation
				for (int d = 0; d < n; ++d)
					secondrealpos[d] = secondposition[d] * delta[d] + min[d];

				
				Finalfunction function = new Finalfunction(secondrealpos, 0.01, 0.0,0.4);
				final double functionvalue = function.Quadfunction();
				final double functionderiv = function.DerivQuadfunction();
				// This gives the distance from realpos to the Normal line
				// constructed at point secondrealpos
				final double distanceline = Finaldistance.disttocurve(secondrealpos, realpos, functionvalue,
						functionderiv);
				if (distanceline <= mindistance) {
					mindistance = distanceline;
					actualposition[0] = secondrealpos[0];
					actualposition[1] = functionvalue;
					distance = Finaldistance.Generalfunctiondist(actualposition, realpos);
					if (distance <= secmindistance)
						secmindistance = distance;
				}
				*/
				intensity = (1 / (sigma * Math.sqrt(2 * Math.PI))) * Math.exp(-distance * distance / (2 * sigmasq));
				outbound.get().setReal(intensity);
			
		}
	}

	public static void main(String[] args) throws FileNotFoundException {

		double[] min = { -40, -10 };
		double[] max = { 40, 40 };

		final double ratio = (max[1] - min[1]) / (max[0] - min[0]);
		final int sizeX = 800;
		final int sizeY = (int) Math.round(sizeX * ratio);

		final Img<FloatType> houghimage = new ArrayImgFactory<FloatType>().create(new long[] { sizeX, sizeY },
				new FloatType());
		long startTime = System.currentTimeMillis();
		pull(houghimage, min, max);
		long endTime = System.currentTimeMillis();
		long totalTime = endTime - startTime;
		System.out.println("Normal line finding time :" + totalTime);

		new ImageJ();
		ImageJFunctions.show(houghimage).setTitle("Exact distance formula");

	}
}

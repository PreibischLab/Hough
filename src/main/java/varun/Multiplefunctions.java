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

public class Multiplefunctions{

	public static <T extends RealType<T>> void pull(RandomAccessibleInterval<T> imgout, double[] min, double[] max) {

		int n = imgout.numDimensions();
		final int[] position = new int[n];
		final double[] realpos = new double[n];
		double sigmasq = 0, sigma, cutoff, distance = 0;
		double[] delta = new double[n];
		double[] center = new double[n];

		for (int d = 0; d < n; ++d) {
			delta[d] = (max[d] - min[d]) / (imgout.dimension(d));
			center[d] = 0; // Center of the circle
			sigmasq += Math.pow(delta[d], 2);

		}
		sigma = Math.sqrt(sigmasq);
		cutoff = 10 * sigma;
		final Cursor<T> inputcursor = Views.iterable(imgout).localizingCursor();
		final RandomAccess<T> outbound = imgout.randomAccess();
		
		
			
		while (inputcursor.hasNext()) {

			inputcursor.fwd();
			inputcursor.localize(position);

			// Forward transformation
			for (int d = 0; d < n; ++d) 
				realpos[d] = position[d] * delta[d] + min[d];
			
			double intensity = 0;
			// for (double radius = 20; radius < 100; radius += 20) {
			for (double slope = 10; slope < 100; slope += 30) {
				// Finalfunction circlefunction = new Finalfunction(realpos,
				// center, radius, 0);
				// distance = circlefunction.Circlefunctiondist();

				Finalfunction linefunction = new Finalfunction(realpos, slope, 4);
				distance = linefunction.Linefunctiondist();

				outbound.setPosition(inputcursor);

				if (Math.abs(distance) < cutoff)
					intensity += Math.exp(-distance * distance / sigmasq);
			}
			outbound.get().setReal(intensity);
		}
	}

	public static void main(String[] args) throws FileNotFoundException {

		double[] min = { -500, -100 };
		double[] max = { 500, 200 };

		final double ratio = (max[1] - min[1]) / (max[0] - min[0]);
		final int sizeX = 1000;
		final int sizeY = (int) Math.round(sizeX * ratio);

		final Img<FloatType> houghquadimage = new ArrayImgFactory<FloatType>().create(new long[] { sizeX, sizeY },
				new FloatType());

		pull(houghquadimage, min, max);

		new ImageJ();
		ImageJFunctions.show(houghquadimage).setTitle("Pull function");

	}
}

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

public class ExactdisttoCircle {
	// function: (x-x0)^2 +(y-y0)^2 = R^2

	public static <T extends RealType<T>> void disttocircle(RandomAccessibleInterval<T> imgout, double[] min, double[] max) {

		int n = imgout.numDimensions();

		final double[] position = new double[n];

		final double[] realpos = new double[n];

		double sigmasq = 0, sigma;

		double[] delta = new double[n];

		for (int d = 0; d < n; ++d) {

			delta[d] = (max[d] - min[d]) / (imgout.dimension(d));

			sigmasq += Math.pow(delta[d], 2);

		}
		sigma = Math.sqrt(sigmasq);
		System.out.println(sigma);
		final Cursor<T> inputcursor = Views.iterable(imgout).localizingCursor();

		final RandomAccess<T> outbound = imgout.randomAccess();

		double[] center = { 0, 5 };
		double radius = 10;
		while (inputcursor.hasNext()) {
			inputcursor.fwd();
			inputcursor.localize(position);

			// Forward transformation
			for (int d = 0; d < n; ++d)
				realpos[d] = position[d] * delta[d] + min[d];
			outbound.setPosition(inputcursor);
			// To set the pixel intensity to the distance from the curve
			double distance = 0;
			double intensity = 0;

			Finalfunction circlefunction = new Finalfunction(realpos, center, radius, 0);

			distance = circlefunction.Circlefunctiondist();

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

        disttocircle(houghimage, min, max);

		new ImageJ();
		ImageJFunctions.show(houghimage).setTitle("Exact distance formula");

	}
}

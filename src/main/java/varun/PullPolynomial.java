package varun;

import java.io.FileNotFoundException;
import ij.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.Localizable;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealCursor;
import net.imglib2.RealLocalizable;
import net.imglib2.RealPointSampleList;
import net.imglib2.RealRandomAccess;
import net.imglib2.RealRandomAccessible;
import net.imglib2.algorithm.region.hypersphere.HyperSphereCursor;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;

public class PullPolynomial {

	public static <T extends RealType<T>> void pull(RandomAccessibleInterval<T> imgout, double[] min, double[] max) {

		int n = imgout.numDimensions();

		final double[] position = new double[n];

		final double[] realpos = new double[n];

		double[] actualposition = new double[n];

		double sigmasq = 0, sigma;

		double[] delta = new double[n];

		for (int d = 0; d < n; ++d) {

			delta[d] = (max[d] - min[d]) / (imgout.dimension(d));

			sigmasq +=  0.5*Math.pow(delta[d], 2);

		}

		sigma = Math.sqrt(sigmasq);

		System.out.println(sigma);
		double ampl = -40;
		final Cursor<T> inputcursor = Views.iterable(imgout).localizingCursor();
		final RandomAccess<T> outbound = imgout.randomAccess();
		final RandomAccess<T> poly = imgout.randomAccess();
		double[] gradient = { 1, 1 };
		while (inputcursor.hasNext()) {
			inputcursor.fwd();
			inputcursor.localize(position);

			for (int d = 0; d < n; ++d)
				realpos[d] = position[d] * delta[d] + min[d];

			outbound.setPosition(inputcursor);

			double mindistance = Double.MAX_VALUE;
			double secmindistance = Double.MAX_VALUE;

			double distance = 0;
			double intensity = 0;
			double step = 0.1;
			double t = min[0];

			while (true) {
				// Gauss function //
				Finalfunction function = new Finalfunction(new double[] { t, 0 }, ampl, 50);
				final double functionvalue = function.Gaussianfunction();
				final double functionderiv = function.DerivGaussianfunction();

				// Sin function //
				/*
				 * Finalfunction function = new Finalfunction(new double[]
				 * {t,0}, ampl, 0); final double functionvalue =
				 * function.Sinfunction(); final double functionderiv =
				 * function.DerivSinfunction();
				 */
				gradient[0] = 1;
				gradient[1] = functionderiv;

				

				poly.setPosition(Math.round(t + gradient[0]), 0);
				poly.setPosition(Math.round(functionvalue + gradient[1]), 1);

				double newx = poly.getDoublePosition(0);
				double newy = poly.getDoublePosition(1);
				// If slope of the tangent is greater than 10 degrees, go slow
				if (Math.abs(Math.toDegrees(Math.atan(functionvalue/t)))>10)
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

	public static void main(String[] args) throws FileNotFoundException {

		double[] min = { -100, -40 };
		double[] max = { 100, 40 };

		final double ratio = (max[1] - min[1]) / (max[0] - min[0]);
		final long sizeX = 200;
		final long sizeY = Math.round(sizeX * ratio);

		final Img<FloatType> houghimage = new ArrayImgFactory<FloatType>().create(new long[] { sizeX, sizeY },
				new FloatType());
		long startTime = System.currentTimeMillis();
		pull(houghimage, min, max);
		long endTime = System.currentTimeMillis();
		long totalTime = endTime - startTime;
		System.out.println("Normal line finding time :" + totalTime);

		// new ImageJ();
		ImageJFunctions.show(houghimage).setTitle("Pull Normal line finding function");

	}
}


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
import net.imglib2.RealRandomAccessibleRealInterval;
import net.imglib2.algorithm.region.hypersphere.HyperSphereCursor;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;

public class Circle {

	public static <T extends RealType<T>> void pull(RandomAccessibleInterval<T> imgout, double[] min, double[] max) {

		int n = imgout.numDimensions();
		final double[] inputposition = new double[n];
		final double[] position = new double[n];
		final double[] realpos = new double[n];
		double[] actualposition = new double[n];
		double sigmasq = 0, sigma;
		double[] delta = new double[n];

		for (int d = 0; d < n; ++d) {
			delta[d] = (max[d] - min[d]) / (imgout.realMax(d) - imgout.realMin(d) + 1);
			sigmasq += Math.pow(delta[d], 2);
		}

		sigma = Math.sqrt(sigmasq);
		System.out.println(sigma);

		double[] center = { 0, 0 };
		double radius = 40;
		final Cursor<T> inputcursor = Views.iterable(imgout).localizingCursor();
		final RandomAccess<T> outbound = imgout.randomAccess();
		final RandomAccess<T> circle = imgout.randomAccess();
		double[] gradient = new double[n];
		// Select a point and assign its pixel value to the distance from the
		// curve
		while (inputcursor.hasNext()) {
			inputcursor.fwd();
			inputcursor.localize(inputposition);
			for (int d = 0; d < n; ++d)
				realpos[d] = inputposition[d] * delta[d] + min[d];
			outbound.setPosition(inputcursor);

			double mindistance = Double.MAX_VALUE;
			double secmindistance = Double.MAX_VALUE;

			double distance = 0;
			double intensity = 0;
			double initheta = 0;
			double theta = initheta;
			double dtheta = 1;
			double[] newpos = new double[n];
			// Compute the Normal to the curve on which the input point lies
			while (true) {

				// Initialize a point on the circle
				position[0] = radius * Math.cos(Math.toRadians(initheta));
				position[1] = radius * Math.sin(Math.toRadians(initheta));

				for (int d = 0; d < n; ++d) {
					circle.setPosition(Math.round(center[d] + position[d]), d);
					newpos[d] = circle.getDoublePosition(d);
				}
				// Compute the gradient for incrementing current position
				
				theta = initheta + dtheta;
				gradient[0] = -radius * Math.sin(Math.toRadians(theta)) * Math.toRadians(dtheta);
				gradient[1] = radius * Math.cos(Math.toRadians(theta)) * Math.toRadians(dtheta);

				// Move the current point along the gradient
				for (int d = 0; d < n; ++d)
					circle.move(Math.round(gradient[d]), d);
				// Get the updated position of the point along the circle
				for (int d = 0; d < n; ++d) {
					circle.setPosition(Math.round(center[d] + position[d] + gradient[d]), d);
					newpos[d] = circle.getDoublePosition(d);
				}

				// Compute the Normal line at the point
				double distanceline = Finaldistance.disttocurve(new double[] { newpos[0], newpos[1] }, realpos,
						newpos[1], -(newpos[0] - center[0]) / (newpos[1] - center[1]));
				// If the point lies on the Normal compute distance bw it and
				// input point
				if (distanceline <= mindistance) {
					mindistance = distanceline;
					for (int d = 0; d < n; ++d)
						actualposition[d] = newpos[d];
					distance = Finaldistance.Generalfunctiondist(actualposition, realpos);
					// Store current distance and for multiple normal case get
					// the closest distance
					if (distance < secmindistance)
						secmindistance = distance;
				}

				intensity = (1 / (sigma * Math.sqrt(2 * Math.PI)))
						* Math.exp(-secmindistance * secmindistance / (2 * sigmasq));
				outbound.get().setReal(intensity);

				initheta = theta;

				if (theta >= 360)
					break;

			}
		}

	}

	public static void main(String[] args) throws FileNotFoundException {

		double[] min = { -100, -100 };
		double[] max = { 100, 100 };

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

		new ImageJ();
		ImageJFunctions.show(houghimage).setTitle("Pull Normal line finding function");

	}
}


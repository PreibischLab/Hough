package varun;

import com.sun.tools.javah.Util.Exit;

import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import varun.GetLocalmaxmin.IntensityType;

public class PushCurves {
	public static void drawCircle(Img<FloatType> imgout, double[] min, double[] max, double[] center, double radius) {
		int n = imgout.numDimensions();
		double[] realpos = new double[n];
		double[] size = new double[n];
		double[] location = new double[n];
		double[] position = new double[n];
		double[] iniposition = new double[n];
		double[] backini = new double[n];
		double[] newpos = new double[n];
		double[] backpos = new double[n];
		double[] sigma = new double[n];
		final RandomAccess<FloatType> outbound = imgout.randomAccess();
		double stepsize = 0.1;
		int[] setpos = new int[n];
		for (int d = 0; d < n; ++d)
			size[d] = imgout.dimension(d);

		Cursor<FloatType> cursor = Views.iterable(imgout).localizingCursor();
		while (cursor.hasNext()) {
			cursor.fwd();
			cursor.localize(location);
			realpos = TransformCordinates.transformfwd(location, size, min, max);

			// To get a starting point on the circle
			if (Math.pow(realpos[0] - center[0], 2) + Math.pow(realpos[1] - center[1], 2)
					- radius * radius <= 1.0E-50) {
				for (int d = 0; d < n; ++d)
					position[d] = realpos[d];
				break;

			}

		}

		for (int d = 0; d < n; ++d)
			iniposition[d] = position[d];

		double initheta = Math.atan2(iniposition[1] - center[1], iniposition[0] - center[0]);
		double increment = Math.acos((2 * radius * radius - stepsize * stepsize) / (2 * radius * radius));

		backini = TransformCordinates.transformback(iniposition, size, min, max);
		sigma[0] = 1;
		sigma[1] = 1;
		while (true) {

			// Move the current point along the curve

			newpos[0] = center[0] + radius * Math.cos((initheta - increment));
			newpos[1] = center[1] + radius * Math.sin((initheta - increment));
			initheta = Math.atan2(newpos[1] - center[1], newpos[0] - center[0]);

			// Transform the co-ordinates back as double[]
			backpos = TransformCordinates.transformback(newpos, size, min, max);

			setpos[0] = (int) Math.round(backpos[0]);
			setpos[1] = (int) Math.round(backpos[1]);

			// To set the pixel intensity
			AddGaussian.addGaussian(imgout, backpos, sigma, false);

			// To make sure that the values transformed back are not out of
			// bounds
			if (backpos[0] < imgout.realMax(0) - imgout.realMin(0) || backpos[0] > imgout.realMin(0)
					|| backpos[1] < imgout.realMax(1) - imgout.realMin(1) || backpos[1] > imgout.realMin(1))

				outbound.setPosition(setpos);

			// Stopping criteria of moving along the circular arc
			if (Math.abs(setpos[0] - (int) Math.round(backini[0])) == 0
					&& Math.abs(setpos[1] - (int) Math.round(backini[1])) == 0)
				break;

			// General Stopping criteria of moving along a curve, when we hit a
			// boundary
			if (newpos[0] >= max[0] || newpos[0] <= min[0] || newpos[1] >= max[1] || newpos[1] <= min[1])

				break;

		}
	}

	public static void DrawSine(RandomAccessibleInterval<FloatType> imgout, double[] min, double[] max,
			double amplitude, double phase) {

		int n = imgout.numDimensions();
		double[] size = new double[n];
		double[] position = new double[n];
		double[] newpos = new double[n];
		double[] backpos = new double[n];
		double[] sigma = new double[n];
		double increment;
		final RandomAccess<FloatType> outbound = imgout.randomAccess();
		double stepsize = 0.1;
		int[] setpos = new int[n];
		for (int d = 0; d < n; ++d)
			size[d] = imgout.dimension(d);

		// Starting position, for explicit curves its easier to choose a
		// starting point
		position[0] = min[0];
		position[1] = amplitude * Math.sin(Math.toRadians(position[0] + phase));
		newpos[0] = position[0];
		newpos[1] = position[1];
		sigma[0] = 1;
		sigma[1] = 1;
		while (true) {
			increment = stepsize * amplitude * Math.cos(Math.toRadians(position[0] + phase));

			for (int d = 0; d < n; ++d)
				position[d] = newpos[d];

			// Transform the co-ordinates back as double[]
			backpos = TransformCordinates.transformback(newpos, size, min, max);

			setpos[0] = (int) Math.round(backpos[0]);
			setpos[1] = (int) Math.round(backpos[1]);

			// To set the pixel intensity
			AddGaussian.addGaussian(imgout, backpos, sigma, false);

			// To make sure that the values transformed back are not out of
			// bounds
			if (backpos[0] < imgout.realMax(0) - imgout.realMin(0) || backpos[0] > imgout.realMin(0)
					|| backpos[1] < imgout.realMax(1) - imgout.realMin(1) || backpos[1] > imgout.realMin(1))

				outbound.setPosition(setpos);
			// Increment from starting position (min) towards max

			newpos[0] = position[0] + stepsize;
			newpos[1] = amplitude * Math.sin(Math.toRadians(newpos[0] + phase));
			// General Stopping criteria of moving along a curve, when we hit a
			// boundary
			if (newpos[0] >= max[0] || newpos[0] <= min[0] || newpos[1] >= max[1] || newpos[1] <= min[1])

				break;

		}

	}

	public static void Drawexactcircle(RandomAccessibleInterval<FloatType> imgout, double[] center, double radius,
			double[] min, double[] max) {

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

		final Cursor<FloatType> inputcursor = Views.iterable(imgout).localizingCursor();

		final RandomAccess<FloatType> outbound = imgout.randomAccess();

		while (inputcursor.hasNext()) {
			inputcursor.fwd();
			inputcursor.localize(position);

			// Forward transformation
			for (int d = 0; d < n; ++d)
				realpos[d] = position[d] * delta[d] + min[d];
			outbound.setPosition(inputcursor);
			// To set the pixel intensity as the shortest distance to the curve
			double distance = 0;
			double intensity = 0;

			Finalfunction circlefunction = new Finalfunction(realpos, center, radius, 0);

			distance = circlefunction.Circlefunctiondist();

			intensity = (1 / (sigma * Math.sqrt(2 * Math.PI))) * Math.exp(-distance * distance / (2 * sigmasq));
			outbound.get().setReal(intensity);

		}
	}

	public static void Drawexactline(RandomAccessibleInterval<FloatType> imgout, double slope, double intercept,
			final IntensityType setintensity) {

		int n = imgout.numDimensions();
		final double[] realpos = new double[n];
		double sigmasq, sigma = 0.5;
		sigmasq = sigma * sigma;
		sigma = Math.sqrt(sigmasq);
		final Cursor<FloatType> inputcursor = Views.iterable(imgout).localizingCursor();
		final RandomAccess<FloatType> outbound = imgout.randomAccess();

		while (inputcursor.hasNext()) {

			inputcursor.fwd();
			inputcursor.localize(realpos);

			// To set the pixel intensity as the shortest distance to the curve
			double distance = 0;
			double intensity = 0;

			Finalfunction linefunction = new Finalfunction(realpos, slope, intercept);
			distance = linefunction.Linefunctiondist();

			outbound.setPosition(inputcursor);

			intensity = (1 / (sigma * Math.sqrt(2 * Math.PI))) * Math.exp(-distance * distance / (2 * sigmasq));
			
			if (intensity < 1.0E-2)
				intensity = 0;
			
			switch (setintensity) {
			case Original:
				outbound.get().setReal(distance);
				break;

			case Gaussian:
				outbound.get().setReal(intensity);
				break;
			case One:
				outbound.get().setReal(intensity);
				break;
			default:
				outbound.get().setReal(intensity);
				break;

			}

		}
	}

}

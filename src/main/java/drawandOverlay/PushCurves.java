package drawandOverlay;

import java.util.ArrayList;
import java.util.Collections;

import houghandWatershed.Finalfunction;
import houghandWatershed.TransformCordinates;
import lut.SinCosinelut;
import net.imglib2.Cursor;
import net.imglib2.Point;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealPoint;
import net.imglib2.RealPointSampleList;
import net.imglib2.RealRandomAccess;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import peakFitter.LengthDetection;
import preProcessing.GetLocalmaxmin;
import preProcessing.GetLocalmaxmin.IntensityType;
import simulateLines.Fakeline;

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
			AddGaussian.addGaussian(imgout, backpos, sigma);

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
		//SinCosinelut.getTable();
		double stepsize = 0.1;
		int[] setpos = new int[n];
		for (int d = 0; d < n; ++d)
			size[d] = imgout.dimension(d);

		// Starting position, for explicit curves its easier to choose a
		// starting point
		// Input angles in degrees for the lut.
		position[0] = min[0];
		position[1] = //amplitude * SinCosinelut.getTable().getSine(position[0] + phase);
				amplitude * Math.sin(Math.toRadians(position[0] + phase));
		newpos[0] = position[0];
		newpos[1] = position[1];
		sigma[0] = 1;
		sigma[1] = 1;
		
		while (true) {
		//	increment = stepsize * amplitude * Math.cos(Math.toRadians(position[0] + phase));

			for (int d = 0; d < n; ++d)
				position[d] = newpos[d];
			
			// Transform the co-ordinates back as double[]
			backpos = TransformCordinates.transformback(newpos, size, min, max);

			setpos[0] = (int) Math.round(backpos[0]);
			setpos[1] = (int) Math.round(backpos[1]);

			// To set the pixel intensity
			AddGaussian.addGaussian(imgout, backpos, sigma);

			// To make sure that the values transformed back are not out of
			// bounds
			if (backpos[0] < imgout.realMax(0) - imgout.realMin(0) || backpos[0] > imgout.realMin(0)
					|| backpos[1] < imgout.realMax(1) - imgout.realMin(1) || backpos[1] > imgout.realMin(1))

				outbound.setPosition(setpos);
			// Increment from starting position (min) towards max
			
			newpos[0] = position[0] + stepsize;
			newpos[1] = //amplitude * SinCosinelut.getTable().getSine(newpos[0] + phase); 
					amplitude * Math.sin(Math.toRadians(newpos[0] + phase));
			// General Stopping criteria of moving along a curve, when we hit a
			// boundary
			if (newpos[0] >= max[0] || newpos[0] <= min[0] || newpos[1] >= max[1] || newpos[1] <= min[1])

				break;
		}

	}

	public static void DrawLine(RandomAccessibleInterval<FloatType> imgout, double[] min, double[] max, double slope,
			double intercept) {

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
		position[1] = slope * min[0] + intercept;
		newpos[0] = position[0];
		newpos[1] = position[1];
		sigma[0] = 1;
		sigma[1] = 1;
		while (true) {

			for (int d = 0; d < n; ++d)
				position[d] = newpos[d];

			// Transform the co-ordinates back as double[]
			backpos = TransformCordinates.transformback(newpos, size, min, max);

			setpos[0] = (int) Math.round(backpos[0]);
			setpos[1] = (int) Math.round(backpos[1]);

			// To set the pixel intensity
			AddGaussian.addGaussian(imgout, backpos, sigma);

			// To make sure that the values transformed back are not out of
			// bounds
			if (backpos[0] < imgout.realMax(0) - imgout.realMin(0) || backpos[0] > imgout.realMin(0)
					|| backpos[1] < imgout.realMax(1) - imgout.realMin(1) || backpos[1] > imgout.realMin(1))

				outbound.setPosition(setpos);
			// Increment from starting position (min) towards max

			newpos[0] = position[0] + stepsize / (1 + slope * slope);
			newpos[1] = position[1] + stepsize * slope / (1 + slope * slope);
			// General Stopping criteria of moving along a curve, when we hit a
			// boundary

			if (newpos[0] >= max[0] || newpos[0] <= min[0] || newpos[1] >= max[1] || newpos[1] <= min[1])

				break;

		}

	}
	public static void DrawDetectedGaussians(RandomAccessibleInterval<FloatType> imgout,
			final ArrayList<double[]> parameters) {

		final int n = imgout.numDimensions();
		ArrayList<Double> Amplitudelist = new ArrayList<Double>();
		ArrayList<double[]> Meanlist = new ArrayList<double[]>();
		ArrayList<double[]> Sigmalist = new ArrayList<double[]>();
		for (int index = 0; index < parameters.size(); ++index) {

			final double Amplitude = parameters.get(index)[0];
			final double Mean[] = new double[n];
			final double Sigma[] = new double[n];

			for (int d = 0; d < n; ++d) {
				Mean[d] = parameters.get(index)[d + 1];
			}
			for (int d = 0; d < n; ++d) {
				Sigma[d] = parameters.get(index)[n + d + 1];

			}

			Amplitudelist.add(Amplitude);

			Meanlist.add(Mean);

			Sigmalist.add(Sigma);

		}

		for (int index = 0; index < parameters.size(); ++index) {

			AddGaussian.addGaussian(imgout, Amplitudelist.get(index), Meanlist.get(index), Sigmalist.get(index));

		}

	}
	
	public static void Drawshortline(RandomAccessibleInterval<FloatType> imgout, ArrayList<Fakeline> linearray,
			double slope, double intercept, final double[] startpos,
			final double[] endpos, final double[] sigma) {

		int ndims = imgout.numDimensions();
		double [] startline = new double[ndims];
		double [] endline = new double[ndims];
		
		final double[] tmppos = new double[ndims];
		
		for (int d = 0; d < ndims; ++d){
			
			
			final double locationdiff = startpos[d] - endpos[d];
			final boolean minsearch = locationdiff > 0;
			tmppos[d] = startpos[d];
			
			startline[d] = minsearch ? endpos[d] : startpos[d];
			endline[d] = minsearch ? tmppos[d] : endpos[d];
			
			
			
		}
		
		final double stepsize = 1;
		final double[] steppos = new double[ndims];
		int count = 0;
		double distance = 0;
		while (true) {
			
			steppos[0] = startline[0] + count * stepsize / Math.sqrt(1 + slope * slope);
			steppos[1] = startline[1] + count * stepsize * slope / Math.sqrt(1 + slope * slope);
			
			AddGaussian.addGaussian(imgout, 1.0,steppos, sigma);

			distance = Distance(startline, steppos);
			
			count++;

			
			if (steppos[0] >= endline[0] || steppos[1] >= endline[1]  )
				break;
		}
		Fakeline singleline = new Fakeline(distance, slope, intercept, startline, endline);
		linearray.add(singleline);
		
	}
	
	// Draw a line between starting and end point
	public static void DrawfinalLine(RandomAccessibleInterval<FloatType> imgout,
			final double[] final_param, final double[] sigma) {

		int ndims = imgout.numDimensions();
		
		double[] startline = new double[ndims];
		double[] endline = new double[ndims];
		
		
		for (int d = 0; d < ndims; ++d){
		startline[d] = final_param[d];
		endline[d] = final_param[ndims +d];
			}
		
		double slope = ( endline[1] - startline[1] ) / (endline[0] - startline[0]);
		final double stepsize = 0.25*(sigma[0] + sigma[1]);
		final double[] steppos = new double[ndims];
		int count = 0;
		
		if (slope >= 0){
		while (true) {
			
			steppos[0] = startline[0] + count * stepsize / Math.sqrt(1 + slope * slope);
			steppos[1] = startline[1] + count * stepsize * slope / Math.sqrt(1 + slope * slope);
			
			AddGaussian.addGaussian(imgout, 1.0,steppos, sigma);

			count++;
			
			if (steppos[0] >= endline[0] || steppos[1] >= endline[1]  )
				break;
		}
		}
		int negcount = 0;
		if (slope < 0){
			while (true) {
				
				steppos[0] = startline[0] + negcount * stepsize / Math.sqrt(1 + slope * slope);
				steppos[1] = startline[1] + negcount * stepsize * slope / Math.sqrt(1 + slope * slope);
				
				AddGaussian.addGaussian(imgout, 1.0,steppos, sigma);

				negcount++;
				
				if (steppos[0] >= endline[0] || steppos[1] <= endline[1]  )
					break;
			}
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
		double sigmasq, sigma = 1;
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
			if (distance < sigma)
				intensity = (1 / (sigma * Math.sqrt(2 * Math.PI))) * Math.exp(-distance * distance / (2 * sigmasq));
			else
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

	// Compute the local maxima of your model image of the detected lines
	// all the points in that image with value > 0 are the guess centroids
	public static void MakeHTguess(RandomAccessibleInterval<FloatType> Modelimg,
			PointSampleList<FloatType> centroidlist) {

		Cursor<FloatType> modelcursor = Views.iterable(Modelimg).localizingCursor();

		while (modelcursor.hasNext()) {
			modelcursor.fwd();
			final FloatType val = new FloatType(0);
			if (modelcursor.get().compareTo(val) > 0) {

				Point centroid = new Point(modelcursor);

				centroidlist.add(centroid, modelcursor.get().copy());

			}

		}

	}

	// Compute the local maxima of your model image of the detected lines
	// all the points in that image with value > 0 are the guess centroids
	public static void Getcentroids(RandomAccessibleInterval<FloatType> Modelimg,
			PointSampleList<FloatType> centroidlist) {

		Cursor<FloatType> modelcursor = Views.iterable(Modelimg).localizingCursor();

		while (modelcursor.hasNext()) {
			modelcursor.fwd();
			final FloatType val = new FloatType(0);
			if (modelcursor.get().compareTo(val) > 0) {

				Point centroid = new Point(modelcursor);

				centroidlist.add(centroid, modelcursor.get().copy());
			}

		}

	}

	public static double MaxIntensitylabel(RandomAccessibleInterval<FloatType> inputimg,
			RandomAccessibleInterval<IntType> intimg, final int label) {
		double maxIntensity = Double.NEGATIVE_INFINITY;
		Cursor<IntType> intcursor = Views.iterable(intimg).localizingCursor();
		RandomAccess<FloatType> ranac = inputimg.randomAccess();

		while (intcursor.hasNext()) {
			intcursor.fwd();

			int i = intcursor.get().get();
			if (i == label) {

				ranac.setPosition(intcursor);

				if (ranac.get().get() > maxIntensity) {

					maxIntensity = ranac.get().get();
				}

			}

		}

		return maxIntensity;
	}

	public static void DrawTruncatedline(
			RandomAccessibleInterval<FloatType> imgout,
			RandomAccessibleInterval<FloatType> inputimg, 
			Img<IntType> intimg, 
			double slope, 
			double intercept, 
			int label) {

		int n = imgout.numDimensions();
		final double[] realpos = new double[n];
		double sigmasq, sigma = 1.0;
		sigmasq = sigma * sigma;
		final Cursor<FloatType> inputcursor = Views.iterable(imgout).localizingCursor();
		double[] newposition = new double[n];
		RandomAccess<IntType> ranac = intimg.randomAccess();
		final RandomAccess<FloatType> ranacinput = inputimg.randomAccess();
		double[] minVal = { Double.MAX_VALUE, Double.MAX_VALUE };
		double[] maxVal = { Double.MIN_VALUE, Double.MIN_VALUE };
		while (inputcursor.hasNext()) {

			inputcursor.fwd();
			inputcursor.localize(realpos);
			ranacinput.setPosition(inputcursor);
		


			ranac.setPosition(inputcursor);
			int i = ranac.get().get();

			if (i == label) {
				inputcursor.localize(newposition);
				long pointonline = (long) (newposition[1] - slope * newposition[0] - intercept);

				// To get the min and max co-rodinates along the line so we have starting points to
				// move on the line smoothly
				
				if (pointonline == 0 ) {
					
						
					for (int d = 0; d < n; ++d) {
						if (inputcursor.getDoublePosition(d) <= minVal[d]) 
							minVal[d] = inputcursor.getDoublePosition(d);
						
						if (inputcursor.getDoublePosition(d) >= maxVal[d]) 
							maxVal[d] = inputcursor.getDoublePosition(d);
						
					}
					

				}

			}
		
		
		}
		
		
		
		
		final double[] steppos = new double[n];
		int count = 0;
		double stepsize = 1;
		if (slope >= 0){
		while (true) {
			
			steppos[0] = minVal[0] + count * stepsize / Math.sqrt(1 + slope * slope);
			steppos[1] = minVal[1] + count * stepsize * slope / Math.sqrt(1 + slope * slope);
			
			AddGaussian.addGaussian(imgout, 1.0,steppos, new double[] {sigma,sigma});

			count++;
			
			if (steppos[0] >= maxVal[0] || steppos[1] >= maxVal[1]  )
				break;
		}
		}
		int negcount = 0;
		if (slope < 0){
			while (true) {
				
				steppos[0] = minVal[0] + negcount * stepsize / Math.sqrt(1 + slope * slope);
				steppos[1] = maxVal[1] + negcount * stepsize * slope / Math.sqrt(1 + slope * slope);
				
				AddGaussian.addGaussian(imgout, 1.0,steppos, new double[] {sigma,sigma});

				negcount++;
				
				if (steppos[0] >= maxVal[0] || steppos[1] <= minVal[1]  )
					break;
			}
			}
		
		
	}

	

	public static void Drawexactline(RandomAccessibleInterval<FloatType> imgout, Img<IntType> intimg, double slope,
			double intercept, int label) {

		int n = imgout.numDimensions();
		final double[] realpos = new double[n];
		double sigmasq, sigma = 1.0;
		sigmasq = sigma * sigma;
		final Cursor<FloatType> inputcursor = Views.iterable(imgout).localizingCursor();

		while (inputcursor.hasNext()) {

			inputcursor.fwd();
			inputcursor.localize(realpos);

			// To set the pixel intensity as the shortest distance to the curve
			double distance = 0;
			double intensity = 0;

			Finalfunction linefunction = new Finalfunction(realpos, slope, intercept);
			distance = linefunction.Linefunctiondist();
			final RandomAccess<FloatType> outbound = imgout.randomAccess();
			outbound.setPosition(inputcursor);

			if (distance < 5 * sigma)
				intensity = (1 / (sigma * Math.sqrt(2 * Math.PI))) * Math.exp(-distance * distance / (2 * sigmasq));
			else
				intensity = 0;

			RandomAccess<IntType> ranac = intimg.randomAccess();
			ranac.setPosition(inputcursor);
			int i = ranac.get().get();

			if (i == label) {

				outbound.get().setReal(intensity);

				final double[] position = new double[n];
				outbound.localize(position);

			}

		}

	}

	
	
	
	
	public static double Distance(final double[] cordone, final double[] cordtwo) {

		double distance = 0;

		for (int d = 0; d < cordone.length; ++d) {

			distance += Math.pow((cordone[d] - cordtwo[d]), 2);

		}
		return Math.sqrt(distance);
	}
}
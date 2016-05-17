package varun;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import com.sun.tools.javac.util.Pair;

import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Line;
import ij.gui.Overlay;
import net.imglib2.Cursor;
import net.imglib2.Point;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealPoint;
import net.imglib2.algorithm.localextrema.RefinedPeak;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import varun.GetLocalmaxmin.IntensityType;

public class OverlayLines {

	// OverlayLines for an input ArrayList<RefinedPeak<Point>>

	public static void Overlay(RandomAccessibleInterval<FloatType> inputimg,
			ArrayList<RefinedPeak<Point>> SubpixelMinlist, double[] sizes, double[] min, double[] max) {

		double[] points = new double[inputimg.numDimensions()];

		ImageStack stack = new ImageStack((int) inputimg.dimension(0), (int) inputimg.dimension(1));

		stack.addSlice(ImageJFunctions.wrap(inputimg, "").getProcessor());

		ImagePlus imp = new ImagePlus("scale space hough", stack);
		imp.show();

		Overlay o = imp.getOverlay();

		if (o == null) {
			o = new Overlay();
			imp.setOverlay(o);
		}

		o.clear();

		for (int index = 0; index < SubpixelMinlist.size(); ++index) {
			points = TransformCordinates.transformfwd(new double[] { SubpixelMinlist.get(index).getDoublePosition(0),
					SubpixelMinlist.get(index).getDoublePosition(1) }, sizes, min, max);
		//	System.out.println(" Found Peaks at :" + "Theta: " + points[0] + " Rho: " + points[1]);

			Line newline = new Line(0, points[1] / Math.sin(Math.toRadians(points[0])), inputimg.dimension(0),
					points[1] / Math.sin(Math.toRadians(points[0]))
							- inputimg.dimension(0) / Math.tan(Math.toRadians(points[0])));

			newline.setStrokeColor(Color.GREEN);
			newline.setStrokeWidth(0.8);

			o.add(newline);
		}
		imp.updateAndDraw();
	}

	public static Pair<Integer, Integer>  IndexofInt(RandomAccessibleInterval<FloatType> inputimg,
			ArrayList<RefinedPeak<Point>> SubpixelMinlist, double[] sizes,
			double[] min, double[] max) {
		RandomAccessibleInterval<FloatType> imgout = new ArrayImgFactory<FloatType>().create(inputimg,
				new FloatType());
		double[] points = new double[imgout.numDimensions()];

		int maxcount = 0;
		int maxindex = 0;
		for (int index = 0; index < SubpixelMinlist.size(); ++index) {
			points = TransformCordinates.transformfwd(new double[] { SubpixelMinlist.get(index).getDoublePosition(0),
					SubpixelMinlist.get(index).getDoublePosition(1) }, sizes, min, max);

			double slope = -1.0 / Math.tan(Math.toRadians(points[0]));
			double intercept = points[1] / Math.sin(Math.toRadians(points[0]));
		//	System.out.println(" Found Peaks at :" + "Theta: " + points[0] + " Rho: " + points[1]);
			PushCurves.Drawexactline(imgout, slope, intercept, IntensityType.Gaussian);

			RandomAccess<FloatType> inran = inputimg.randomAccess();
			Cursor<FloatType> outcursor = Views.iterable(imgout).localizingCursor();

			int count = 0;
			while (outcursor.hasNext()) {
				outcursor.fwd();

				
				if (outcursor.get().get() > 0) {
					inran.setPosition(outcursor);

					if (inran.get().get() > 0)
						count++;

					if (count > maxcount) {
						maxcount = count;
						maxindex = index;

					}

				}

			}
		}
		Pair<Integer, Integer> pair = new Pair<Integer, Integer>(maxcount, maxindex);
		return pair;
	}

	public static ArrayList<RefinedPeak<Point>>  ReducedList(RandomAccessibleInterval<FloatType> inputimg,
			ArrayList<RefinedPeak<Point>> SubpixelMinlist, double[] sizes,
			double[] min, double[] max) {
		RandomAccessibleInterval<FloatType> imgout = new ArrayImgFactory<FloatType>().create(inputimg,
				new FloatType());
		double[] points = new double[imgout.numDimensions()];

		int maxcount = 0;
		int maxindex = 0;
		for (int index = 0; index < SubpixelMinlist.size(); ++index) {
			points = TransformCordinates.transformfwd(new double[] { SubpixelMinlist.get(index).getDoublePosition(0),
					SubpixelMinlist.get(index).getDoublePosition(1) }, sizes, min, max);

			double slope = -1.0 / Math.tan(Math.toRadians(points[0]));
			double intercept = points[1] / Math.sin(Math.toRadians(points[0]));
		//	System.out.println(" Found Peaks at :" + "Theta: " + points[0] + " Rho: " + points[1]);
			PushCurves.Drawexactline(imgout, slope, intercept, IntensityType.Gaussian);

			RandomAccess<FloatType> inran = inputimg.randomAccess();
			Cursor<FloatType> outcursor = Views.iterable(imgout).localizingCursor();

			int count = 0;
			while (outcursor.hasNext()) {
				outcursor.fwd();

				
				if (outcursor.get().get() > 0) {
					inran.setPosition(outcursor);

					if (inran.get().get() > 0)
						count++;

					if (count > maxcount) {
						maxcount = count;
						maxindex = index;

					}

				}

			}
		}
	//	System.out.println("Main file: "+ maxcount + " " + maxindex);
		ArrayList<RefinedPeak<Point>> MainMinlist = new ArrayList<RefinedPeak<Point>>(inputimg.numDimensions());
		if (maxcount > 0) {
		MainMinlist.add(SubpixelMinlist.get(maxindex));
		}
		return MainMinlist;
	}
	
	public static void OverlayExactline(RandomAccessibleInterval<FloatType> inputimg,
			ArrayList<RefinedPeak<Point>> SubpixelMinlist, double[] sizes,
			double[] min, double[] max) {
		RandomAccessibleInterval<FloatType> imgout = new ArrayImgFactory<FloatType>().create(inputimg,
				new FloatType());
		double[] points = new double[imgout.numDimensions()];

		int maxcount = 0;
		int maxindex = 0;
		for (int index = 0; index < SubpixelMinlist.size(); ++index) {
			points = TransformCordinates.transformfwd(new double[] { SubpixelMinlist.get(index).getDoublePosition(0),
					SubpixelMinlist.get(index).getDoublePosition(1) }, sizes, min, max);

			double slope = -1.0 / Math.tan(Math.toRadians(points[0]));
			double intercept = points[1] / Math.sin(Math.toRadians(points[0]));
			
			PushCurves.Drawexactline(imgout, slope, intercept, IntensityType.Gaussian);

			RandomAccess<FloatType> inran = inputimg.randomAccess();
			Cursor<FloatType> outcursor = Views.iterable(imgout).localizingCursor();

			int count = 0;
			while (outcursor.hasNext()) {
				outcursor.fwd();

				if (outcursor.get().get() > 0) {
					inran.setPosition(outcursor);

					if (inran.get().get() > 0)
						count++;

					if (count > maxcount) {
						maxcount = count;
						maxindex = index;

					}

				}

			}

		}

	//	System.out.println(maxcount + " " + maxindex);
		if (maxcount > 0) {

			ImageStack stack = new ImageStack((int) inputimg.dimension(0), (int) inputimg.dimension(1));

			stack.addSlice(ImageJFunctions.wrap(inputimg, "").getProcessor());

			ImagePlus imp = new ImagePlus("scale space hough", stack);
			imp.show();

			Overlay o = imp.getOverlay();

			if (o == null) {
				o = new Overlay();
				imp.setOverlay(o);
			}

			o.clear();
			points = TransformCordinates.transformfwd(new double[] { SubpixelMinlist.get(maxindex).getDoublePosition(0),
					SubpixelMinlist.get(maxindex).getDoublePosition(1) }, sizes, min, max);
		//	System.out.println("All: "+ " Found Peaks at :" + "Theta: " + points[0] + " Rho: " + points[1]);
			Line newline = new Line(0, points[1] / Math.sin(Math.toRadians(points[0])), inputimg.dimension(0),
					points[1] / Math.sin(Math.toRadians(points[0]))
							- inputimg.dimension(0) / Math.tan(Math.toRadians(points[0])));

			newline.setStrokeColor(Color.GREEN);
			newline.setStrokeWidth(0.2);

			o.add(newline);

			imp.updateAndDraw();
		}

	}

	
	

}

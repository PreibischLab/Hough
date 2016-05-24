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
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealCursor;
import net.imglib2.RealPoint;
import net.imglib2.RealPointSampleList;
import net.imglib2.algorithm.localextrema.RefinedPeak;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.labeling.NativeImgLabeling;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import varun.GetLocalmaxmin.IntensityType;
import varun.PerformWatershedding.Lineobjects;

public class OverlayLines {
	public static final class Simulatedline {
		final int Label;
		final double[] point;
		final FloatType Value;

		protected Simulatedline(final int Label, final double[] point, final FloatType Value) {
			this.Label = Label;
			this.point = point;
			this.Value = Value;
			

		}
	}

	public static void OverlayObject(RandomAccessibleInterval<FloatType> inputimg, ArrayList<Lineobjects> linelist) {

		double rho;
		double theta;
		long minX, maxX, minY, maxY;

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
		for (int index = 0; index < linelist.size(); ++index) {
			rho = linelist.get(index).Rho;
			theta = linelist.get(index).Theta;
			minX = linelist.get(index).boxXmin;
			maxX = linelist.get(index).boxXmax;

			Line newline = new Line(minX,
					rho / Math.sin(Math.toRadians(theta)) - minX / Math.tan(Math.toRadians(theta)), maxX,
					rho / Math.sin(Math.toRadians(theta)) - maxX / Math.tan(Math.toRadians(theta)));

			newline.setStrokeColor(Color.GREEN);
			newline.setStrokeWidth(0.8);

			o.add(newline);
		}
		imp.updateAndDraw();
	}

	// OverlayLines for an input ArrayList<RefinedPeak<Point>>

	public static void Overlay(RandomAccessibleInterval<FloatType> inputimg,
			ArrayList<RefinedPeak<Point>> SubpixelMinlist, double[] sizes, double[] min, double[] max, int label,
			long minX, long maxX) {

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
			// System.out.println(" Found Peaks at :" + "Theta: " + points[0] +
			// " Rho: " + points[1]);

			Line newline = new Line(minX,
					points[1] / Math.sin(Math.toRadians(points[0])) - minX / Math.tan(Math.toRadians(points[0])), maxX,
					points[1] / Math.sin(Math.toRadians(points[0])) - (maxX) / Math.tan(Math.toRadians(points[0])));

			newline.setStrokeColor(Color.GREEN);
			newline.setStrokeWidth(0.8);

			o.add(newline);
		}
		imp.updateAndDraw();
	}

	public static Pair<Integer, Integer> IndexofInt(RandomAccessibleInterval<FloatType> inputimg,
			ArrayList<RefinedPeak<Point>> SubpixelMinlist, double[] sizes, double[] min, double[] max) {
		RandomAccessibleInterval<FloatType> imgout = new ArrayImgFactory<FloatType>().create(inputimg, new FloatType());
		double[] points = new double[imgout.numDimensions()];

		int maxcount = 0;
		int maxindex = 0;
		for (int index = 0; index < SubpixelMinlist.size(); ++index) {
			points = TransformCordinates.transformfwd(new double[] { SubpixelMinlist.get(index).getDoublePosition(0),
					SubpixelMinlist.get(index).getDoublePosition(1) }, sizes, min, max);

			double slope = -1.0 / Math.tan(Math.toRadians(points[0]));
			double intercept = points[1] / Math.sin(Math.toRadians(points[0]));
			// System.out.println(" Found Peaks at :" + "Theta: " + points[0] +
			// " Rho: " + points[1]);
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

	public static ArrayList<RefinedPeak<Point>> ReducedList(RandomAccessibleInterval<FloatType> inputimg,
			ArrayList<RefinedPeak<Point>> SubpixelMinlist, double[] sizes, double[] min, double[] max,
			double minlength) {

		RandomAccessibleInterval<FloatType> imgout = new ArrayImgFactory<FloatType>().create(inputimg, new FloatType());
		double[] points = new double[imgout.numDimensions()];

		int maxcount = 0;
		int maxindex = 0;
		for (int index = 0; index < SubpixelMinlist.size(); ++index) {
			points = TransformCordinates.transformfwd(new double[] { SubpixelMinlist.get(index).getDoublePosition(0),
					SubpixelMinlist.get(index).getDoublePosition(1) }, sizes, min, max);

			double slope = -1.0 / Math.tan(Math.toRadians(points[0]));
			double intercept = points[1] / Math.sin(Math.toRadians(points[0]));
			// System.out.println(" Found Peaks at :" + "Theta: " + points[0] +
			// " Rho: " + points[1]);
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
		// System.out.println("Main file: "+ maxcount + " " + maxindex);
		ArrayList<RefinedPeak<Point>> MainMinlist = new ArrayList<RefinedPeak<Point>>(inputimg.numDimensions());
		if (maxcount > minlength) {
			MainMinlist.add(SubpixelMinlist.get(maxindex));
		}
		return MainMinlist;
	}

	public static double[] GetRhoTheta(ArrayList<RefinedPeak<Point>> MainMinlist, double[] sizes, double[] min,
			double[] max, double minlength) {

		double[] points = new double[sizes.length];
		for (int index = 0; index < MainMinlist.size(); ++index) {

			points = TransformCordinates.transformfwd(new double[] { MainMinlist.get(index).getDoublePosition(0),
					MainMinlist.get(index).getDoublePosition(1) }, sizes, min, max);
		}
		return points;
	}

	public static void GetAlllines(
			RandomAccessibleInterval<FloatType> imgout,
			ArrayList<Simulatedline> totalsimline,
			Img<IntType> intimg, 
			ArrayList<Lineobjects> linelist,
			double cutoff) {

		for (int index = 0; index < linelist.size(); ++index) {

			final int label = linelist.get(index).Label;
			final double rho = linelist.get(index).Rho;
			final double theta = linelist.get(index).Theta;
			
			ArrayList<Simulatedline> simline = new ArrayList<Simulatedline>();
			double slope = -1.0 / Math.tan(Math.toRadians(theta));
			double intercept = rho / Math.sin(Math.toRadians(theta));

			PushCurves.Drawexactline(imgout,simline,intimg, slope, intercept, label, cutoff);
			for (int simindex = 0; simindex< simline.size(); ++simindex)
			totalsimline.add(simline.get(simindex));
		}
	}

}

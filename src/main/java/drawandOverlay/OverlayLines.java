package drawandOverlay;

import java.awt.Color;
import java.util.ArrayList;
import com.sun.tools.javac.util.Pair;

import houghandWatershed.TransformCordinates;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Line;
import ij.gui.Overlay;
import labeledObjects.Lineobjects;
import labeledObjects.Simpleobject;
import net.imglib2.Cursor;
import net.imglib2.Point;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealPointSampleList;
import net.imglib2.algorithm.localextrema.RefinedPeak;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import preProcessing.GetLocalmaxmin;
import preProcessing.GetLocalmaxmin.IntensityType;

public class OverlayLines {
	

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
			minX = linelist.get(index).boxmin[0];
			maxX = linelist.get(index).boxmax[0];

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

	

	public static ArrayList<RefinedPeak<Point>> ReducedList(RandomAccessibleInterval<FloatType> inputimg,
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
		// System.out.println("Main file: "+ maxcount + " " + maxindex);
		ArrayList<RefinedPeak<Point>> MainMinlist = new ArrayList<RefinedPeak<Point>>(inputimg.numDimensions());
		if (maxcount > 0) {
			MainMinlist.add(SubpixelMinlist.get(maxindex));
		}
		return MainMinlist;
	}

	public static double[] GetRhoTheta(ArrayList<RefinedPeak<Point>> MainMinlist, double[] sizes, double[] min,
			double[] max) {

		double[] points = new double[sizes.length];
		for (int index = 0; index < MainMinlist.size(); ++index) {

			points = TransformCordinates.transformfwd(new double[] { MainMinlist.get(index).getDoublePosition(0),
					MainMinlist.get(index).getDoublePosition(1) }, sizes, min, max);
		}
		return points;
	}

	public static void GetAlllines(
			RandomAccessibleInterval<FloatType> imgout,
			RandomAccessibleInterval<FloatType> inputimg,
			Img<IntType> intimg, 
			ArrayList<Lineobjects> linelist,
			ArrayList<Simpleobject> lineobject,
			final long radius) {

		for (int index = 0; index < linelist.size(); ++index) {

			final int label = linelist.get(index).Label;
			final double rho = linelist.get(index).Rho;
			final double theta = linelist.get(index).Theta;
			
			
			double slope = -1.0 / (Math.tan(Math.toRadians(theta)));
			double intercept = rho / Math.sin(Math.toRadians(theta));
			if (Math.abs(slope)!=Double.POSITIVE_INFINITY || Math.abs(intercept)!=Double.POSITIVE_INFINITY  ){
			final Simpleobject simpleobj = new Simpleobject(label, slope, intercept);
			lineobject.add(simpleobj);
			}
		//	System.out.println(slope +"  "+ theta);
			//PushCurves.Drawexactline(testimgout,intimg, slope, intercept, label);
			
			if (Math.abs(slope)!=Double.POSITIVE_INFINITY )
			PushCurves.DrawTruncatedline(imgout, inputimg, intimg, slope, intercept, label);
			
			
			
			
		}
		
	}
	public static void GetCurrentlines(
			RandomAccessibleInterval<FloatType> imgout,
			RandomAccessibleInterval<FloatType> maximgout,
			Img<IntType> intimg, 
			ArrayList<Lineobjects> linelist,
			int currentlabel) {

		for (int index = 0; index < linelist.size(); ++index) {

			final int label = linelist.get(index).Label;
			final double rho = linelist.get(index).Rho;
			final double theta = linelist.get(index).Theta;
		
		if (label == currentlabel){
			
			double slope = -1.0 / Math.tan(Math.toRadians(theta));
			double intercept = rho / Math.sin(Math.toRadians(theta));

			PushCurves.Drawexactline(imgout,intimg, slope,intercept, label);
			maximgout = GetLocalmaxmin.FindandDisplayLocalMaxima(imgout,
					IntensityType.Original, new double[]{1,1});
			
		}
		}
	}

}
package varun;

import java.awt.Color;
import java.util.ArrayList;

import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Line;
import ij.gui.Overlay;
import net.imglib2.Point;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealPoint;
import net.imglib2.algorithm.localextrema.RefinedPeak;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;

public class OverlayLines {
	
	
	// OverlayLines for an input ArrayList<RefinedPeak<Point>>

		public static void Overlay(RandomAccessibleInterval<FloatType> inputimg,
				ArrayList<RefinedPeak<Point>> SubpixelMinlist, double[] sizes, double[] min, double[] max) {

			double[] points = new double[inputimg.numDimensions()];

			ImageStack stack = new ImageStack((int) inputimg.dimension(0), (int) inputimg.dimension(1));

			stack.addSlice(ImageJFunctions.wrap(inputimg, "").getProcessor());
			new ImageJ();

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
				//System.out.println(" Found Peaks at :" + "Theta: " + points[0] + " Rho: " + points[1]);

				Line newline = new Line(0, points[1] / Math.sin(Math.toRadians(points[0])), inputimg.dimension(0),
						points[1] / Math.sin(Math.toRadians(points[0]))
								- inputimg.dimension(0) / Math.tan(Math.toRadians(points[0])));

				newline.setStrokeColor(Color.RED);
				newline.setStrokeWidth(0.8);

				o.add(newline);
			}
			imp.updateAndDraw();
		}

		
	
	
	
	
	// OverlayLines for an input ArrayList<RealPoint>

		public static void Overlaysecond(RandomAccessibleInterval<FloatType> inputimg, ArrayList<RealPoint> Minlist,
				double[] sizes, double[] min, double[] max) {

			double[] points = new double[inputimg.numDimensions()];

			ImageStack stack = new ImageStack((int) inputimg.dimension(0), (int) inputimg.dimension(1));

			stack.addSlice(ImageJFunctions.wrap(inputimg, "").getProcessor());
			new ImageJ();

			ImagePlus imp = new ImagePlus("scale space hough", stack);
			imp.show();

			Overlay o = imp.getOverlay();

			if (o == null) {
				o = new Overlay();
				imp.setOverlay(o);
			}

			o.clear();

			for (int index = 0; index < Minlist.size(); ++index) {
				points = TransformCordinates.transformfwd(
						new double[] { Minlist.get(index).getDoublePosition(0), Minlist.get(index).getDoublePosition(1) },
						sizes, min, max);
				System.out.println(" Found Peaks at :" + "Theta: " + points[0] + " Rho: " + points[1]);

				Line newline = new Line(0, points[1] / Math.sin(Math.toRadians(points[0])), inputimg.dimension(0),
						points[1] / Math.sin(Math.toRadians(points[0]))
								- inputimg.dimension(0) / Math.tan(Math.toRadians(points[0])));

				newline.setStrokeColor(Color.RED);
				newline.setStrokeWidth(0.8);

				o.add(newline);
			}
			imp.updateAndDraw();
		}

}

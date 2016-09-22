package velocityanalyser;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import com.sun.tools.javac.util.Pair;

import drawandOverlay.OverlayLines;
import drawandOverlay.PushCurves;
import houghandWatershed.PerformWatershedding;
import ij.ImageJ;
import ij.ImagePlus;
import labeledObjects.Lineobjects;
import labeledObjects.Simpleobject;
import net.imglib2.Point;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.stats.Normalize;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.IntervalView;
import net.imglib2.view.Views;
import peakFitter.Linefitter;
import preProcessing.Kernels;

public class Velocitydetector {

	public static void main(String[] args) throws Exception {
		new ImageJ();
		// Load the stack of images
		final RandomAccessibleInterval<FloatType> img = util.ImgLib2Util
				.openAs32Bit(new File("src/main/resources/Fake_moving.tif"), new ArrayImgFactory<FloatType>());
		int ndims = img.numDimensions();
		new Normalize();

		FloatType minval = new FloatType(0);
		FloatType maxval = new FloatType(1);
		Normalize.normalize(Views.iterable(img), minval, maxval);

		// Declare all the constants needed by the program here:

		final double[] psf = { 1.7, 1.8 };
		final long radius = (long) Math.ceil(Math.sqrt(psf[0] * psf[0] + psf[1] * psf[1]));
		final int minlength = 5;
		double[] final_param = new double[2 * ndims + 2];

		// Show the stack
		ImagePlus imp = ImageJFunctions.show(img);

		// Do Hough transform on the Fisrt seed image

		IntervalView<FloatType> groundframe = Views.hyperSlice(img, ndims - 1, 0);
		RandomAccessibleInterval<FloatType> preprocessedimg = new ArrayImgFactory<FloatType>().create(groundframe,
				new FloatType());
		preprocessedimg = Kernels.Meanfilterandsupress(groundframe, radius);

		System.out.println("Doing Hough transform in labels: ");

		PerformWatershedding Houghobject = new PerformWatershedding(groundframe, preprocessedimg, minlength);

		Pair<Img<IntType>, ArrayList<Lineobjects>> linepair = Houghobject.DowatersheddingandHough();

		// Display the Hough detection
		RandomAccessibleInterval<FloatType> imgout = new ArrayImgFactory<FloatType>().create(groundframe,
				new FloatType());
		final ArrayList<Simpleobject> simpleobject = new ArrayList<Simpleobject>();
		OverlayLines.GetAlllines(imgout, groundframe, linepair.fst, linepair.snd, simpleobject, radius);
		ImageJFunctions.show(imgout).setTitle("Rough-Reconstruction");

		Linefitter MTline = new Linefitter(groundframe, linepair.fst);

		RandomAccessibleInterval<FloatType> gaussimg = new ArrayImgFactory<FloatType>().create(groundframe,
				new FloatType());

		double distance = 0;
		ArrayList<double[]> final_paramlist = new ArrayList<double[]>();
		for (int index = 0; index < simpleobject.size(); ++index) {

			// Do gradient descent to improve the Hough detected lines
			final_param = MTline.Getfinallineparam(simpleobject.get(index).Label, simpleobject.get(index).slope,
					simpleobject.get(index).intercept, psf, minlength, false);

			if (final_param != null) {
				final_paramlist.add(final_param);
				final double[] cordone = { final_param[0], final_param[1] };
				final double[] cordtwo = { final_param[2], final_param[3] };

				distance = MTline.Distance(cordone, cordtwo);

				PushCurves.DrawfinalLine(gaussimg, final_param, psf);

				System.out.println("Fits :" + "StartX:" + final_param[0] + " StartY:" + final_param[1] + " " + "EndX:"
						+ final_param[2] + "EndY: " + final_param[3]);
				System.out.println("Length: " + distance);

				
			}
		}

		for (int i = 1; i < img.dimension(ndims - 1); ++i) {

			
			IntervalView<FloatType> currentframe = Views.hyperSlice(img, ndims - 1, i);
			
			PerformWatershedding Watershedobject = new PerformWatershedding(currentframe, minlength);
			
			RandomAccessibleInterval<IntType> currentlabelledimg = Watershedobject.Dowatersheddingonly();
			
			Linefitter currentline = new Linefitter(currentframe, currentlabelledimg);
			
			
			for (int index = 0; index < final_paramlist.size(); ++index) {

				Point linepoint = new Point(ndims - 1);
				linepoint.setPosition(
						new long[] { (long) final_paramlist.get(index)[0], (long) final_paramlist.get(index)[1] });
				
				
				
				 int currentlabel = currentline.Getlabel(linepoint);
				
				 double[] final_trackparam =
						currentline.Getfinaltrackparam(final_paramlist.get(index),
								currentlabel, psf, i, false);
				System.out.println("Frame:" + i + " " +  "Fits :" + currentlabel + " "+ "StartX:" + final_trackparam[0] 
						+ " StartY:" + final_trackparam[1] + " " + "EndX:"
						+ final_trackparam[2] + "EndY: " + final_trackparam[3]);
				
				
				
				
			}

		}

	}

}

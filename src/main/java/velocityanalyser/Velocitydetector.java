package velocityanalyser;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.SimpleWeightedGraph;

import com.sun.tools.javac.util.Pair;

import drawandOverlay.DisplayGraph;
import drawandOverlay.OverlayLines;
import drawandOverlay.PushCurves;
import graphconstructs.Staticproperties;
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
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.IntervalView;
import net.imglib2.view.Views;
import peakFitter.Linefitter;
import preProcessing.GetLocalmaxmin;
import preProcessing.GlobalThresholding;
import preProcessing.Kernels;
import preProcessing.MedianFilter2D;

public class Velocitydetector {

	public static void main(String[] args) throws Exception {

		/***
		 * Hough Transform to detect Microtubules and track the growth at
		 * Sub-pixel accuracy. Optimizers used: Levenberg-Marqurat solver and
		 * Weighted centre of mass fits. Program reqires PSF of the microscope
		 * to be computed and analysed and takes the determined Sigmas as the
		 * input. @ Varun Kapoor
		 */

		new ImageJ();

		// Load the stack of images
		final RandomAccessibleInterval<FloatType> img = util.ImgLib2Util
				.openAs32Bit(new File("../res/2016-05-26-test-brighter.tif"), new ArrayImgFactory<FloatType>());
		int ndims = img.numDimensions();

		// Normalize the intensity of the whole stack to be between min and max
		// value specified here
		new Normalize();

		FloatType minval = new FloatType(0);
		FloatType maxval = new FloatType(1);
		Normalize.normalize(Views.iterable(img), minval, maxval);

		// Declare all the constants needed by the program here:

		final double[] psf = { 1.65, 1.47 };
		final long radius = (long) Math.ceil(Math.sqrt(psf[0] * psf[0] + psf[1] * psf[1]));
		final int minlength = 5;
		double[] final_param = new double[2 * ndims + 3];

		// Show the stack
		ImagePlus imp = ImageJFunctions.show(img);
		ArrayList<ArrayList<Staticproperties>> Allstartandend = new ArrayList<ArrayList<Staticproperties>>();

		// Do Hough transform on the First seed image

		IntervalView<FloatType> groundframe = Views.hyperSlice(img, ndims - 1, 0);

		RandomAccessibleInterval<FloatType> preprocessedimg = new ArrayImgFactory<FloatType>().create(groundframe,
				new FloatType());

		preprocessedimg = Kernels.Supressthresh(groundframe);
		// Kernels.Meanfilterandsupress(groundframe, radius);

		ImageJFunctions.show(preprocessedimg);

		
		final Float ThresholdValue = GlobalThresholding.AutomaticThresholding(preprocessedimg);
		RandomAccessibleInterval<BitType> bitimg = new ArrayImgFactory<BitType>().create(preprocessedimg, new BitType());
		GetLocalmaxmin.ThresholdingBit(preprocessedimg, bitimg, ThresholdValue);
		System.out.println("Doing Hough transform in labels: ");

		PerformWatershedding Houghobject = new PerformWatershedding(preprocessedimg, bitimg,  minlength);

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
		// Run LM solver optimizer and mask fits to improve the Hough detected
		// lines
		for (int index = 0; index < simpleobject.size(); ++index) {

			final_param = MTline.Getfinallineparam(simpleobject.get(index).Label, simpleobject.get(index).slope,
					simpleobject.get(index).intercept, psf, minlength, true);

			if (final_param != null) {
				final_paramlist.add(final_param);
				final double[] cordone = { final_param[0], final_param[1] };
				final double[] cordtwo = { final_param[2], final_param[3] };

				distance = MTline.Distance(cordone, cordtwo);

				System.out.println("Fits :" + "StartX:" + final_param[0] + " StartY:" + final_param[1] + " " + "EndX:"
						+ final_param[2] + "EndY: " + final_param[3]);
				System.out.println("Length: " + distance);

			}
		}

		// Draw detected lines from the seed image
		PushCurves.DrawallLine(gaussimg, final_paramlist, psf);
		ImageJFunctions.show(gaussimg);
		// Write down the line parameters for the seed image

		for (int listindex = 0; listindex < final_paramlist.size(); ++listindex) {
			final double[] cordone = { final_paramlist.get(listindex)[0], final_paramlist.get(listindex)[1] };
			final double[] cordtwo = { final_paramlist.get(listindex)[2], final_paramlist.get(listindex)[3] };

			distance = MTline.Distance(cordone, cordtwo);

			try {
				FileWriter writer = new FileWriter("../res/2016-05-26-test-brighter.txt", true);
				writer.write("Frame:" + " " + 0 + "StartX:" + final_paramlist.get(listindex)[0] + " StartY:"
						+ final_paramlist.get(listindex)[1] + " " + "EndX:" + final_paramlist.get(listindex)[2]
						+ "EndY: " + final_paramlist.get(listindex)[3] + " "

						+ "Length: " + " " + distance);
				writer.write("\r\n");
				writer.close();

			} catch (IOException e) {
				e.printStackTrace();
			}

		}

		ArrayList<double[]> PrevFrameparam = final_paramlist;

		// Now start tracking the moving ends of the Microtubule and make
		// seperate graph for both ends

		for (int i = 1; i < img.dimension(ndims - 1); ++i) {

			IntervalView<FloatType> currentframe = Views.hyperSlice(img, ndims - 1, i);
			
			RandomAccessibleInterval<FloatType> precurrent = new ArrayImgFactory<FloatType>().create(currentframe,
					new FloatType());

			precurrent = Kernels.Supressthresh(currentframe);
			// Kernels.Meanfilterandsupress(currentframeframe, radius);

			ImageJFunctions.show(preprocessedimg);

			
			final Float currThresholdValue = GlobalThresholding.AutomaticThresholding(precurrent);
			RandomAccessibleInterval<BitType> currbitimg = new ArrayImgFactory<BitType>().create(precurrent, new BitType());
			GetLocalmaxmin.ThresholdingBit(precurrent, currbitimg, currThresholdValue);
			final Trackgrowth growthtracker = new Trackgrowth(currentframe, currbitimg, minlength, PrevFrameparam, i, psf, true);

			Pair<ArrayList<double[]>, ArrayList<Staticproperties>> pair = growthtracker.Updatetrackpoints();

			// Update the list of line parameters with the current frame
			// detection
			PrevFrameparam = pair.fst;

			// Append the object static properties with the current frame
			// detection
			Allstartandend.add(pair.snd);

			// Draw the lines detected in the current frame
			RandomAccessibleInterval<FloatType> newgaussimg = new ArrayImgFactory<FloatType>().create(groundframe,
					new FloatType());
			PushCurves.DrawallLine(newgaussimg, pair.fst, psf);
			ImageJFunctions.show(newgaussimg);
			// Write down the line parameters for the current frame image
			for (int listindex = 0; listindex < PrevFrameparam.size(); ++listindex) {
				final double[] cordone = { PrevFrameparam.get(listindex)[0], PrevFrameparam.get(listindex)[1] };
				final double[] cordtwo = { PrevFrameparam.get(listindex)[2], PrevFrameparam.get(listindex)[3] };

				distance = MTline.Distance(cordone, cordtwo);

				try {
					FileWriter writer = new FileWriter("../res/2016-05-26-test-brighter.txt", true);
					writer.write("Frame: " + i + " " + "StartX:" + PrevFrameparam.get(listindex)[0] + " StartY:"
							+ PrevFrameparam.get(listindex)[1] + " " + "EndX:" + PrevFrameparam.get(listindex)[2]
							+ "EndY: " + PrevFrameparam.get(listindex)[3] + " "

							+ "Length: " + " " + distance);
					writer.write("\r\n");
					writer.close();

				} catch (IOException e) {
					e.printStackTrace();
				}

			}

		}

		// Overlay the graphs on the stack

		// Make graph to track the start and the end point

		final int maxframe = (int) img.dimension(ndims - 1);
		final Trackstart trackerstart = new Trackstart(Allstartandend, maxframe);
		final Trackend trackerend = new Trackend(Allstartandend, maxframe);
		trackerstart.process();
		SimpleWeightedGraph<double[], DefaultWeightedEdge> graphstart = trackerstart.getResult();

		ImagePlus impstart = imp.duplicate();

		DisplayGraph displaytracksstart = new DisplayGraph(impstart, graphstart);
		displaytracksstart.getImp();
		impstart.draw();

		trackerend.process();
		SimpleWeightedGraph<double[], DefaultWeightedEdge> graphend = trackerend.getResult();

		ImagePlus impend = imp.duplicate();

		DisplayGraph displaytracksend = new DisplayGraph(impend, graphend);
		displaytracksend.getImp();
		impend.draw();

	}

}

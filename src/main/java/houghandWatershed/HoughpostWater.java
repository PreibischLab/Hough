
package houghandWatershed;

import java.io.File;
import java.util.ArrayList;

import com.sun.tools.javac.util.Pair;

import drawandOverlay.OverlayLines;
import drawandOverlay.PushCurves;
import ij.ImageJ;
import labeledObjects.Lineobjects;
import labeledObjects.Simpleobject;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.stats.Normalize;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import peakFitter.SubpixelLength;
import preProcessing.GetLocalmaxmin;
import preProcessing.GlobalThresholding;
import preProcessing.Kernels;
import preProcessing.MedianFilter2D;
import util.ImgLib2Util;

public class HoughpostWater {

	public static void main(String[] args) throws Exception {

		RandomAccessibleInterval<FloatType> biginputimg = ImgLib2Util.openAs32Bit(new File("../res/Pnonoise1.tif"));

		RandomAccessibleInterval<FloatType> processedimg = ImgLib2Util.openAs32Bit(new File("../res/Pnonoise1.tif"));
		
		new ImageJ();

		new Normalize();
		FloatType minval = new FloatType(0);
		FloatType maxval = new FloatType(1);

		Normalize.normalize(Views.iterable(biginputimg), minval, maxval);

		ImageJFunctions.show(biginputimg);
		// Define the psf of the microscope
		double[] psf = { 1.4, 1.5 };
		final long radius = (long) Math.ceil(Math.sqrt(psf[0] * psf[0] + psf[1] * psf[1]));
		// Declare minimum length of the line(in pixels) to be detected
		double minlength = 2;
		// Initialize empty images to be used later

		RandomAccessibleInterval<FloatType> inputimg = new ArrayImgFactory<FloatType>().create(biginputimg,
				new FloatType());
		RandomAccessibleInterval<FloatType> preinputimg = new ArrayImgFactory<FloatType>().create(biginputimg,
				new FloatType());
		RandomAccessibleInterval<FloatType> imgout = new ArrayImgFactory<FloatType>().create(biginputimg,
				new FloatType());
		RandomAccessibleInterval<FloatType> gaussimg = new ArrayImgFactory<FloatType>().create(biginputimg,
				new FloatType());

		
		// Preprocess image using Median Filter and suppress background
		final MedianFilter2D<FloatType> medfilter = new MedianFilter2D<FloatType>( processedimg, (int)radius );
		medfilter.process();
	    inputimg = medfilter.getResult();
		Normalize.normalize(Views.iterable(inputimg), minval, maxval);
		

		ImageJFunctions.show(inputimg).setTitle("Preprocessed image");

		// Create the Bit image for distance transform and seed image for watershedding
		final Float ThresholdValue = GlobalThresholding.AutomaticThresholding(inputimg);
		RandomAccessibleInterval<BitType> bitimg = new ArrayImgFactory<BitType>().create(inputimg, new BitType());
		GetLocalmaxmin.ThresholdingBit(inputimg, bitimg, ThresholdValue);

		// Do watershedding and Hough
		System.out.println("Doing Hough transform in labels: ");

		HoughTransform2D Houghobject = new HoughTransform2D(inputimg, bitimg, minlength);
		
		Houghobject.checkInput();
		Houghobject.process();
		Pair<RandomAccessibleInterval<IntType>, ArrayList<Lineobjects>> linepair  = Houghobject.getResult();


		// Overlay detected lines on the image
		final ArrayList<Simpleobject> simpleobject = new ArrayList<Simpleobject>();
		OverlayLines.GetAlllines(imgout, biginputimg, linepair.fst, linepair.snd, simpleobject, radius);

		ImageJFunctions.show(imgout).setTitle("Rough-Reconstruction");

		// Input the image on which you want to do the fitting (original noisy image, along with the
		// labelled image (watershedded image) 

		SubpixelLength MTline = new SubpixelLength(biginputimg, linepair.fst, simpleobject, psf, minlength);
		MTline.checkInput();
		MTline.process();
		ArrayList<double[]> final_paramlist = MTline.getResult();
		
		// Draw the detected lines
		PushCurves.DrawallLine(gaussimg, final_paramlist, psf);
		ImageJFunctions.show(gaussimg).setTitle("Exact-line");

	}
}

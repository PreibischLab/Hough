package houghandWatershed;

import java.io.File;
import java.util.ArrayList;

import com.sun.tools.javac.util.Pair;

import Spindles.Boundingbox.Objectproperties;
import drawandOverlay.OverlayLines;
import drawandOverlay.PushCurves;
import ij.ImageJ;
import labeledObjects.Lineobjects;
import net.imglib2.Cursor;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealCursor;
import net.imglib2.RealPointSampleList;
import net.imglib2.algorithm.stats.Normalize;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import peakFitter.LengthDetection;
import preProcessing.Kernels;
import preProcessing.Kernels.ProcessingType;
import util.ImgLib2Util;

public class HoughpostWater {

	public static void main(String[] args) throws Exception {

		RandomAccessibleInterval<FloatType> biginputimg = ImgLib2Util
				.openAs32Bit(new File("src/main/resources/small_mt.tif"));
		// small_mt.tif image to be used for testing
		// 2015-01-14_Seeds-1.tiff for actual
		// mt_experiment.tif for big testing
		new ImageJ();

		new Normalize();
		FloatType minval = new FloatType(0);
		FloatType maxval = new FloatType(1);
		Normalize.normalize(Views.iterable(biginputimg), minval, maxval);
		ImageJFunctions.show(biginputimg);
		final int n = biginputimg.numDimensions();
		double[] sigma = new double[biginputimg.numDimensions()];

		for (int d = 0; d < sigma.length; ++d)
			sigma[d] = 2;

		// Initialize empty images to be used later
		RandomAccessibleInterval<FloatType> inputimg = new ArrayImgFactory<FloatType>().create(biginputimg,
				new FloatType());
		
		RandomAccessibleInterval<FloatType> imgout = new ArrayImgFactory<FloatType>().create(biginputimg,
				new FloatType());
		RandomAccessibleInterval<FloatType> gaussimg = new ArrayImgFactory<FloatType>().create(biginputimg,
				new FloatType());

		// Preprocess image
		
		
		
		
		inputimg = Kernels.Preprocess(biginputimg, ProcessingType.CannyEdge);

		ImageJFunctions.show(inputimg).setTitle("Preprocessed image");

		ArrayList<Lineobjects> linelist = new ArrayList<Lineobjects>(biginputimg.numDimensions());
		Img<IntType> Intimg = new ArrayImgFactory<IntType>().create(biginputimg, new IntType());
		Pair<Img<IntType>, ArrayList<Lineobjects>> linepair = new Pair<Img<IntType>, ArrayList<Lineobjects>>(Intimg,
				linelist);

		// Do watershedding and Hough
		
		// Declare minimum length of the line(in pixels) to be detected
				double minlength = 10;
				
		linepair = PerformWatershedding.DowatersheddingandHough(biginputimg, inputimg, minlength);

		// Overlay detected lines on the image

		PointSampleList<FloatType> centroidlist = new PointSampleList<FloatType>(n);

		

		final int ndims = biginputimg.numDimensions();
		double[] final_param = new double[2 * ndims + 2];
		double[] psf = new double[ndims];
		final long radius = 2;
		psf[0] = 1.7;
		psf[1] = 1.54;
		// Input the psf-sigma here to be used for convolving Gaussians on a
		// line, will not change during iteration.

		
		

		ArrayList<double[]> totalgausslist = new ArrayList<double[]>();

		// Get a rough reconstruction of the line and the list of centroids where psf of the image has to be convolved
		
		OverlayLines.GetAlllines(imgout, biginputimg, linepair.fst, centroidlist, linepair.snd, radius);

		ImageJFunctions.show(imgout).setTitle("Rough-Reconstruction");


		// Input the image on which you want to do the fitting, along with the labelled image
		
		LengthDetection MTlength = new LengthDetection(biginputimg, linepair.fst);

		final double[] listpoint = new double[n];

		// Choose the noise level of the image, 0 for pre-processed image and >0 for original image
		
		
		Cursor<FloatType> listcursor = centroidlist.localizingCursor();
		while (listcursor.hasNext()) {
			listcursor.fwd();
			listcursor.localize(listpoint);
			final_param = MTlength.Getfinalparam(listcursor, radius, psf);

			// Choosing values above the nosie level of the image 
			
			if ( final_param[3] > 0 && final_param[4] > 0){
		
				totalgausslist.add(final_param);

				System.out.println(" Amplitude: " + final_param[0] + " " + "Mean X: " + final_param[1] + " "
						+ "Mean Y: " + final_param[2] + " " + "SigmaX: " + Math.sqrt(1.0/final_param[3]) + " " + "SigmaY: "
						+ Math.sqrt(1.0/final_param[4]) + " " + "Poisson: "+ final_param[5]);
		}
		}

		// Draw the Gaussian convolved line fitted with the original data
		PushCurves.DrawDetectedGaussians(gaussimg, totalgausslist);
		ImageJFunctions.show(gaussimg).setTitle("Iterated Result");

	}
}
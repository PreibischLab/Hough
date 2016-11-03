package velocityanalyser;

import java.io.File;
import java.util.ArrayList;

import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.SimpleWeightedGraph;

import com.sun.tools.javac.util.Pair;

import drawandOverlay.DisplayGraph;
import drawandOverlay.DisplaysubGraphend;
import drawandOverlay.DisplaysubGraphstart;
import drawandOverlay.OverlayLines;
import drawandOverlay.PushCurves;
import graphconstructs.Staticproperties;
import houghandWatershed.HoughTransformandwatershed2D;
import houghandWatershed.WatershedDistimg;
import ij.ImageJ;
import ij.ImagePlus;
import ij.io.FileSaver;
import labeledObjects.LabelledImg;
import labeledObjects.Lineobjects;
import labeledObjects.Simpleobject;
import labeledObjects.Subgraphs;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.stats.Normalize;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.IntervalView;
import net.imglib2.view.Views;
import peakFitter.SubpixelLength;
import peakFitter.SubpixelVelocity;
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
				.openAs32Bit(
				//		new File("../res/10frame_moving.tif"),
				//			new ArrayImgFactory<FloatType>());
						
						new File("/Users/varunkapoor/Hough/src/main/resources/2015-01-14_Seeds-1.tiff"),
						new ArrayImgFactory<FloatType>());
		int ndims = img.numDimensions();

		// Normalize the intensity of the whole stack to be between min and max
		// value specified here
		new Normalize();

		FloatType minval = new FloatType(0);
		FloatType maxval = new FloatType(1);
		Normalize.normalize(Views.iterable(img), minval, maxval);

		// Declare all the constants needed by the program here:

		final double[] psf = { 1.65, 1.47 };
		//final double[] psf = { 1.7, 1.8 };
		
		final long radius = (long) Math.ceil(Math.sqrt(psf[0] * psf[0] + psf[1] * psf[1]));
		final int minlength = 4;

		// Show the stack
		ImagePlus impstart = ImageJFunctions.show(img);
		ImagePlus impend = ImageJFunctions.show(img);
		ArrayList<ArrayList<Staticproperties>> Allstartandend = new ArrayList<ArrayList<Staticproperties>>();

		
		if (ndims == 2 ){
			
			
			RandomAccessibleInterval<FloatType> inputimg = new ArrayImgFactory<FloatType>().create(img,
					new FloatType());
			RandomAccessibleInterval<FloatType> preinputimg = new ArrayImgFactory<FloatType>().create(img,
					new FloatType());
			RandomAccessibleInterval<FloatType> imgout = new ArrayImgFactory<FloatType>().create(img,
					new FloatType());
			RandomAccessibleInterval<FloatType> gaussimg = new ArrayImgFactory<FloatType>().create(img,
					new FloatType());

			
			// Preprocess image using Median Filter and suppress background
			final MedianFilter2D<FloatType> medfilter = new MedianFilter2D<FloatType>( img, 1 );
			medfilter.process();
			preinputimg = medfilter.getResult();
			inputimg = Kernels.Supressthresh(preinputimg);
			Normalize.normalize(Views.iterable(inputimg), minval, maxval);
			

			ImageJFunctions.show(inputimg).setTitle("Preprocessed image");

			// Create the Bit image for distance transform and seed image for watershedding
			final Float ThresholdValue = GlobalThresholding.AutomaticThresholding(inputimg);
			RandomAccessibleInterval<BitType> bitimg = new ArrayImgFactory<BitType>().create(inputimg, new BitType());
			GetLocalmaxmin.ThresholdingBit(inputimg, bitimg, ThresholdValue);

		
			
			
			// Do watershedding and Hough
			System.out.println("Doing Hough transform in labels: ");

			HoughTransformandwatershed2D Houghobject = new HoughTransformandwatershed2D(inputimg, bitimg, minlength);
			
			Houghobject.checkInput();
			Houghobject.process();
			Pair<RandomAccessibleInterval<IntType>, ArrayList<Lineobjects>> linepair  = Houghobject.getResult();


			// Overlay detected lines on the image
			final ArrayList<Simpleobject> simpleobject = new ArrayList<Simpleobject>();
			OverlayLines.GetAlllines(imgout,  linepair.fst, linepair.snd, simpleobject, radius);

			ImageJFunctions.show(imgout).setTitle("Rough-Reconstruction");

			// Input the image on which you want to do the fitting (original noisy image, along with the
			// labelled image (watershedded image) 

			SubpixelLength MTline = new SubpixelLength(img, linepair.fst, simpleobject, psf, minlength);
			MTline.checkInput();
			MTline.process();
			ArrayList<double[]> final_paramlist = MTline.getResult();
			
			// Draw the detected lines
			PushCurves.DrawallLine(gaussimg, final_paramlist, psf);
			ImageJFunctions.show(gaussimg).setTitle("Exact-line");

		}
			
		
		
		
		
		if (ndims > 2){
		// Do Hough transform on the First seed image

		IntervalView<FloatType> groundframe = Views.hyperSlice(img, ndims - 1, 0);

		
		
		
		System.out.println("Applying Median filter to the first image.");
		// Preprocess image using Median Filter and suppress background
				final MedianFilter2D<FloatType> medfilter = new MedianFilter2D<FloatType>( groundframe, 1);
				medfilter.process();
				RandomAccessibleInterval<FloatType> groundframepre = medfilter.getResult();
				Normalize.normalize(Views.iterable(groundframepre), minval, maxval);
		System.out.println("Median Filter applied sucessfully.");
		
		
		// for thresholding extremly noisy data, if non noisy set the value to 1
		RandomAccessibleInterval<FloatType> inputimg = Kernels.Supressthresh(groundframepre);
				
		ImageJFunctions.show(inputimg);
		
		final Float ThresholdValue = GlobalThresholding.AutomaticThresholding(inputimg);
		RandomAccessibleInterval<BitType> bitimg = new ArrayImgFactory<BitType>().create(inputimg, new BitType());
		GetLocalmaxmin.ThresholdingBit(inputimg, bitimg, ThresholdValue);
		System.out.println("Doing Hough transform in labels: ");

        HoughTransformandwatershed2D Houghobject = new HoughTransformandwatershed2D(inputimg, bitimg, minlength);
		Houghobject.checkInput();
		Houghobject.process();
		Pair<RandomAccessibleInterval<IntType>, ArrayList<Lineobjects>> linepair  = Houghobject.getResult();

		

		// Display the Hough detection
		RandomAccessibleInterval<FloatType> imgout = new ArrayImgFactory<FloatType>().create(groundframe,
				new FloatType());
		final ArrayList<Simpleobject> simpleobject = new ArrayList<Simpleobject>();

		OverlayLines.GetAlllines(imgout, linepair.fst, linepair.snd, simpleobject, radius);

		ImageJFunctions.show(imgout).setTitle("Rough-Reconstruction");


		RandomAccessibleInterval<FloatType> gaussimg = new ArrayImgFactory<FloatType>().create(groundframe,
				new FloatType());

		SubpixelLength MTline = new SubpixelLength(groundframe, linepair.fst, simpleobject, psf, minlength);
		MTline.checkInput();
		MTline.process();
		ArrayList<double[]> final_paramlist = MTline.getResult();

		

		// Draw detected lines from the seed image
		PushCurves.DrawallLine(gaussimg, final_paramlist, psf);
		ImageJFunctions.show(gaussimg);
		

		ArrayList<double[]> PrevFrameparam = final_paramlist;

		// Now start tracking the moving ends of the Microtubule and make
		// seperate graph for both ends

		for (int frame = 1; frame < img.dimension(ndims - 1); ++frame) {

			IntervalView<FloatType> currentframe = Views.hyperSlice(img, ndims - 1, frame);
			
			System.out.println("Applying Median filter to current frame.");
			// Preprocess image using Median Filter and suppress background
					final MedianFilter2D<FloatType> medfiltercurr = new MedianFilter2D<FloatType>( currentframe, 1);
					medfiltercurr.process();
					RandomAccessibleInterval<FloatType> precurrent = medfiltercurr.getResult();
					Normalize.normalize(Views.iterable(precurrent), minval, maxval);
			System.out.println("Median Filter applied sucessfully.");
		
			
			RandomAccessibleInterval<FloatType> inputimgpre = Kernels.Supressthresh(precurrent);
					//Kernels.CannyEdge(precurrent, psf, true, Heavythresholdvalue);
	
			final Float currThresholdValue = GlobalThresholding.AutomaticThresholding(inputimgpre);
			RandomAccessibleInterval<BitType> currbitimg = new ArrayImgFactory<BitType>().create(inputimgpre, new BitType());
			GetLocalmaxmin.ThresholdingBit(inputimgpre, currbitimg, currThresholdValue);
			
			WatershedDistimg Watershedobject = new WatershedDistimg(inputimgpre, currbitimg);
			Watershedobject.checkInput();
			Watershedobject.process();
			RandomAccessibleInterval<IntType> watershedimg = Watershedobject.getResult();
			
			
			final SubpixelVelocity growthtracker = new SubpixelVelocity(currentframe, watershedimg, PrevFrameparam, psf, frame);
			growthtracker.checkInput();
			growthtracker.process();
			ArrayList<double[]> NewFrameparam = growthtracker.getResult();
			ArrayList<Staticproperties> StateVectors = growthtracker.getStateVectors();
			// Update the list of line parameters with the current frame
			// detectionx
			PrevFrameparam = NewFrameparam;
			// Append the object static properties with the current frame
			// detection
			Allstartandend.add(StateVectors);

			// Draw the lines detected in the current frame
			RandomAccessibleInterval<FloatType> newgaussimg = new ArrayImgFactory<FloatType>().create(groundframe,
					new FloatType());
			PushCurves.DrawallLine(newgaussimg, NewFrameparam, psf);
			ImageJFunctions.show(newgaussimg);
			

		}

		// Overlay the graphs on the stack

		// Make graph to track the start and the end point

		final int maxframe = (int) img.dimension(ndims - 1);
		final Trackstart trackerstart = new Trackstart(Allstartandend, maxframe);
		final Trackend trackerend = new Trackend(Allstartandend, maxframe);
		trackerstart.process();
		SimpleWeightedGraph<double[], DefaultWeightedEdge> graphstart = trackerstart.getResult();
		ArrayList<Subgraphs> subgraphstart = trackerstart.getFramedgraph();
		

		DisplaysubGraphstart displaytrackstart = new DisplaysubGraphstart(impstart, subgraphstart);
		displaytrackstart.getImp();
		impstart.draw();
		
		trackerend.process();
		SimpleWeightedGraph<double[], DefaultWeightedEdge> graphend = trackerend.getResult();
		ArrayList<Subgraphs> subgraphend = trackerend.getFramedgraph();
		

		DisplaysubGraphend displaytrackend = new DisplaysubGraphend(impend, subgraphend);
		displaytrackend.getImp();
		impend.draw();
		
	}
	}
}

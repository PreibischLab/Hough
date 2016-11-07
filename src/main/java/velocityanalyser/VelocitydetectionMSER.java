package velocityanalyser;

import java.io.File;
import java.util.ArrayList;

import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.SimpleWeightedGraph;

import drawandOverlay.DisplaysubGraphend;
import drawandOverlay.DisplaysubGraphstart;
import drawandOverlay.OverlayLines;
import drawandOverlay.PushCurves;
import getRoi.RoiforMSER;
import graphconstructs.Staticproperties;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.gui.Overlay;
import labeledObjects.LabelledImg;
import labeledObjects.Simpleobject;
import labeledObjects.Subgraphs;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.stats.Normalize;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.IntervalView;
import net.imglib2.view.Views;
import peakFitter.SubpixelLengthMSER;
import peakFitter.SubpixelVelocityMSER;
import preProcessing.Biggify;
import preProcessing.MedianFilter2D;

public class VelocitydetectionMSER {

	public static void main(String[] args) throws Exception {

		/***
		 * MSER and optionally Hough Transform to detect Microtubules and track the growth at
		 * Sub-pixel accuracy. Optimizers used: Levenberg-Marqurat solver and
		 * Weighted centre of mass fits. Program reqires PSF of the microscope
		 * to be computed and analysed and takes the determined Sigmas as the
		 * input. @ Varun Kapoor
		 */

		new ImageJ();

		// Load the stack of images
		 RandomAccessibleInterval<FloatType> img = util.ImgLib2Util
				.openAs32Bit(
						new File("../res/2016-09-28_bovine_cy5seeds_cy3tub_6uM_seeds.tif"),
							new ArrayImgFactory<FloatType>());
					
					//	new File("../res/2015-01-14_Seeds-1.tiff"),
					//	new ArrayImgFactory<FloatType>());
						
					//	new File("../res/10frame_moving.tif"),
					//	new ArrayImgFactory<FloatType>());
		int ndims = img.numDimensions();
		
		// Normalize the intensity of the whole stack to be between min and max
		// value specified here
		new Normalize();

		final boolean darktoBright = false;
		FloatType minval = new FloatType(0);
		FloatType maxval = new FloatType(1);
		Normalize.normalize(Views.iterable(img), minval, maxval);
		final double[] psf = { 1.65, 1.47 };
		final int extendBorderpixels = 10;
		img = Biggify.biggifyimage(img, extendBorderpixels);
		// Declare all the constants needed by the program here:

		
		//final double[] psf = { 1.4, 1.5 };
		
		// minimum length of the lines to be detected, the smallest possible number is 2.
		final double minlength = 5;

		// Show the stack
		ImagePlus impstart = ImageJFunctions.show(img);
		ImagePlus impend = ImageJFunctions.show(img);
		ArrayList<ArrayList<Staticproperties>> Allstartandend = new ArrayList<ArrayList<Staticproperties>>();
		final int delta = 5;
		final long minSize = 10;
		final long maxSize = Long.MAX_VALUE;
		final double maxVar = 0.2;
		final double minDiversity = 0;
		

		if (ndims == 2 ){
			
			
			RandomAccessibleInterval<FloatType> inputimg = new ArrayImgFactory<FloatType>().create(img,
					new FloatType());
			RandomAccessibleInterval<FloatType> preinputimg = new ArrayImgFactory<FloatType>().create(img,
					new FloatType());
			
			
			    
			
			// Preprocess image using Median Filter and suppress background
			final MedianFilter2D<FloatType> medfilter = new MedianFilter2D<FloatType>( img, 2 );
			medfilter.process();
			inputimg = medfilter.getResult();
		
			Normalize.normalize(Views.iterable(inputimg), minval, maxval);
			
			 
			
	        
			ImageJFunctions.show(inputimg).setTitle("Preprocessed extended image");
	        
			RandomAccessibleInterval<FloatType> imgout = new ArrayImgFactory<FloatType>().create(img,
					new FloatType());
			RandomAccessibleInterval<FloatType> gaussimg = new ArrayImgFactory<FloatType>().create(img,
					new FloatType());
			
			
			RoiforMSER Roiobject = new RoiforMSER(inputimg, img, delta, minSize, maxSize, maxVar, minDiversity, darktoBright);
			Roiobject.checkInput();
			Roiobject.process();
			ArrayList<LabelledImg> arrayimg = Roiobject.getResult();
			Overlay ov = Roiobject.getOverlay();
			ImageJFunctions.show(inputimg);
			ImagePlus imp = IJ.getImage();
			//ImagePlus imp = ImageJFunctions.wrap(inputimg, "curr");
			imp.setOverlay(ov);
			/**
			 * To see the overlay 
			 * 
			 * ImageJFunctions.show(inputimg);
			 * ImagePlus imp = IJ.getImage();
			 * 
			 */
			// Overlay detected lines on the image
				final ArrayList<Simpleobject> simpleobject = new ArrayList<Simpleobject>();
				OverlayLines.Getmserlines(imgout, arrayimg, simpleobject);

						ImageJFunctions.show(imgout).setTitle("Rough-Reconstruction");
						
						SubpixelLengthMSER MTline = new SubpixelLengthMSER(img, arrayimg, simpleobject, psf, minlength);
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

		RandomAccessibleInterval<FloatType> inputimg = new ArrayImgFactory<FloatType>().create(groundframe,
				new FloatType());
		
		
		
	
		
		System.out.println("Applying Median filter to the first image.");
		// Preprocess image using Median Filter and suppress background
				final MedianFilter2D<FloatType> medfilter = new MedianFilter2D<FloatType>( groundframe, 1);
				medfilter.process();
				inputimg = medfilter.getResult();
				Normalize.normalize(Views.iterable(inputimg), minval, maxval);
				
		System.out.println("Median Filter applied sucessfully.");
		
		
		ImageJFunctions.show(inputimg);
		RandomAccessibleInterval<FloatType> imgout = new ArrayImgFactory<FloatType>().create(groundframe,
				new FloatType());
		RandomAccessibleInterval<FloatType> gaussimg = new ArrayImgFactory<FloatType>().create(groundframe,
				new FloatType());
		
		System.out.println("Running MSER: ");

		RoiforMSER Roiobject = new RoiforMSER(inputimg, groundframe, delta, minSize, maxSize, maxVar, minDiversity, false);
		Roiobject.checkInput();
		Roiobject.process();
		ArrayList<LabelledImg> arrayimg = Roiobject.getResult();
		Overlay ov = Roiobject.getOverlay();
		ImageJFunctions.show(inputimg);
		ImagePlus imp = IJ.getImage();
		//ImagePlus imp = ImageJFunctions.wrap(inputimg, "curr");
		/**
		 * To see the overlay 
		 * 
		 * ImageJFunctions.show(inputimg);
		 * ImagePlus imp = IJ.getImage();
		 * 
		 */
		
		imp.setOverlay(ov);
		// Overlay detected lines on the image
		final ArrayList<Simpleobject> simpleobject = new ArrayList<Simpleobject>();
		OverlayLines.Getmserlines(imgout, arrayimg, simpleobject);

				ImageJFunctions.show(imgout).setTitle("Rough-Reconstruction");
				
				SubpixelLengthMSER MTline = new SubpixelLengthMSER(groundframe, arrayimg, simpleobject, psf, minlength);
				MTline.checkInput();
				MTline.process();
				ArrayList<double[]> final_paramlist = MTline.getResult();
				
				// Draw the detected lines
				PushCurves.DrawallLine(gaussimg, final_paramlist, psf);
				ImageJFunctions.show(gaussimg).setTitle("Exact-line");
				
	

		ArrayList<double[]> PrevFrameparam = final_paramlist;

		// Now start tracking the moving ends of the Microtubule and make
		// seperate graph for both ends

		for (int frame = 1; frame < img.dimension(ndims - 1); ++frame) {

			IntervalView<FloatType> currentframe = Views.hyperSlice(img, ndims - 1, frame);
			
			
			System.out.println("Applying Median filter to current frame.");
			// Preprocess image using Median Filter and suppress background
					final MedianFilter2D<FloatType> medfiltercurr = new MedianFilter2D<FloatType>( currentframe, 1);
					medfiltercurr.process();
					RandomAccessibleInterval<FloatType> inputimgpre = medfiltercurr.getResult();
					Normalize.normalize(Views.iterable(inputimgpre), minval, maxval);

					ImageJFunctions.show(inputimgpre);
			System.out.println("Median Filter applied sucessfully.");
		
			
			
			
			
			RoiforMSER Roiobjectframe = new RoiforMSER(inputimgpre, currentframe, delta, minSize, maxSize, maxVar, minDiversity, false);
			Roiobjectframe.checkInput();
			Roiobjectframe.process();
			ArrayList<LabelledImg> arrayimgframe = Roiobjectframe.getResult();
			Overlay ovframe = Roiobjectframe.getOverlay();
		
			ImagePlus impframe = IJ.getImage();
			//ImagePlus impframe = ImageJFunctions.wrap(inputimgpre, "curr");
			impframe.setOverlay(ovframe);
			/**
			 * To see the overlay 
			 * 
			 * ImageJFunctions.show(inputimgpre);
			 * ImagePlus imp = IJ.getImage();
			 * 
			 */
			
			
			final SubpixelVelocityMSER growthtracker = new SubpixelVelocityMSER(currentframe, arrayimgframe, PrevFrameparam, psf, frame);
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
			ImageJFunctions.show(newgaussimg).setTitle("Exact-line");
			

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

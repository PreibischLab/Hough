package houghandWatershed;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import com.sun.tools.javac.util.Pair;

import drawandOverlay.OverlayLines;
import drawandOverlay.PushCurves;
import ij.ImageJ;
import labeledObjects.Finalobject;
import labeledObjects.Indexedlength;
import labeledObjects.LabelMax;
import labeledObjects.Lineobjects;
import labeledObjects.PreFinalobject;
import net.imglib2.Cursor;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealCursor;
import net.imglib2.RealPoint;
import net.imglib2.RealPointSampleList;
import net.imglib2.algorithm.stats.Normalize;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import pSF.Boundingbox.Objectproperties;
import peakFitter.LengthDetection;
import preProcessing.GetLocalmaxmin.IntensityType;
import preProcessing.Kernels;
import preProcessing.Kernels.ProcessingType;
import util.ImgLib2Util;

public class HoughpostWater {

	public static void main(String[] args) throws Exception {

		RandomAccessibleInterval<FloatType> biginputimg = ImgLib2Util
				.openAs32Bit(new File("src/main/resources/Fake_datasmall.tif"));
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
		

		// Initialize empty images to be used later
	
		RandomAccessibleInterval<FloatType> inputimg = new ArrayImgFactory<FloatType>().create(biginputimg,
				new FloatType());
		RandomAccessibleInterval<FloatType> preinputimg = new ArrayImgFactory<FloatType>().create(biginputimg,
				new FloatType());
		RandomAccessibleInterval<FloatType> imgout = new ArrayImgFactory<FloatType>().create(biginputimg,
				new FloatType());
		RandomAccessibleInterval<FloatType> gaussimg = new ArrayImgFactory<FloatType>().create(biginputimg,
				new FloatType());
		

		// Preprocess image
		
		
		
		preinputimg = Kernels.Meanfilterandsupress(biginputimg, 1.0);
		inputimg = Kernels.Preprocess(preinputimg, ProcessingType.CannyEdge);

		ImageJFunctions.show(inputimg).setTitle("Preprocessed image");

		ArrayList<Lineobjects> linelist = new ArrayList<Lineobjects>(biginputimg.numDimensions());
		Img<IntType> Intimg = new ArrayImgFactory<IntType>().create(biginputimg, new IntType());
		Pair<Img<IntType>, ArrayList<Lineobjects>> linepair = new Pair<Img<IntType>, ArrayList<Lineobjects>>(Intimg,
				linelist);

		// Do watershedding and Hough
		
		// Declare minimum length of the line(in pixels) to be detected
				double minlength = 0;
				
		linepair = PerformWatershedding.DowatersheddingandHough(biginputimg, inputimg, minlength);

		// Overlay detected lines on the image

		PointSampleList<FloatType> centroidlist = new PointSampleList<FloatType>(n);

		

		final int ndims = biginputimg.numDimensions();
		double[] final_param = new double[2 * ndims + 2];
		
		double[] noise_param = new double[ndims];
		double[] psf = new double[ndims];
		
		final double SNR = 4000/240;
		psf[0] = 1.7;
		psf[1] = 1.8;
		final long radius = (long) Math.ceil(  Math.sqrt( psf[0] * psf[0] +  psf[1] * psf[1]));
		// Input the psf-sigma here to be used for convolving Gaussians on a
		// line, will not change during iteration.

		
		

		ArrayList<double[]> totalgausslist = new ArrayList<double[]>();

		// Rough reconstruction of line after Hough detection

		final ArrayList<PreFinalobject> prefinalparamlist = new ArrayList<PreFinalobject>();
		
		OverlayLines.GetAlllines(imgout, biginputimg, linepair.fst, centroidlist,prefinalparamlist, linepair.snd, radius);

		ImageJFunctions.show(imgout).setTitle("Rough-Reconstruction");


		// Input the image on which you want to do the fitting, along with the labelled image
		
		LengthDetection MTlength = new LengthDetection(biginputimg, linepair.fst);
		
		
		final ArrayList<Finalobject> finalparamlist = new ArrayList<Finalobject>();
		
		for (int index = 0; index < prefinalparamlist.size(); ++index){
			
			
			// Do gradient descent to improve the Hough detected lines
			final_param = MTlength.Getfinalparam(prefinalparamlist.get(index).centroid, radius, psf);
			noise_param = MTlength.Getnoiseparam(prefinalparamlist.get(index).centroid, radius);
			
		if (final_param[0] > 1.0E-10){
			
				
				
				final double[] newmeans = {final_param[1], final_param[2]};
				 
				RealPoint newcentroid = new RealPoint(newmeans);
				
				// Update the Intensity, Mean and Variance of the lines
				Finalobject finalline = new Finalobject(prefinalparamlist.get(index).Label, newcentroid, final_param[0],
						1.0/Math.sqrt(final_param[3]), 1.0/Math.sqrt(final_param[4]),
						prefinalparamlist.get(index).slope, prefinalparamlist.get(index).intercept);
				
				finalparamlist.add(finalline);
			/*	
				System.out.println(" Amp: " + finalparamlist.get(index).Intensity + " " + "Mu X: " 
				+ finalparamlist.get(index).centroid.getDoublePosition(0) + " "
						+ "Mu Y: " + finalparamlist.get(index).centroid.getDoublePosition(1) +
						" " + "Sig X: " + finalparamlist.get(index).sigmaX + " " + "Sig Y: "
						+ finalparamlist.get(index).sigmaY );
				
			
			*/
		}
		
		}
		ArrayList<Indexedlength> finallength = new ArrayList<Indexedlength>();
		
		
		// Since the centroid positions changed, update the slope and intercept of Hough lines
		 ArrayList<Finalobject> updateparamlist = new ArrayList<Finalobject>();
		 ArrayList<Finalobject> correctparamlist = new ArrayList<Finalobject>();
		
		
		updateparamlist = MTlength.Updateslopeandintercept(finalparamlist);
		ArrayList<LabelMax> labelmaxlist = new ArrayList<LabelMax>();
		correctparamlist = MTlength.Removepoints(updateparamlist, labelmaxlist);
	
		for (int index = 0; index < correctparamlist.size(); ++index){
			
			double[] reduced_param = new double[2 * ndims + 2];
			System.out.println(" Label : " + correctparamlist.get(index).Label  
					+ " Amp: " + correctparamlist.get(index).Intensity + " " + "Mu X: " + 
			correctparamlist.get(index).centroid.getDoublePosition(0) + " "
					+ "Mu Y: " + correctparamlist.get(index).centroid.getDoublePosition(1) 
					+ " " + "Sig X: " + correctparamlist.get(index).sigmaX+ " " + "Sig Y: "
					+  correctparamlist.get(index).sigmaY  );
			
			reduced_param[0] = correctparamlist.get(index).Intensity;
		    
			for (int d = 0; d < n; ++d){
				reduced_param[d + 1] = correctparamlist.get(index).centroid.getDoublePosition(d); 
			}
           
			reduced_param[3] =   1.0 / Math.pow(correctparamlist.get(index).sigmaX, 2);
			
			reduced_param[4] =   1.0 / Math.pow(correctparamlist.get(index).sigmaY, 2);
			
			
			
			totalgausslist.add(reduced_param);
			
		}
	
		
		// Draw the detected and iterated lines
		PushCurves.DrawDetectedGaussians(gaussimg, totalgausslist);
		ImageJFunctions.show(gaussimg).setTitle("Iterated Result");
		
		// Determine starting and end points by giving the start and end positions along the line and then doing a half Gaussian mask fit
		
		MTlength.Returnlengths(correctparamlist , finallength, labelmaxlist, psf);

		
		
		
		for (int index = 0; index< finallength.size(); ++index){
			
		if (finallength.get(index).length > 0 )	
			try {
	            FileWriter writer = new FileWriter("finallengths.txt", true);
	            writer.write("Label: " + finallength.get(index).Label + " " +
	           		 "Length: "+ finallength.get(index).length + " " + "Slope: " + finallength.get(index).slope +  
	        		 " Intercept :" + finallength.get(index).intercept + " StartposX: " + finallength.get(index).startpos[0] 
	        				 +" StartPosY: " + finallength.get(index).startpos[1]+ " EndposX: " + finallength.get(index).endpos[0] 
	    	        				 +" EndPosY: " + finallength.get(index).endpos[1] );
	            writer.write("\r\n"); 
	            writer.close();
	        } catch (IOException e) {
	            e.printStackTrace();
	        }
			
		
		}
		
		

	}
}
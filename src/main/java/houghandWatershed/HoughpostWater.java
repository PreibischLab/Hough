package houghandWatershed;

import java.io.File;
import java.util.ArrayList;

import com.sun.tools.javac.util.Pair;

import drawandOverlay.OverlayLines;
import drawandOverlay.PushCurves;
import ij.ImageJ;
import labeledObjects.Lineobjects;
import net.imglib2.Cursor;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccessibleInterval;
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
				.openAs32Bit(new File("src/main/resources/Fresh_data/psf_488_01.tif"));
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
			sigma[d] = 1;

		// Initialize empty images to be used later
		RandomAccessibleInterval<FloatType> inputimg = new ArrayImgFactory<FloatType>().create(biginputimg,
				new FloatType());
		RandomAccessibleInterval<FloatType> imgout = new ArrayImgFactory<FloatType>().create(biginputimg,
				new FloatType());
		RandomAccessibleInterval<FloatType> localmaximgout = new ArrayImgFactory<FloatType>().create(biginputimg,
				new FloatType());
		RandomAccessibleInterval<FloatType> gaussimg = new ArrayImgFactory<FloatType>().create(biginputimg,
				new FloatType());
		
		// Preprocess image
		//inputimg = Kernels.Preprocess(biginputimg, ProcessingType.Meanfilter);
		 inputimg = Kernels.Preprocess(biginputimg, ProcessingType.Meanfilter);

		ImageJFunctions.show(inputimg).setTitle("Preprocessed image");

		ArrayList<Lineobjects> linelist = new ArrayList<Lineobjects>(biginputimg.numDimensions());
		Img<IntType> Intimg = new ArrayImgFactory<IntType>().create(biginputimg, new IntType());
		Pair<Img<IntType>, ArrayList<Lineobjects>> linepair = new Pair<Img<IntType>, ArrayList<Lineobjects>>(Intimg,
				linelist);

		// Do watershedding and Hough
		linepair = PerformWatershedding.DowatersheddingandHough(biginputimg, inputimg);

		// Overlay detected lines on the image
		
		PointSampleList<FloatType> centroidlist = new PointSampleList<FloatType>(n);
		// Model lines
		localmaximgout = OverlayLines.GetAlllines(imgout,linepair.fst, linepair.snd);

		ImageJFunctions.show(imgout);
		ImageJFunctions.show(linepair.fst);
		ImageJFunctions.show(localmaximgout);
		
		PushCurves.MakeHTguess(localmaximgout, centroidlist);
		
		// Do Gaussian Fit
		final int ndims = biginputimg.numDimensions();
		 double[] final_param= new double[2*ndims+1];
		 final double [] point_spread_sigma = new double[ndims];
		 // Input the psf-sigma here to be used as a replacment for very large sigma in the solver
		 for (int d = 0; d < ndims; ++d)
			 point_spread_sigma[d] = 1;
		 ArrayList<double[]> totalgausslist = new ArrayList<double[]>();
		 
			
			
	
		 LengthDetection MTlength = new LengthDetection(inputimg,linepair.fst);
		
		 Cursor<FloatType> listcursor = centroidlist.localizingCursor();
		 while(listcursor.hasNext()){
			 listcursor.fwd();
			 final_param = MTlength.Getfinalparam(listcursor, listcursor.get().get(), point_spread_sigma);
			 
			 if (final_param[3] > 0  && final_param[4]>0 ){
			 totalgausslist.add(final_param);
			 
			 System.out.println( " Amplitude: " + final_param[0] + " " + "Mean X: "
						+ final_param[1] + " " + "Mean Y: " + final_param[2] + " " + "1/SigmaX^2: "
						+ final_param[3] + " " + "1/SigmaY^2: "
						+ final_param[4]);
			 }
		 }
		 
		
			 
		 PushCurves.DrawDetectedGaussians(gaussimg, totalgausslist);	
			ImageJFunctions.show(gaussimg).setTitle("Iterated Result");
			
		
	}
}
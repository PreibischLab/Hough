package varun;

import java.io.File;
import java.util.ArrayList;

import com.sun.tools.javac.util.Pair;

import ij.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.Interval;
import net.imglib2.Point;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.gauss3.Gauss3;
import net.imglib2.algorithm.localextrema.RefinedPeak;
import net.imglib2.algorithm.stats.Normalize;
import net.imglib2.exception.IncompatibleTypeException;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Intervals;
import net.imglib2.view.Views;
import util.ImgLib2Util;
import varun.GetLocalmaxmin.IntensityType;
import varun.Kernels.ProcessingType;
import varun.LengthDetection.Labelparam;
import varun.OverlayLines.Simulatedline;
import varun.PerformWatershedding.Lineobjects;

public class HoughpostWater {

	public static void main(String[] args) throws Exception {

		RandomAccessibleInterval<FloatType> biginputimg = ImgLib2Util
				.openAs32Bit(new File("src/main/resources/small_mt.tif"));
		// small_mt.tif image to be used for testing
		// 2015-01-14_Seeds-1.tiff for actual
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

		RandomAccessibleInterval<FloatType> inputimg = new ArrayImgFactory<FloatType>().create(biginputimg,
				new FloatType());
		RandomAccessibleInterval<FloatType> imgout = new ArrayImgFactory<FloatType>().create(biginputimg,
				new FloatType());
		RandomAccessibleInterval<FloatType> localmaximgout = new ArrayImgFactory<FloatType>().create(biginputimg,
				new FloatType());
		RandomAccessibleInterval<FloatType> gaussimg = new ArrayImgFactory<FloatType>().create(biginputimg,
				new FloatType());
		RandomAccessibleInterval<FloatType> initialgaussimg = new ArrayImgFactory<FloatType>().create(biginputimg,
				new FloatType());
		// Preprocess image
		inputimg = Kernels.Preprocess(biginputimg, ProcessingType.Meanfilter);
		// inputimg = Kernels.Preprocess(tmpinputimg, ProcessingType.NaiveEdge);

		ImageJFunctions.show(inputimg).setTitle("Preprocessed image");

		ArrayList<Lineobjects> linelist = new ArrayList<Lineobjects>(biginputimg.numDimensions());
		Img<IntType> Intimg = new ArrayImgFactory<IntType>().create(biginputimg, new IntType());
		Pair<Img<IntType>, ArrayList<Lineobjects>> linepair = new Pair<Img<IntType>, ArrayList<Lineobjects>>(Intimg,
				linelist);

		// Do watershedding and Hough
		linepair = PerformWatershedding.DowatersheddingandHough(biginputimg, inputimg);

		// Overlay detected lines on the image
		ArrayList<Simulatedline> simline = new ArrayList<Simulatedline>();
		
		ArrayList<Labelparam> params = new ArrayList<Labelparam>();
		ArrayList<Labelparam> lineparams = new ArrayList<Labelparam>();
		PointSampleList<FloatType> centroidlist = new PointSampleList<FloatType>(n);
		// Model lines
		OverlayLines.GetAlllines(imgout, simline, params, linepair.fst, linepair.snd);

		ImageJFunctions.show(imgout);
		localmaximgout = GetLocalmaxmin.FindandDisplayLocalMaxima(imgout,
				IntensityType.Original, sigma);
		ImageJFunctions.show(localmaximgout);
		
		PushCurves.MakeHTguess(localmaximgout, centroidlist);
		
		// Do Gaussian Fit
		final int ndims = biginputimg.numDimensions();
		 double[] final_param= new double[2*ndims+1];
		 ArrayList<double[]> totalgausslist = new ArrayList<double[]>();
		 final double[] psf = new double[ndims];
			for (int d = 0; d < ndims; ++d)
				psf[d] = 2;
			double[] pad_size = new double[ndims];
			
			for (int d = 0; d < ndims; d++) {
				pad_size[d] =  (2 * psf[d] + 1);
			}
			
			


			int[] size = new int[ndims];
			
			for (int i = 0; i < ndims; i++) {
			size[i] = (int) (2 * pad_size[i] + 1);
			}
			
			
			int span = 0;
			for (int i = 0; i < ndims; i++) {
				span = size[i]/ndims;
				}
			
			Interval interval = Intervals.expand(inputimg, -span);

			// create a view on the source with this interval
			inputimg = Views.interval(inputimg, interval);
			
	
		 LengthDetection MTlength = new LengthDetection(inputimg, linepair.fst, centroidlist);
		
		 Cursor<FloatType> listcursor = centroidlist.localizingCursor();
		 while(listcursor.hasNext()){
			 listcursor.fwd();
			 final_param = MTlength.Getfinalparam(psf, listcursor, listcursor.get().get());
			 totalgausslist.add(final_param);
			 System.out.println( " Amplitude: " + final_param[0] + " " + "Mean X: "
						+ final_param[1] + " " + "Mean Y: " + final_param[2] + " " + "SigmaX: "
						+ 1.0 / Math.sqrt(final_param[3]) + " " + "SigmaY: "
						+ 1.0 / Math.sqrt(final_param[4]));
		 }
		 PushCurves.DrawDetectedGaussians(gaussimg, totalgausslist);	
			ImageJFunctions.show(gaussimg).setTitle("Iterated Result");
			
		 /*
		
		// Do Gaussian Fit
				final int ndims = biginputimg.numDimensions();
				 double[] final_param= new double[2*ndims+1];
				int firstlabel = lineparams.get(0).Label;
				int lastlabel = lineparams.get(lineparams.size() - 1).Label;
				final double[] psf = new double[ndims];
				for (int d = 0; d < ndims; ++d)
					psf[d] = 2;
		LengthDetection MTlength = new LengthDetection(inputimg, linepair.fst, lineparams);
		ArrayList<double[]> totalgausslist = new ArrayList<double[]>();
		for (int label = firstlabel; label <= lastlabel; ++label) {
			
			for (int listindex = 0; listindex< lineparams.size();++listindex){
				
			final_param = MTlength.Getfinalparam(psf, listindex, ndims, label);
				
			totalgausslist.add(final_param);
			System.out.println( " Label: "+ label + "   "+ " Listindex: " + listindex   + " Amplitude: " + final_param[0] + " " + "Mean X: "
					+ final_param[1] + " " + "Mean Y: " + final_param[2] + " " + "SigmaX: "
					+ 1.0 / Math.sqrt(final_param[3]) + " " + "SigmaY: "
					+ 1.0 / Math.sqrt(final_param[4]));
			
			}
			
		}
		PushCurves.DrawDetectedGaussians(gaussimg, totalgausslist);	
		ImageJFunctions.show(gaussimg).setTitle("Iterated Result");
		
/*
		// Do Gaussian Fit
		final int ndims = biginputimg.numDimensions();
		int firstlabel = simline.get(0).Label;
		int lastlabel = simline.get(simline.size() - 1).Label;
		double[] typical_sigma = new double[ndims];
		for (int d = 0; d < ndims; ++d)
			typical_sigma[d] = 100;

		ArrayList<double[]> totalgausslist = new ArrayList<double[]>();
		for (int label = firstlabel; label <= lastlabel; ++label) {

			ArrayList<double[]> gausslist = new ArrayList<double[]>();
			int setlength = simline.size()-1;

			gausslist = LengthDetection.makeBeads(inputimg, linepair.fst, simline, typical_sigma, ndims, label,
					setlength);
			for (int index = 0; index < gausslist.size(); ++index)
				totalgausslist.add(gausslist.get(index));
			
			PushCurves.DrawDetectedGaussians(gaussimg, gausslist);

			ImageJFunctions.show(gaussimg).setTitle("Iterated Result");

		}
		for (int index = 0; index < totalgausslist.size(); ++index)
			System.out.println("Amplitude: " + totalgausslist.get(index)[0] + " " + "Mean X: "
					+ totalgausslist.get(index)[1] + " " + "Mean Y: " + totalgausslist.get(index)[2] + " " + "SigmaX: "
					+ 1.0 / Math.sqrt(totalgausslist.get(index)[3]) + " " + "SigmaY: "
					+ 1.0 / Math.sqrt(totalgausslist.get(index)[4]));

	//	PushCurves.DrawDetectedGaussians(gaussimg, totalgausslist);

//		ImageJFunctions.show(gaussimg).setTitle("Iterated Result");
*/
	}
}
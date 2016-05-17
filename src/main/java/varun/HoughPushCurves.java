package varun;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;

import net.imglib2.algorithm.labeling.AllConnectedComponents;
import net.imglib2.algorithm.labeling.ConnectedComponents;
import net.imglib2.algorithm.localextrema.RefinedPeak;
import ij.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.IterableInterval;
import net.imglib2.Point;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.stats.Normalize;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.labeling.Labeling;
import net.imglib2.labeling.LabelingType;
import net.imglib2.labeling.NativeImgLabeling;
import net.imglib2.roi.labeling.ImgLabeling;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.IntegerType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import util.ImgLib2Util;
import varun.GetLocalmaxmin.IntensityType;
import varun.PerformWatershedding.InverseType;

public class HoughPushCurves {

	public static void Houghspace(RandomAccessibleInterval<FloatType> inputimage,
			RandomAccessibleInterval<FloatType> imgout, double[] min, double[] max, Float threshold) {

		int n = inputimage.numDimensions();

		final long[] position = new long[n];
		double Amplitude, Phase;

		final Cursor<FloatType> inputcursor = Views.iterable(inputimage).localizingCursor();

		// for every function (as defined by an individual pixel)
		while (inputcursor.hasNext()) {

			inputcursor.fwd();
			inputcursor.localize(position);
			if (inputcursor.get().get() > threshold) {
				Amplitude = Math.sqrt(Math.pow(position[0], 2) + Math.pow(position[1], 2));
				Phase = Math.toDegrees(Math.atan2(position[0], position[1]));

				// draw the function into the hough space

				PushCurves.DrawSine(imgout, min, max, Amplitude, Phase);
			}

		}
	}

	public static void main(String[] args) {

		RandomAccessibleInterval<FloatType> biginputimg = ImgLib2Util
				.openAs32Bit(new File("src/main/resources/2015-01-14_Porcine_Tubulin009-1.tiff"));

		new ImageJ();
		new Normalize();
		FloatType minval = new FloatType(0);
		FloatType maxval = new FloatType(1);
		Normalize.normalize(Views.iterable(biginputimg), minval, maxval);

		RandomAccessibleInterval<FloatType> inputimg = new ArrayImgFactory<FloatType>().create(biginputimg,
				new FloatType());
		RandomAccessibleInterval<FloatType> testinputimg = new ArrayImgFactory<FloatType>().create(biginputimg,
				new FloatType());
		final double[] sigma = { 1, 1 };
		
		  
		  ImageJFunctions.show(biginputimg).setTitle("Original image");
		  testinputimg = Kernels.NaiveEdge(biginputimg, sigma); 
		  ImageJFunctions.show(testinputimg).setTitle("Conditional Max image"); 
		  inputimg = Kernels.CannyEdge(biginputimg,new ArrayImgFactory<FloatType>(), sigma);
		  
		  // Automatic threshold determination for doing the Hough transform
		  final Float val = GlobalThresholding.AutomaticThresholding(inputimg);
		  
		  
		  ImageJFunctions.show(inputimg).setTitle("Input image"); 
		  int mintheta = 0;
		  
		  // Usually is 180 but to allow for detection of vertical lines,allowing // a few more degrees
		  int maxtheta = 200; 
		  double size = Math.sqrt((inputimg.dimension(0) * inputimg.dimension(0) +
		  inputimg.dimension(1) * inputimg.dimension(1))); int minRho = (int)
		  -Math.round(size); 
		  int maxRho = -minRho; 
		  // Set size of pixels in Hough space 
		  double thetaPerPixel = 0.1; 
		  double rhoPerPixel = 0.1;
		  double[] min = { mintheta, minRho }; 
		  double[] max = { maxtheta, maxRho };
		  
		  int pixelsTheta = (int) Math.round((maxtheta - mintheta) /
		  thetaPerPixel); int pixelsRho = (int) Math.round((maxRho - minRho) /
		  rhoPerPixel);
		  
		  double ratio = (max[0] - min[0]) / (max[1] - min[1]);
		  
		  // Size of Hough space 
		  FinalInterval interval = new FinalInterval(new long[] { pixelsTheta, (long) (pixelsRho * ratio) }); 
		  final Img<FloatType> houghimage = new
		  ArrayImgFactory<FloatType>().create(interval, new FloatType());
		  
		  ArrayList<RefinedPeak<Point>> SubpixelMinlist = new ArrayList<RefinedPeak<Point>>(inputimg.numDimensions()); 
		  final double[] sizes = new double[inputimg.numDimensions()]; 
		  for (int d = 0; d < houghimage.numDimensions(); ++d) 
			  sizes[d] = houghimage.dimension(d);
		  
		  // Do the Hough transform 
		  Houghspace(inputimg, houghimage, min, max,val); 
		  ImageJFunctions.show(houghimage).setTitle("Hough transform of input image"); 
		  final Float houghval = GlobalThresholding.AutomaticThresholding(houghimage);
		  // Get local Minima in scale space to get Max rho-theta points 
		  double minPeakValue = houghval; //0.09/(thetaPerPixel*rhoPerPixel); 
		  double smallsigma = 1;
		  double bigsigma = 1.1; 
		  SubpixelMinlist = GetLocalmaxmin.ScalespaceMinima(houghimage, interval, thetaPerPixel,
		  rhoPerPixel, minPeakValue, smallsigma, bigsigma);
		 


		// Reconstruct lines and overlay on the input image

		  OverlayLines.Overlay(biginputimg, SubpixelMinlist, sizes, min, max);

	}
}

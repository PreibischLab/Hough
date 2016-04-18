package varun;

import java.awt.Color;
import java.awt.Image;
import java.io.File;
import java.util.ArrayList;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import net.imglib2.algorithm.dog.DifferenceOfGaussian;
import net.imglib2.algorithm.dog.DogDetection;
import net.imglib2.algorithm.edge.Edgel;
import net.imglib2.algorithm.gauss3.Gauss3;
import net.imglib2.algorithm.localextrema.RefinedPeak;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Line;
import ij.gui.Overlay;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.IterableInterval;
import net.imglib2.Point;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealCursor;
import net.imglib2.RealPoint;
import net.imglib2.RealPointSampleList;
import net.imglib2.algorithm.stats.Normalize;
import net.imglib2.exception.IncompatibleTypeException;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import util.ImgLib2Util;

public class HoughPushCurves {

	public static void Houghspace(RandomAccessibleInterval<FloatType> inputimage,
			RandomAccessibleInterval<FloatType> imgout, double[] min, double[] max, FloatType threshold) {

		int n = inputimage.numDimensions();

		final long[] position = new long[n];
		double Amplitude, Phase;

		final Cursor<FloatType> inputcursor = Views.iterable(inputimage).localizingCursor();

		// for every function (as defined by an individual pixel)
		while (inputcursor.hasNext()) {

			inputcursor.fwd();
			inputcursor.localize(position);
			if (inputcursor.get().compareTo(threshold) > 0) {
				Amplitude = Math.sqrt(Math.pow(position[0], 2) + Math.pow(position[1], 2));
				Phase = Math.toDegrees(Math.atan2(position[0], position[1]));
				
				// draw the function into the hough space

				PushCurves.DrawSine(imgout, min, max, Amplitude, Phase);
			}

		}
	}

	public static void main(String[] args) {

		final RandomAccessibleInterval<FloatType> biginputimg = ImgLib2Util
				.openAs32Bit(new File("src/main/resources/Horizontal_line.tif"));
		new Normalize();
		 FloatType minval = new FloatType(0);
		 FloatType maxval = new FloatType(255);
		 Normalize.normalize(Views.iterable(biginputimg), minval, maxval);
		 
		
		// Preprocess the image and store it as new imputimg
		RandomAccessibleInterval<FloatType> inputimg = new ArrayImgFactory<FloatType>().create(biginputimg,
				new FloatType());
		
		 FloatType ThresholdValue = new FloatType(10);
		// Kernels.Edgedetector(biginputimg);
		// double filterradius = 1.0;
		// GetLocalmaxmin.MeanFilter(biginputimg,tmpinputimg,filterradius, ThresholdValue);
		// Thresholding the inputimage
		// GetLocalmaxmin.Thresholding(biginputimg, inputimg, ThresholdValue);
			inputimg = biginputimg;
		 
		 
		 Normalize.normalize(Views.iterable(inputimg), minval, maxval);
		 
		new ImageJ();
		ImageJFunctions.show(inputimg);
		// Usually is 180 but to allow for detection of vertical lines allowing for 20 more degrees
		int maxtheta = 200; 
		double size = Math
				.sqrt((inputimg.dimension(0) * inputimg.dimension(0) + inputimg.dimension(1) * inputimg.dimension(1)));
		int minRho = (int) -Math.round(size);
		int maxRho = -minRho;
		// Set size of pixels in Hough space
		double thetaPerPixel = 0.1;
		double rhoPerPixel = 0.1;

		

		double[] min = { 0, minRho };
		double[] max = { maxtheta, maxRho };

		int pixelsTheta = (int) Math.round((maxtheta) / thetaPerPixel);
		int pixelsRho = (int) Math.round((maxRho - minRho) / rhoPerPixel);

		double ratio = (max[0] - min[0]) / (max[1] - min[1]);

		// Size of Hough space
		FinalInterval interval = new FinalInterval(new long[] { pixelsTheta, (long) (pixelsRho * ratio) });
		final Img<FloatType> houghimage = new ArrayImgFactory<FloatType>().create(interval, new FloatType());
		ArrayList<RefinedPeak<Point>> SubpixelMinlist = new ArrayList<RefinedPeak<Point>>(inputimg.numDimensions());
		
		
		
		FloatType val = new FloatType(100);
		final double[] sizes = new double[inputimg.numDimensions()];

		for (int d = 0; d < houghimage.numDimensions(); ++d)
			sizes[d] = houghimage.dimension(d);

		// Do the Hough transform
		Houghspace(inputimg, houghimage, min, max, val);
		
		
		Kernels.BigButterflyKernel(houghimage);
		Normalize.normalize(houghimage, minval, maxval);
		
		ImageJFunctions.show(houghimage);
		
		// Get local Minima in sclae space
		SubpixelMinlist = GetLocalmaxmin.ScalespaceMinima(houghimage, interval, thetaPerPixel, rhoPerPixel, 
				0.5, 0.8, 1.2);
       // Reconstruct line and overlay on the input image
		GetLocalmaxmin.Overlaylines(inputimg,  SubpixelMinlist, sizes, 
				 min,  max);
		
	}

}

package houghandWatershed;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import com.sun.tools.javac.util.Pair;

import drawandOverlay.HoughPushCurves;
import drawandOverlay.OverlayLines;
import labeledObjects.Lineobjects;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.Interval;
import net.imglib2.KDTree;
import net.imglib2.Point;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealPoint;
import net.imglib2.RealPointSampleList;
import net.imglib2.algorithm.labeling.AllConnectedComponents;
import net.imglib2.algorithm.labeling.ConnectedComponents;
import net.imglib2.algorithm.labeling.Watershed;
import net.imglib2.algorithm.localextrema.RefinedPeak;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.labeling.DefaultROIStrategyFactory;
import net.imglib2.labeling.Labeling;
import net.imglib2.labeling.LabelingROIStrategy;
import net.imglib2.labeling.NativeImgLabeling;
import net.imglib2.neighborsearch.NearestNeighborSearchOnKDTree;
import net.imglib2.roi.labeling.ImgLabeling;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import preProcessing.GetLocalmaxmin;
import preProcessing.GlobalThresholding;

@SuppressWarnings("deprecation")
public class PerformWatershedding {
 

	public static enum InverseType {
		Straight, Inverse
	}

	public static Pair<Img<IntType>, ArrayList<Lineobjects>> DowatersheddingandHough(
			final RandomAccessibleInterval<FloatType> biginputimg,
			final RandomAccessibleInterval<FloatType> processedimg,
			final double minlength) 
	{

		// Prepare seed image for watershedding
		NativeImgLabeling<Integer, IntType> oldseedLabeling = new NativeImgLabeling<Integer, IntType>(
				new ArrayImgFactory<IntType>().create(biginputimg, new IntType()));

		oldseedLabeling = PrepareSeedImage(biginputimg);

		// Get maximum labels on the watershedded image

		final int Maxlabel = GetMaxlabelsseeded(oldseedLabeling.getStorageImg());

		ArrayList<RefinedPeak<Point>> ReducedMinlist = new ArrayList<RefinedPeak<Point>>(biginputimg.numDimensions());
		ArrayList<RefinedPeak<Point>> MainMinlist = new ArrayList<RefinedPeak<Point>>(biginputimg.numDimensions());

		// Perform the distance transform
		final Img<FloatType> distimg = new ArrayImgFactory<FloatType>().create(biginputimg, new FloatType());

		PerformWatershedding.DistanceTransformImage(biginputimg, distimg, InverseType.Straight);

		// Do watershedding on the distance transformed image

		NativeImgLabeling<Integer, IntType> outputLabeling = new NativeImgLabeling<Integer, IntType>(
				new ArrayImgFactory<IntType>().create(biginputimg, new IntType()));

		outputLabeling = GetlabeledImage(distimg, oldseedLabeling);
		
		final double[] sizes = new double[biginputimg.numDimensions()];

		// Automatic threshold determination for doing the Hough transform
		final Float val = GlobalThresholding.AutomaticThresholding(processedimg);
		
		ArrayList<Lineobjects> linelist = new ArrayList<Lineobjects>(biginputimg.numDimensions());
		
		
		for (int label = 1; label < Maxlabel-1; label++) {

			System.out.println("Label Number:" +label);
			

			RandomAccessibleInterval<FloatType> outimg = new ArrayImgFactory<FloatType>()
							.create(biginputimg, new FloatType());
		
			outimg = CurrentLabelImage(outputLabeling.getStorageImg(), processedimg,label);
			

			// Set size of pixels in Hough space
			int mintheta = 0;
			// Usually is 180 but to allow for detection of vertical
			// lines,allowing a few more degrees
			int maxtheta = 200;
			double size = Math
					.sqrt((outimg.dimension(0) * outimg.dimension(0) + outimg.dimension(1) * outimg.dimension(1)));
			int minRho = (int) -Math.round(size);
			int maxRho = -minRho;
			double thetaPerPixel = 1;
			double rhoPerPixel = 1;
			double[] min = { mintheta, minRho };
			double[] max = { maxtheta, maxRho };
			int pixelsTheta = (int) Math.round((maxtheta - mintheta) / thetaPerPixel);
			int pixelsRho = (int) Math.round((maxRho - minRho) / rhoPerPixel);

			double ratio = (max[0] - min[0]) / (max[1] - min[1]);
			FinalInterval interval = new FinalInterval(new long[] { pixelsTheta, (long) (pixelsRho * ratio) });
			final Img<FloatType> houghimage = new ArrayImgFactory<FloatType>().create(interval, new FloatType());
			long[] minCorner =  new long[biginputimg.numDimensions()];
			long[] maxCorner =  new long[biginputimg.numDimensions()];
			minCorner = GetMincorners(outputLabeling.getStorageImg(), label);
			maxCorner = GetMaxcorners(outputLabeling.getStorageImg(), label);
			
			
			FinalInterval intervalsmall = new FinalInterval( minCorner, maxCorner );
			
			RandomAccessibleInterval<FloatType> outimgview = Views.interval(outimg, intervalsmall);
			HoughPushCurves.Houghspace(outimgview, houghimage, min, max, val);


			for (int d = 0; d < houghimage.numDimensions(); ++d)
				sizes[d] = houghimage.dimension(d);
			ArrayList<RefinedPeak<Point>> SubpixelMinlist = new ArrayList<RefinedPeak<Point>>(
					biginputimg.numDimensions());
			SubpixelMinlist = GetLocalmaxmin.HoughspaceMaxima(houghimage, interval, sizes, thetaPerPixel, rhoPerPixel);

			ReducedMinlist = OverlayLines.ReducedList(outimg, SubpixelMinlist, sizes, min, max, minlength);
			
			double[] points = new double[biginputimg.numDimensions()];
			
			
				for (int index = 0; index < ReducedMinlist.size(); ++index)
				MainMinlist.add(ReducedMinlist.get(index));
			
			points = OverlayLines.GetRhoTheta( ReducedMinlist, sizes, min, max, minlength);
			 
			// This object has rho, theta, min and max dimensions of the watershedded image along x 	
			final Lineobjects line = new Lineobjects(label, points[1], points[0], minCorner, maxCorner);

			linelist.add(line);
			
		}
		Pair<Img<IntType>, ArrayList<Lineobjects>> linepair = new Pair<Img<IntType>, ArrayList<Lineobjects>>(outputLabeling.getStorageImg(), linelist);
		

		return linepair;
	}

	public static RandomAccessibleInterval<IntType> Labelobjects(
			final RandomAccessibleInterval<FloatType> biginputimg) 
	{

		// Prepare seed image for watershedding
		NativeImgLabeling<Integer, IntType> oldseedLabeling = new NativeImgLabeling<Integer, IntType>(
				new ArrayImgFactory<IntType>().create(biginputimg, new IntType()));

		oldseedLabeling = PrepareSeedImage(biginputimg);

		
		
		return oldseedLabeling.getStorageImg();
	}
	
	
	public static RandomAccessibleInterval<IntType> Dowatersheddingonly(
			final RandomAccessibleInterval<FloatType> biginputimg) 
	{

		// Prepare seed image for watershedding
		NativeImgLabeling<Integer, IntType> oldseedLabeling = new NativeImgLabeling<Integer, IntType>(
				new ArrayImgFactory<IntType>().create(biginputimg, new IntType()));

		oldseedLabeling = PrepareSeedImage(biginputimg);

		// Perform the distance transform
		final Img<FloatType> distimg = new ArrayImgFactory<FloatType>().create(biginputimg, new FloatType());

		PerformWatershedding.DistanceTransformImage(biginputimg, distimg, InverseType.Straight);

		// Do watershedding on the distance transformed image

		NativeImgLabeling<Integer, IntType> outputLabeling = new NativeImgLabeling<Integer, IntType>(
				new ArrayImgFactory<IntType>().create(biginputimg, new IntType()));

		outputLabeling = GetlabeledImage(distimg, oldseedLabeling);
		
		return outputLabeling.getStorageImg();
	}
	
	
	public static Lineobjects Getlabelobject(
			final RandomAccessibleInterval<FloatType> biginputimg,
			final RandomAccessibleInterval<FloatType> processedimg,
			final int currentLabel) 
	{

		// Prepare seed image for watershedding
		NativeImgLabeling<Integer, IntType> oldseedLabeling = new NativeImgLabeling<Integer, IntType>(
				new ArrayImgFactory<IntType>().create(biginputimg, new IntType()));

		oldseedLabeling = PrepareSeedImage(biginputimg);

		// Get maximum labels on the watershedded image

		ArrayList<RefinedPeak<Point>> ReducedMinlist = new ArrayList<RefinedPeak<Point>>(biginputimg.numDimensions());

		// Perform the distance transform
		final Img<FloatType> distimg = new ArrayImgFactory<FloatType>().create(biginputimg, new FloatType());

		PerformWatershedding.DistanceTransformImage(biginputimg, distimg, InverseType.Straight);

		// Do watershedding on the distance transformed image

		NativeImgLabeling<Integer, IntType> outputLabeling = new NativeImgLabeling<Integer, IntType>(
				new ArrayImgFactory<IntType>().create(biginputimg, new IntType()));

		outputLabeling = GetlabeledImage(distimg, oldseedLabeling);
		
		final double[] sizes = new double[biginputimg.numDimensions()];

		// Automatic threshold determination for doing the Hough transform
		final Float val = GlobalThresholding.AutomaticThresholding(processedimg);
		
		ArrayList<Lineobjects> linelist = new ArrayList<Lineobjects>(biginputimg.numDimensions());
		
		// Declare minimum length of the line(in pixels) to be detected
		double minlength = 0;
		int label = currentLabel;

			System.out.println("Label Number:" +label);
			

			RandomAccessibleInterval<FloatType> outimg = new ArrayImgFactory<FloatType>()
							.create(biginputimg, new FloatType());
		
			outimg = CurrentLabelImage(outputLabeling.getStorageImg(), processedimg,label);
			

			// Set size of pixels in Hough space
			int mintheta = 0;
			// Usually is 180 but to allow for detection of vertical
			// lines,allowing a few more degrees
			int maxtheta = 200;
			double size = Math
					.sqrt((biginputimg.dimension(0) * biginputimg.dimension(0) + biginputimg.dimension(1) * biginputimg.dimension(1)));
			int minRho = (int) -Math.round(size);
			int maxRho = -minRho;
			double thetaPerPixel = 0.4;
			double rhoPerPixel = 0.4;
			double[] min = { mintheta, minRho };
			double[] max = { maxtheta, maxRho };
			int pixelsTheta = (int) Math.round((maxtheta - mintheta) / thetaPerPixel);
			int pixelsRho = (int) Math.round((maxRho - minRho) / rhoPerPixel);

			double ratio = (max[0] - min[0]) / (max[1] - min[1]);
			FinalInterval interval = new FinalInterval(new long[] { pixelsTheta, (long) (pixelsRho * ratio) });
			final Img<FloatType> houghimage = new ArrayImgFactory<FloatType>().create(interval, new FloatType());
			
			long[] minCorner =  new long[biginputimg.numDimensions()];
			long[] maxCorner =  new long[biginputimg.numDimensions()];
			minCorner = GetMincorners(outputLabeling.getStorageImg(), label);
			maxCorner = GetMaxcorners(outputLabeling.getStorageImg(), label);
			
			
			FinalInterval intervalsmall = new FinalInterval( minCorner, maxCorner );
			
			RandomAccessibleInterval<FloatType> outimgview = Views.interval(outimg, intervalsmall);
			HoughPushCurves.Houghspace(outimgview, houghimage, min, max, val);


			for (int d = 0; d < houghimage.numDimensions(); ++d)
				sizes[d] = houghimage.dimension(d);
			ArrayList<RefinedPeak<Point>> SubpixelMinlist = new ArrayList<RefinedPeak<Point>>(
					biginputimg.numDimensions());
			SubpixelMinlist = GetLocalmaxmin.HoughspaceMaxima(houghimage, interval, sizes, thetaPerPixel, rhoPerPixel);

			ReducedMinlist = OverlayLines.ReducedList(outimg, SubpixelMinlist, sizes, min, max, minlength);
			
			double[] points = new double[biginputimg.numDimensions()];
			
			
			points = OverlayLines.GetRhoTheta( ReducedMinlist, sizes, min, max, minlength);
			 
			// This object has rho, theta, min and max dimensions of the watershedded image along x 	
			final Lineobjects line = new Lineobjects(label, points[1], points[0], minCorner, maxCorner);

			linelist.add(line);
			
		
		
		

		return line;
	}

	public static void DistanceTransformImage(
			RandomAccessibleInterval<FloatType> inputimg,
			RandomAccessibleInterval<FloatType> outimg, 
			final InverseType invtype) 
	{
		int n = inputimg.numDimensions();

		final Img<BitType> bitimg = new ArrayImgFactory<BitType>().create(inputimg, new BitType());
		// make an empty list
		final RealPointSampleList<BitType> list = new RealPointSampleList<BitType>(n);

		final Float threshold = GlobalThresholding.AutomaticThresholding(inputimg);
		GetLocalmaxmin.ThresholdingBit(inputimg, bitimg, threshold);

		// cursor on the binary image
		final Cursor<BitType> cursor = bitimg.localizingCursor();

		// for every pixel that is 1, make a new RealPoint at that location
		while (cursor.hasNext())
			if (cursor.next().getInteger() == 1)
				list.add(new RealPoint(cursor), cursor.get());

		// build the KD-Tree from the list of points that == 1
		final KDTree<BitType> tree = new KDTree<BitType>(list);

		// Instantiate a nearest neighbor search on the tree (does not modifiy
		// the tree, just uses it)
		final NearestNeighborSearchOnKDTree<BitType> search = new NearestNeighborSearchOnKDTree<BitType>(tree);

		// randomaccess on the output
		final RandomAccess<FloatType> ranac = outimg.randomAccess();

		// reset cursor for the input (or make a new one)
		cursor.reset();

		// for every pixel of the binary image
		while (cursor.hasNext()) {
			cursor.fwd();

			// set the randomaccess to the same location
			ranac.setPosition(cursor);

			// if value == 0, look for the nearest 1-valued pixel
			if (cursor.get().getInteger() == 0) {
				// search the nearest 1 to the location of the cursor (the
				// current 0)
				search.search(cursor);

				// get the distance (the previous call could return that, this
				// for generality that it is two calls)
				switch (invtype) {

				case Straight:
					ranac.get().setReal(search.getDistance());
					break;
				case Inverse:
					ranac.get().setReal(-search.getDistance());
					break;

				}

			} else {
				// if value == 1, no need to search
				ranac.get().setZero();
			}
		}

	}

	public static NativeImgLabeling<Integer, IntType> PrepareSeedImage(
			RandomAccessibleInterval<FloatType> inputimg) 
	{

		// Preparing the seed image
		RandomAccessibleInterval<BitType> maximgBit = new ArrayImgFactory<BitType>().create(inputimg, new BitType());
		final Float threshold = GlobalThresholding.AutomaticThresholding(inputimg);
		GetLocalmaxmin.ThresholdingBit(inputimg, maximgBit, threshold);

		// New Labeling type
		final ImgLabeling<Integer, IntType> seedLabeling = new ImgLabeling<Integer, IntType>(
				new ArrayImgFactory<IntType>().create(inputimg, new IntType()));

		// Old Labeling type
		final NativeImgLabeling<Integer, IntType> oldseedLabeling = new NativeImgLabeling<Integer, IntType>(
				new ArrayImgFactory<IntType>().create(inputimg, new IntType()));

		// The label generator for both new and old type
		final Iterator<Integer> labelGenerator = AllConnectedComponents.getIntegerNames(0);

		// Getting unique labelled image (new version)
		ConnectedComponents.labelAllConnectedComponents(maximgBit, seedLabeling, labelGenerator,
				ConnectedComponents.StructuringElement.EIGHT_CONNECTED);

		// Getting unique labelled image (old version)
		AllConnectedComponents.labelAllConnectedComponents(oldseedLabeling, maximgBit, labelGenerator,
				AllConnectedComponents.getStructuringElement(inputimg.numDimensions()));

		return oldseedLabeling;
	}

      public static long[] GetMaxcorners(
    		  Img<IntType> inputimg, 
    		  int label){
		

		Cursor<IntType> intCursor = inputimg.localizingCursor();
		int n = inputimg.numDimensions();
		long[] maxVal = { Long.MIN_VALUE, Long.MIN_VALUE };

		while (intCursor.hasNext()) {
			intCursor.fwd();
			int i = intCursor.get().get();
			if (i == label) {

				for (int d = 0; d < n; ++d) {
					
					final long p = intCursor.getLongPosition( d );
					if ( p > maxVal[ d ] )
						maxVal[ d ] = p;
					
				}

			}
		}
         
		
		return maxVal;
		
	}
	public static long[] GetMincorners(
			Img<IntType> inputimg, 
			int label){
		

		Cursor<IntType> intCursor = inputimg.localizingCursor();
		int n = inputimg.numDimensions();
		long[] minVal = { Long.MAX_VALUE, Long.MAX_VALUE };

		while (intCursor.hasNext()) {
			intCursor.fwd();
			int i = intCursor.get().get();
			if (i == label) {

				for (int d = 0; d < n; ++d) {
					
					final long p = intCursor.getLongPosition( d );
					if ( p < minVal[ d ] )
						minVal[ d ] = p;
				}

			}
		}
		
		
		return minVal;
		
	}
	
	public static Pair<long[], long[]> GetBoundingbox(
			Img<IntType> inputimg, 
			int label) {

		Cursor<IntType> intCursor = inputimg.localizingCursor();
		int n = inputimg.numDimensions();
		long[] position = new long[n];
		long[] minVal = { Long.MAX_VALUE, Long.MAX_VALUE };
		long[] maxVal = { Long.MIN_VALUE, Long.MIN_VALUE };

		while (intCursor.hasNext()) {
			intCursor.fwd();
			int i = intCursor.get().get();
			if (i == label) {

				intCursor.localize(position);
				for (int d = 0; d < n; ++d) {
					if (position[d] < minVal[d]) {
						minVal[d] = position[d];
					}
					if (position[d] > maxVal[d]) {
						maxVal[d] = position[d];
					}

				}

			}
		}

		Pair<long[], long[]> boundingBox = new Pair<long[], long[]>(minVal, maxVal);
		return boundingBox;
	}

	

	public static int GetMaxlabelsseeded(
			RandomAccessibleInterval<IntType> intimg) 
	{
		
		// To get maximum Labels on the image
		Cursor<IntType> intCursor = Views.iterable(intimg).cursor();
		int currentLabel = 1;
		boolean anythingFound = true;
		while (anythingFound) {
			anythingFound = false;
			intCursor.reset();
			while (intCursor.hasNext()) {
				intCursor.fwd();
				int i = intCursor.get().get();
				if (i == currentLabel) {

					anythingFound = true;

				}
			}
			currentLabel++;
		}

		return currentLabel;

	}
	public static NativeImgLabeling<Integer, IntType> GetlabeledImage(
			RandomAccessibleInterval<FloatType> inputimg,
			NativeImgLabeling<Integer, IntType> seedLabeling) 
	{

		int n = inputimg.numDimensions();
		long[] dimensions = new long[n];

		for (int d = 0; d < n; ++d)
			dimensions[d] = inputimg.dimension(d);
		final NativeImgLabeling<Integer, IntType> outputLabeling = new NativeImgLabeling<Integer, IntType>(
				new ArrayImgFactory<IntType>().create(inputimg, new IntType()));

		final Watershed<FloatType, Integer> watershed = new Watershed<FloatType, Integer>();

		watershed.setSeeds(seedLabeling);
		watershed.setIntensityImage(inputimg);
		watershed.setStructuringElement(AllConnectedComponents.getStructuringElement(2));
		watershed.setOutputLabeling(outputLabeling);
		watershed.process();
		DefaultROIStrategyFactory<Integer> deffactory = new DefaultROIStrategyFactory<Integer>();
		LabelingROIStrategy<Integer, Labeling<Integer>> factory = deffactory
				.createLabelingROIStrategy(watershed.getResult());
		outputLabeling.setLabelingCursorStrategy(factory);

		return outputLabeling;

	}

	public static RandomAccessibleInterval<FloatType> CurrentLabelImage(
			Img<IntType> Intimg,
			RandomAccessibleInterval<FloatType> originalimg, 
			int currentLabel) 
	{

		RandomAccess<FloatType> inputRA = originalimg.randomAccess();

		Cursor<IntType> intCursor = Intimg.cursor();
		
 		RandomAccessibleInterval<FloatType> outimg = new ArrayImgFactory<FloatType>()
				.create(originalimg, new FloatType());
		RandomAccess<FloatType> imageRA = outimg.randomAccess();

		// Go through the whole image and add every pixel, that belongs to
		// the currently processed label

		while (intCursor.hasNext()) {
			intCursor.fwd();
			inputRA.setPosition(intCursor);
			imageRA.setPosition(inputRA);
			int i = intCursor.get().get();
			if (i == currentLabel) {
				
				imageRA.get().set(inputRA.get());
				
			}

		}

		return outimg;

	}

}
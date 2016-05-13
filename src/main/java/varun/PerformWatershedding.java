package varun;

import java.awt.List;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Set;

import com.sun.tools.javac.util.Pair;

import ij.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
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
import net.imglib2.converter.Converter;
import net.imglib2.converter.RealARGBConverter;
import net.imglib2.converter.read.ConvertedRandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.NativeImg;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.basictypeaccess.IntAccess;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.labeling.Labeling;
import net.imglib2.labeling.LabelingType;
import net.imglib2.labeling.NativeImgLabeling;
import net.imglib2.multithreading.SimpleMultiThreading;
import net.imglib2.neighborsearch.NearestNeighborSearchOnKDTree;
import net.imglib2.roi.labeling.ImgLabeling;
import net.imglib2.roi.labeling.LabelRegions;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.IntegerType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.integer.ShortType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Fraction;
import net.imglib2.view.Views;
import varun.GetLocalmaxmin.IntensityType;
import varun.PerformWatershedding.InverseType;

public class PerformWatershedding {

	public static enum InverseType {
		Straight, Inverse
	}

	public static void DistanceTransformImage(RandomAccessibleInterval<FloatType> inputimg,
			RandomAccessibleInterval<FloatType> outimg, final InverseType invtype) {
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

	@SuppressWarnings("deprecation")
	public static void Dowatershedding(RandomAccessibleInterval<FloatType> inputimg) {
		final Img<FloatType> distimg = new ArrayImgFactory<FloatType>().create(inputimg, new FloatType());

		// Distance transfrom the original image
		PerformWatershedding.DistanceTransformImage(inputimg, distimg, InverseType.Straight);

		// Preparing the seed image for wateshedding
		RandomAccessibleInterval<BitType> maximgBit = new ArrayImgFactory<BitType>().create(inputimg, new BitType());
		final Float threshold = GlobalThresholding.AutomaticThresholding(inputimg);
		GetLocalmaxmin.ThresholdingBit(inputimg, maximgBit, threshold);

		// For the New Labeling type
		final ImgLabeling<Integer, IntType> seedLabeling = new ImgLabeling<Integer, IntType>(
				new ArrayImgFactory<IntType>().create(maximgBit, new IntType()));
		// For the Old Labeling type
		final NativeImgLabeling<Integer, IntType> oldseedLabeling = new NativeImgLabeling<Integer, IntType>(
				new ArrayImgFactory<IntType>().create(maximgBit, new IntType()));
		// The label generator for both new and old type
		final Iterator<Integer> labelGenerator = AllConnectedComponents.getIntegerNames(0);
		// Getting unique labelled image (new version)
		ConnectedComponents.labelAllConnectedComponents(maximgBit, seedLabeling, labelGenerator,
				ConnectedComponents.StructuringElement.EIGHT_CONNECTED);
		// Getting unique labelled image (old version)
		AllConnectedComponents.labelAllConnectedComponents(oldseedLabeling, maximgBit, labelGenerator,
				AllConnectedComponents.getStructuringElement(inputimg.numDimensions()));

		// ImageJFunctions.show(maximgBit).setTitle("Bit Seed Image for
		// watershed");
		// ImageJFunctions.show(seedLabeling.getIndexImg()).setTitle("New-Method");
		// ImageJFunctions.show(oldseedLabeling.getStorageImg()).setTitle("Old-Method");
		OldWatersherImage(distimg, oldseedLabeling, inputimg);
	}

	@SuppressWarnings("deprecation")
	public static void DowatersheddingandHough(RandomAccessibleInterval<FloatType> inputimg) {

		// Do the mean filtering and edge detection to make a pre-processed
		// image
		RandomAccessibleInterval<FloatType> meanfilterimg = new ArrayImgFactory<FloatType>().create(inputimg,
				new FloatType());
		

		meanfilterimg = Kernels.Meanfilterandsupress(inputimg, 1.0);

		// We can choose to do the Hough transform on the edge image or the
		// original image
		
		ImageJFunctions.show(meanfilterimg);

		final Img<FloatType> distimg = new ArrayImgFactory<FloatType>().create(inputimg, new FloatType());

		PerformWatershedding.DistanceTransformImage(inputimg, distimg, InverseType.Straight);

		// Preparing the seed image
		RandomAccessibleInterval<BitType> maximgBit = new ArrayImgFactory<BitType>().create(inputimg, new BitType());
		final Float threshold = GlobalThresholding.AutomaticThresholding(inputimg);
		GetLocalmaxmin.ThresholdingBit(inputimg, maximgBit, threshold);

		// New Labeling type
		final ImgLabeling<Integer, IntType> seedLabeling = new ImgLabeling<Integer, IntType>(
				new ArrayImgFactory<IntType>().create(maximgBit, new IntType()));
		// Old Labeling type
		final NativeImgLabeling<Integer, IntType> oldseedLabeling = new NativeImgLabeling<Integer, IntType>(
				new ArrayImgFactory<IntType>().create(maximgBit, new IntType()));
		// The label generator for both new and old type
		final Iterator<Integer> labelGenerator = AllConnectedComponents.getIntegerNames(0);
		// Getting unique labelled image (new version)
		ConnectedComponents.labelAllConnectedComponents(maximgBit, seedLabeling, labelGenerator,
				ConnectedComponents.StructuringElement.EIGHT_CONNECTED);
		// Getting unique labelled image (old version)
		AllConnectedComponents.labelAllConnectedComponents(oldseedLabeling, maximgBit, labelGenerator,
				AllConnectedComponents.getStructuringElement(inputimg.numDimensions()));

		RandomAccessibleInterval<FloatType> outimg = new ArrayImgFactory<FloatType>().create(distimg, new FloatType());

		Cursor<IntType> intCursor = oldseedLabeling.getStorageImg().cursor();
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

		ArrayList<RefinedPeak<Point>> MainMinlist = new ArrayList<RefinedPeak<Point>>(inputimg.numDimensions());
		final double[] sizes = new double[inputimg.numDimensions()];
		int mintheta = 0;
		// Usually is 180 but to allow for detection of vertical lines,allowing
		// // a few more degrees
		int maxtheta = 200;
		double size = Math
				.sqrt((outimg.dimension(0) * outimg.dimension(0) + outimg.dimension(1) * outimg.dimension(1)));
		int minRho = (int) -Math.round(size);
		int maxRho = -minRho;
		// Set size of pixels in Hough space
		double thetaPerPixel = 1;
		double rhoPerPixel = 1;
		double[] min = { mintheta, minRho };
		double[] max = { maxtheta, maxRho };
		NativeImgLabeling<Integer, IntType> outputLabeling = new NativeImgLabeling<Integer, IntType>(
				new ArrayImgFactory<IntType>().create(inputimg, new IntType()));
		outputLabeling = GetlabeledImage(distimg, oldseedLabeling);
		// int testlabel = 10;
		for (int label = 1; label < currentLabel - 1; label++) {
			// for (int label = 1; label<testlabel; label++){
			// choose the inputimage to be original or after passing through
			// edge detector
			outimg = CurrentLabelImage(outputLabeling, meanfilterimg, label);

			// Automatic threshold determination for doing the Hough transform
			final Float val = GlobalThresholding.AutomaticThresholding(outimg);

			int pixelsTheta = (int) Math.round((maxtheta - mintheta) / thetaPerPixel);
			int pixelsRho = (int) Math.round((maxRho - minRho) / rhoPerPixel);

			double ratio = (max[0] - min[0]) / (max[1] - min[1]);
			FinalInterval interval = new FinalInterval(new long[] { pixelsTheta, (long) (pixelsRho * ratio) });
			final Img<FloatType> houghimage = new ArrayImgFactory<FloatType>().create(interval, new FloatType());
			// Do the Hough transform

			HoughPushCurves.Houghspace(outimg, houghimage, min, max, val);
			// ImageJFunctions.show(houghimage);
			final Float houghval = GlobalThresholding.AutomaticThresholding(houghimage);
			ArrayList<RefinedPeak<Point>> SubpixelMinlist = new ArrayList<RefinedPeak<Point>>(outimg.numDimensions());

			for (int d = 0; d < houghimage.numDimensions(); ++d)
				sizes[d] = houghimage.dimension(d);
			// Get local Minima in scale space to get Max rho-theta points
			double minPeakValue = 1.2*houghval; // 0.09/(thetaPerPixel*rhoPerPixel);
			double smallsigma = 1;
			double bigsigma = 1.1;
			SubpixelMinlist = GetLocalmaxmin.ScalespaceMinima(houghimage, interval, thetaPerPixel, rhoPerPixel,
					minPeakValue, smallsigma, bigsigma);

			for (int index = 0; index < SubpixelMinlist.size(); ++index) {

				MainMinlist.add(SubpixelMinlist.get(index));

			}
			
		}

		// Reconstruct lines and overlay on the input image

		OverlayLines.Overlay(inputimg, MainMinlist, sizes, min, max);
		OverlayLines.Overlay(meanfilterimg, MainMinlist, sizes, min, max);
		
	}

	@SuppressWarnings("deprecation")
	public static NativeImgLabeling<Integer, IntType> GetlabeledImage(RandomAccessibleInterval<FloatType> inputimg,
			NativeImgLabeling<Integer, IntType> seedLabeling) {

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
		watershed.getResult();
		watershed.process();

		return outputLabeling;

	}

	@SuppressWarnings("deprecation")
	public static RandomAccessibleInterval<FloatType> CurrentLabelImage(
			NativeImgLabeling<Integer, IntType> outputLabeling, RandomAccessibleInterval<FloatType> originalimg,
			int currentLabel) {
		assert currentLabel > 0;
		RandomAccess<FloatType> inputRA = originalimg.randomAccess();
		Cursor<IntType> intCursor = outputLabeling.getStorageImg().cursor();
		RandomAccessibleInterval<FloatType> outimg = new ArrayImgFactory<FloatType>().create(originalimg,
				new FloatType());
		RandomAccess<FloatType> imageRA = outimg.randomAccess();

		// Go through the whole image and add every pixel, that belongs to
		// the currently processed label

		while (intCursor.hasNext()) {
			intCursor.fwd();

			imageRA.setPosition(intCursor);
			inputRA.setPosition(intCursor);
			int i = intCursor.get().get();
			if (i == currentLabel) {
				imageRA.get().set(inputRA.get());

			}

		}
		return outimg;

	}

	@SuppressWarnings("deprecation")
	public static void OldWatersherImage(RandomAccessibleInterval<FloatType> inputimg,
			NativeImgLabeling<Integer, IntType> seedLabeling, RandomAccessibleInterval<FloatType> originalimg) {

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
		watershed.getResult();
		watershed.process();
		ImageJFunctions.show(outputLabeling.getStorageImg()).setTitle("labeling storage image");

		RandomAccess<FloatType> inputRA = originalimg.randomAccess();
		Cursor<IntType> intCursor = outputLabeling.getStorageImg().cursor();

		int currentLabel = 1;
		boolean anythingFound = true;
		while (anythingFound) {
			anythingFound = false;
			intCursor.reset();
			RandomAccessibleInterval<FloatType> outimg = new ArrayImgFactory<FloatType>().create(inputimg,
					new FloatType());
			RandomAccess<FloatType> imageRA = outimg.randomAccess();
			// Go through the whole image and add every pixel, that belongs to
			// the currently processed label
			int count = 0;
			while (intCursor.hasNext()) {
				intCursor.fwd();

				imageRA.setPosition(intCursor);
				inputRA.setPosition(intCursor);
				int i = intCursor.get().get();
				if (i == currentLabel) {
					imageRA.get().set(inputRA.get());
					anythingFound = true;
					count++;

				}

			}
			System.out.println("Number of input pixels in label " + currentLabel + ": " + count);
			currentLabel++;
			ImageJFunctions.show(outimg).setTitle("Watershed Images");

		}
	}

	/*
	 * public static void WatershedImage(RandomAccessibleInterval<FloatType>
	 * inputimg,RandomAccessibleInterval<FloatType> seedimg){
	 * 
	 * int n = inputimg.numDimensions(); long[] dimensions = new long[n]; for
	 * (int d = 0; d < n; ++d ) dimensions[d] = inputimg.dimension(d);
	 * 
	 * RandomAccessibleInterval<FloatType> outputimg = new
	 * ArrayImgFactory<FloatType>().create(inputimg, new FloatType());
	 * ImgLabeling<Double, ShortType> seeds = new ImgLabeling<Double,
	 * ShortType>(new ArrayImgFactory<ShortType>().create(seedimg, new
	 * ShortType())); ImgLabeling< Double, ShortType > seedLabeling = new
	 * ImgLabeling< Double, ShortType >( new ArrayImgFactory< ShortType
	 * >().create( dimensions, new ShortType() ) );
	 * RandomAccessibleInterval<LabelingType<Float>> watershedResult = new
	 * ImgLabeling<Float, ShortType>(new
	 * ArrayImgFactory<ShortType>().create(inputimg, new ShortType()));
	 * 
	 * final Cursor<LabelingType<Double>> seedcursor = seeds.localizingCursor();
	 * 
	 * 
	 * }
	 * 
	 * 
	 * Converter<LabelingType<Integer>, FloatType> converter = new
	 * Converter<LabelingType<Integer>, FloatType>() {
	 * 
	 * @Override
	 * 
	 * public void convert(LabelingType<Integer> input, FloatType output) {
	 * 
	 * java.util.List<Integer> labelNames = input.getMapping().getLabels();
	 * 
	 * 
	 * }
	 * 
	 * };
	 * 
	 * ConvertedRandomAccessibleInterval<LabelingType<Integer>, FloatType>
	 * converted = new ConvertedRandomAccessibleInterval<LabelingType<Integer>,
	 * FloatType>( outputLabeling.copy(), converter, new FloatType());
	 * 
	 * ImageJFunctions.show(converted);
	 */

}

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

	@SuppressWarnings("deprecation")
	public static void DowatersheddingandHough(RandomAccessibleInterval<FloatType> inputimg, double[] min, double[] max,
			int pixelsTheta, int pixelsRho, double thetaPerPixel, double rhoPerPixel) {

		// Prepare seed image for watershedding
		NativeImgLabeling<Integer, IntType> oldseedLabeling = new NativeImgLabeling<Integer, IntType>(
				new ArrayImgFactory<IntType>().create(inputimg, new IntType()));

		oldseedLabeling = PrepareSeedImage(inputimg);

		RandomAccessibleInterval<FloatType> outimg = new ArrayImgFactory<FloatType>().create(inputimg, new FloatType());

		// Get maximum labels on the watershedded image

		final int currentLabel = GetMaxlabels(oldseedLabeling);

		ArrayList<RefinedPeak<Point>> ReducedMinlist = new ArrayList<RefinedPeak<Point>>(inputimg.numDimensions());
		ArrayList<RefinedPeak<Point>> MainMinlist = new ArrayList<RefinedPeak<Point>>(inputimg.numDimensions());

		// Perform the distance transform
		final Img<FloatType> distimg = new ArrayImgFactory<FloatType>().create(inputimg, new FloatType());

		PerformWatershedding.DistanceTransformImage(inputimg, distimg, InverseType.Straight);

		// Do watershedding on the distance transformed image

		NativeImgLabeling<Integer, IntType> outputLabeling = new NativeImgLabeling<Integer, IntType>(
				new ArrayImgFactory<IntType>().create(inputimg, new IntType()));

		outputLabeling = GetlabeledImage(distimg, oldseedLabeling);
		final double[] sizes = new double[inputimg.numDimensions()];
		for (int label = 1; label < currentLabel - 1; label++) {

			outimg = CurrentLabelImage(outputLabeling, inputimg, label);

			// Do the Hough transform on the outimg

			// Automatic threshold determination for doing the Hough transform
			final Float val = GlobalThresholding.AutomaticThresholding(outimg);

			double ratio = (max[0] - min[0]) / (max[1] - min[1]);
			FinalInterval interval = new FinalInterval(new long[] { pixelsTheta, (long) (pixelsRho * ratio) });
			final Img<FloatType> houghimage = new ArrayImgFactory<FloatType>().create(interval, new FloatType());

			HoughPushCurves.Houghspace(outimg, houghimage, min, max, val);
			
			for (int d = 0; d < houghimage.numDimensions(); ++d)
				sizes[d] = houghimage.dimension(d);
			ArrayList<RefinedPeak<Point>> SubpixelMinlist = new ArrayList<RefinedPeak<Point>>(inputimg.numDimensions());
			SubpixelMinlist = GetLocalmaxmin.HoughspaceMaxima(houghimage, interval, sizes, thetaPerPixel, rhoPerPixel);

			

			ReducedMinlist = OverlayLines.ReducedList(outimg, SubpixelMinlist, sizes, min, max);
			MainMinlist.add(ReducedMinlist.get(0));
			OverlayLines.Overlay(outimg, ReducedMinlist, sizes, min, max);
		}

		// Reconstruct lines and overlay on the input image
		OverlayLines.Overlay(inputimg, MainMinlist, sizes, min, max);

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
	public static NativeImgLabeling<Integer, IntType> PrepareSeedImage(RandomAccessibleInterval<FloatType> inputimg) {

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

	@SuppressWarnings("deprecation")
	public static int GetMaxlabels(NativeImgLabeling<Integer, IntType> oldseedLabeling) {

		// To get maximum Labels on the image
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

		return currentLabel;

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

}

package houghandWatershed;

import java.util.Iterator;

import com.sun.tools.javac.util.Pair;

import houghandWatershed.PerformWatershedding.InverseType;
import net.imglib2.Cursor;
import net.imglib2.KDTree;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealPoint;
import net.imglib2.RealPointSampleList;
import net.imglib2.algorithm.BenchmarkAlgorithm;
import net.imglib2.algorithm.OutputAlgorithm;
import net.imglib2.algorithm.labeling.AllConnectedComponents;
import net.imglib2.algorithm.labeling.ConnectedComponents;
import net.imglib2.algorithm.labeling.Watershed;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.labeling.DefaultROIStrategyFactory;
import net.imglib2.labeling.Labeling;
import net.imglib2.labeling.LabelingROIStrategy;
import net.imglib2.labeling.NativeImgLabeling;
import net.imglib2.neighborsearch.NearestNeighborSearchOnKDTree;
import net.imglib2.roi.labeling.ImgLabeling;
import net.imglib2.type.NativeType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Util;
import net.imglib2.view.Views;

/**
 * Hough transform of images that operates on 2D images.
 * <p>
 * For 3D images, the Hough transform is done only in 2D XY slices. 
 * 
 * @author Varun Kapoor - 2016
 *
 * @param <T>
 *            the type of the source image.
 */

public class HoughTransform2D< T extends RealType< T > & NativeType< T >> extends BenchmarkAlgorithm implements OutputAlgorithm< RandomAccessibleInterval< T >> {

	

	private final RandomAccessibleInterval<T> source;
	private final RandomAccessibleInterval<BitType> bitimg;
	private final double minlength;
	private static final String BASE_ERROR_MSG = "[HoughTransform2D] ";
	
	/**
	 * Instantiate a new Hough Transform object that does Hough
	 * Transform on 2D images or slices of a 3D image.
	 * 
	 * @param source
	 *              the source 2D image on which the Hough transform has to be done
	 * @param bitimg
	 *              the bitimg of the source to be used for doing the distance transform and
	 *              watershedding, created using user defined threshold
	 * @param minlength
	 *               the minimum length of the lines to be detected
	 *               this is done to avoid dots which appear in microscopy images.  
	 * 
	 */

	public HoughTransform2D(final RandomAccessibleInterval<T> source,
			final RandomAccessibleInterval<BitType> bitimg, final double minlength) {

		this.source = source;
		this.minlength = minlength;
		this.bitimg = bitimg;

	}
	
	
	
	
	@Override
	public boolean checkInput() {
		if ( source.numDimensions() > 2 )
		{
			errorMessage = BASE_ERROR_MSG + " Can only operate on 1D, 2D, make slices of your stack . Got " + source.numDimensions() + "D.";
			return false;
		}
		if ( minlength < 2 )
		{
			errorMessage = BASE_ERROR_MSG + "minimum length of line to be detected cannot be smaller than 2. Got " + minlength + ".";
			return false;
		}
		return true;
	}

	@Override
	public boolean process() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public RandomAccessibleInterval<T> getResult() {
		// TODO Auto-generated method stub
		return null;
	}

	
	/**
	 * 
	 * Method the get the seed image used for watershedding
	 * 
	 * @return An IntType image used as seeds to perform the watershedding.
	 */
	private RandomAccessibleInterval<IntType> Labelobjects() {

		// Prepare seed image for watershedding
		NativeImgLabeling<Integer, IntType> oldseedLabeling = new NativeImgLabeling<Integer, IntType>(
				new ArrayImgFactory<IntType>().create(source, new IntType()));

		oldseedLabeling = PrepareSeedImage(source);

		return oldseedLabeling.getStorageImg();
	}
	/***
	 * 
	 * Do the distance transform of the input image using the bit image
	 * provided.
	 * 
	 * @param inputimg
	 *            The pre-processed input image as RandomAccessibleInterval
	 *            <T>
	 * @param outimg
	 *            The distance transormed image having the same dimensions as
	 *            the input image.
	 * @param invtype
	 *            Straight: The intensity value is set to the distance, gives
	 *            white on black background. Inverse: The intensity is set to
	 *            the negative of the distance, gives black on white background.
	 */

	private void DistanceTransformImage(RandomAccessibleInterval<T> inputimg,
			RandomAccessibleInterval<T> outimg) {
		int n = inputimg.numDimensions();

		// make an empty list
		final RealPointSampleList<BitType> list = new RealPointSampleList<BitType>(n);

		// cursor on the binary image
		final Cursor<BitType> cursor = Views.iterable(bitimg).localizingCursor();

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
		final RandomAccess<T> ranac = outimg.randomAccess();

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
				
					ranac.get().setReal(search.getDistance());
					

			} else {
				// if value == 1, no need to search
				ranac.get().setZero();
			}
		}

	}

	private NativeImgLabeling<Integer, IntType> PrepareSeedImage(RandomAccessibleInterval<T> inputimg) {

		// New Labeling type
		final ImgLabeling<Integer, IntType> seedLabeling = new ImgLabeling<Integer, IntType>(
				new ArrayImgFactory<IntType>().create(inputimg, new IntType()));

		// Old Labeling type
		final NativeImgLabeling<Integer, IntType> oldseedLabeling = new NativeImgLabeling<Integer, IntType>(
				new ArrayImgFactory<IntType>().create(inputimg, new IntType()));

		// The label generator for both new and old type
		final Iterator<Integer> labelGenerator = AllConnectedComponents.getIntegerNames(0);

		// Getting unique labelled image (new version)
		ConnectedComponents.labelAllConnectedComponents(bitimg, seedLabeling, labelGenerator,
				ConnectedComponents.StructuringElement.EIGHT_CONNECTED);

		// Getting unique labelled image (old version)
		AllConnectedComponents.labelAllConnectedComponents(oldseedLabeling, bitimg, labelGenerator,
				AllConnectedComponents.getStructuringElement(inputimg.numDimensions()));

		return oldseedLabeling;
	}
	public static long[] GetMaxcorners(RandomAccessibleInterval<IntType> inputimg, int label) {

		Cursor<IntType> intCursor = Views.iterable(inputimg).localizingCursor();
		int n = inputimg.numDimensions();
		long[] maxVal = { Long.MIN_VALUE, Long.MIN_VALUE };

		while (intCursor.hasNext()) {
			intCursor.fwd();
			int i = intCursor.get().get();
			if (i == label) {

				for (int d = 0; d < n; ++d) {

					final long p = intCursor.getLongPosition(d);
					if (p > maxVal[d])
						maxVal[d] = p;

				}

			}
		}

		return maxVal;

	}

	public static long[] GetMincorners(RandomAccessibleInterval<IntType> inputimg, int label) {

		Cursor<IntType> intCursor = Views.iterable(inputimg).localizingCursor();
		int n = inputimg.numDimensions();
		long[] minVal = { Long.MAX_VALUE, Long.MAX_VALUE };

		while (intCursor.hasNext()) {
			intCursor.fwd();
			int i = intCursor.get().get();
			if (i == label) {

				for (int d = 0; d < n; ++d) {

					final long p = intCursor.getLongPosition(d);
					if (p < minVal[d])
						minVal[d] = p;
				}

			}
		}

		return minVal;

	}

	public static double GetBoundingbox(RandomAccessibleInterval<IntType> inputimg, int label) {

		Cursor<IntType> intCursor = Views.iterable(inputimg).localizingCursor();
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

		double boxsize = Distance(minVal, maxVal);

		Pair<long[], long[]> boundingBox = new Pair<long[], long[]>(minVal, maxVal);
		return boxsize;
	}

	public static int GetMaxlabelsseeded(RandomAccessibleInterval<IntType> intimg) {

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

	public NativeImgLabeling<Integer, IntType> GetlabeledImage(RandomAccessibleInterval<T> inputimg,
			NativeImgLabeling<Integer, IntType> seedLabeling) {

		int n = inputimg.numDimensions();
		long[] dimensions = new long[n];

		for (int d = 0; d < n; ++d)
			dimensions[d] = inputimg.dimension(d);
		final NativeImgLabeling<Integer, IntType> outputLabeling = new NativeImgLabeling<Integer, IntType>(
				new ArrayImgFactory<IntType>().create(inputimg, new IntType()));

		final Watershed<T, Integer> watershed = new Watershed<T, Integer>();

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

	public  RandomAccessibleInterval<T> CurrentLabelImage(RandomAccessibleInterval<IntType> Intimg,
			RandomAccessibleInterval<T> originalimg, int currentLabel) {

		RandomAccess<T> inputRA = originalimg.randomAccess();

		Cursor<IntType> intCursor = Views.iterable(Intimg).cursor();
		final T type = originalimg.randomAccess().get().createVariable();
		final ImgFactory< T > factory = Util.getArrayOrCellImgFactory( source, type );
		RandomAccessibleInterval<T> outimg = factory.create(originalimg, type);
		RandomAccess<T> imageRA = outimg.randomAccess();

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

	public static double Distance(final long[] minCorner, final long[] maxCorner) {

		double distance = 0;

		for (int d = 0; d < minCorner.length; ++d) {

			distance += Math.pow((minCorner[d] - maxCorner[d]), 2);

		}
		return Math.sqrt(distance);
	}
}

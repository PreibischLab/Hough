package Spindles;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;

import ij.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.labeling.AllConnectedComponents;
import net.imglib2.algorithm.stats.Normalize;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.labeling.NativeImgLabeling;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import preProcessing.GetLocalmaxmin;
import preProcessing.GlobalThresholding;
import util.ImgLib2Util;

// This class computes the min and max bounds of a labelled space and returns the maxextent between these co-ordinates
@SuppressWarnings("deprecation")
public class Boundingbox {

	public static final class Objectproperties {

		final int Label;
		final double maxextent;
		final int Area;

		protected Objectproperties(final int Label, final double maxextent, final int Area) {
			this.Label = Label;
			this.maxextent = maxextent;
			this.Area = Area;

		}
	}

	// This is the connected component bit, all the objects that are connected
	// in the image are given a unique label
	public static NativeImgLabeling<Integer, IntType> PrepareSeedImage(RandomAccessibleInterval<FloatType> inputimg) {

		// Preparing the seed image
		RandomAccessibleInterval<BitType> maximgBit = new ArrayImgFactory<BitType>().create(inputimg, new BitType());
		final Float threshold = GlobalThresholding.AutomaticThresholding(inputimg);
		GetLocalmaxmin.ThresholdingBit(inputimg, maximgBit, threshold);

		// Old Labeling type
		final NativeImgLabeling<Integer, IntType> oldseedLabeling = new NativeImgLabeling<Integer, IntType>(
				new ArrayImgFactory<IntType>().create(inputimg, new IntType()));

		// The label generator, I label the background to be 0
		final Iterator<Integer> labelGenerator = AllConnectedComponents.getIntegerNames(0);

		// Getting unique labelled image (old version)
		AllConnectedComponents.labelAllConnectedComponents(oldseedLabeling, maximgBit, labelGenerator,
				AllConnectedComponents.getStructuringElement(inputimg.numDimensions()));

		return oldseedLabeling;
	}

	// Once we have the label we get bounding box (BB for each of the labels, we
	// can only choose to get the BB for the largest region
	// after neglecting the background which carries the label 0
	public static void Getobjectprops(RandomAccessibleInterval<IntType> intimg, ArrayList<Objectproperties> listprops) {

		Cursor<IntType> intCursor = Views.iterable(intimg).localizingCursor();

		long[] position = new long[intimg.numDimensions()];

		// The background is labelled 0
		int currentLabel = 0;

		boolean anythingFound = true;
		while (anythingFound) {
			anythingFound = false;
			intCursor.reset();
			// Go through the whole image and add every pixel, that belongs to
			// the currently processed label

			long[] minVal = { Long.MAX_VALUE, Long.MAX_VALUE };
			long[] maxVal = { Long.MIN_VALUE, Long.MIN_VALUE };
			int count = 0;

			while (intCursor.hasNext()) {
				intCursor.fwd();

				int i = intCursor.get().get();
				// Now we are going inside a labelled space
				if (i == currentLabel) {
					anythingFound = true;
					count++;
					intCursor.localize(position);

					// Here we compute the bounding box for each labelled space
					for (int d = 0; d < intimg.numDimensions(); ++d) {
						if (position[d] < minVal[d]) {
							minVal[d] = position[d];
						}
						if (position[d] > maxVal[d]) {
							maxVal[d] = position[d];
						}

					}

				}

			}
			// Get the maxextent between min and max co-ordinates of the box
			// computed

			double maxextent = 0;
			for (int d = 0; d < intimg.numDimensions(); ++d)
				maxextent += (maxVal[d] - minVal[d]) * (maxVal[d] - minVal[d]);

			// Store all object properties in the java object and arraylist of
			// that object
			final Objectproperties props = new Objectproperties(currentLabel, Math.sqrt(maxextent), count);
			listprops.add(props);

			// Go to the next label
			currentLabel++;

		}

	}

	public static void main(String[] args) {

		RandomAccessibleInterval<FloatType> testimg = ImgLib2Util
				.openAs32Bit(new File("src/main/resources/test_10.tif"));

		new ImageJ();

		new Normalize();
		FloatType minval = new FloatType(0);
		FloatType maxval = new FloatType(1);
		Normalize.normalize(Views.iterable(testimg), minval, maxval);

		// Prepare seed image
		NativeImgLabeling<Integer, IntType> oldseedLabeling = new NativeImgLabeling<Integer, IntType>(
				new ArrayImgFactory<IntType>().create(testimg, new IntType()));

		oldseedLabeling = PrepareSeedImage(testimg);

		ImageJFunctions.show(oldseedLabeling.getStorageImg());
		ArrayList<Objectproperties> listprops = new ArrayList<Objectproperties>();

		// Execute the method to get object properties in the labelled space
		Getobjectprops(oldseedLabeling.getStorageImg(), listprops);

		// Display all object properties in the console
		for (int listindex = 0; listindex < listprops.size(); ++listindex) {
			final int Area = listprops.get(listindex).Area;
			final double maxextent = listprops.get(listindex).maxextent;
			final int Label = listprops.get(listindex).Label;
			System.out.println("Label: " + Label + " " + "Area: " + Area + " " + " maxextent: " + maxextent);

		}

	}

}

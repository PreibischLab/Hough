package varun;

import java.io.File;
import java.util.Iterator;

import com.sun.tools.javac.util.Pair;

import ij.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.labeling.AllConnectedComponents;
import net.imglib2.algorithm.labeling.ConnectedComponents;
import net.imglib2.algorithm.stats.Normalize;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.labeling.NativeImgLabeling;
import net.imglib2.roi.labeling.ImgLabeling;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import util.ImgLib2Util;
import varun.Kernels.ProcessingType;

// This class computes the min and max bounds of a labelled space and returns the distance between these co-ordinates
public class Boundingbox {
	
	// This is the connected component bit, all the objects that are connected in the image are given a unique label
	@SuppressWarnings("deprecation")
	public static NativeImgLabeling<Integer, IntType> PrepareSeedImage(
			RandomAccessibleInterval<FloatType> inputimg) 
	{

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

	// Once we have the label we get bounding box (BB for each of the labels, we can only choose to get the BB for the largest region
	// after neglecting the background which carries the label 0
	public static int LabelofChoice(RandomAccessibleInterval<IntType> intimg){
	
		Cursor<IntType> intCursor = Views.iterable(intimg).localizingCursor();
		
		long[] position = new long[intimg.numDimensions()];
		
   // The background is labelled 0 and is ignored here
	int currentLabel = 1;
	int maxcount = 0;
	int labelofinterest = 1;
	double finaldistance = 0;
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
		
		if (count>maxcount){
			maxcount =count;
			labelofinterest = currentLabel;
		// Get the distance between min and max co-ordinates of the box computed
		
			 double distance = 0;
			for (int d = 0; d < intimg.numDimensions(); ++d)
				distance+= (maxVal[d] - minVal[d]) * (maxVal[d] - minVal[d]);
			
		finaldistance = distance;
			
		}
		
		// Go to the next label
			currentLabel++;
			
		
		}

	System.out.println("distance: "+ Math.sqrt(finaldistance)+ " Label: "+ labelofinterest);
	
	return labelofinterest;
	}
	
	
	
	@SuppressWarnings("deprecation")
	public static void main(String[] args)  {

		RandomAccessibleInterval<FloatType> testimg = ImgLib2Util
				.openAs32Bit(new File("src/main/resources/test_10.tif"));
		
		new ImageJ();
	      
		new Normalize();
		FloatType minval = new FloatType(0);
		FloatType maxval = new FloatType(1);
		Normalize.normalize(Views.iterable(testimg), minval, maxval);
		
	
		
		
		// Prepare seed image for watershedding
				NativeImgLabeling<Integer, IntType> oldseedLabeling = new NativeImgLabeling<Integer, IntType>(
						new ArrayImgFactory<IntType>().create(testimg, new IntType()));

				oldseedLabeling = PrepareSeedImage(testimg);
				
				ImageJFunctions.show(oldseedLabeling.getStorageImg());
		final int label = LabelofChoice(oldseedLabeling.getStorageImg());
		
	
	
	}
	
}

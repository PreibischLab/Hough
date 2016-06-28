package preProcessing;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import net.imglib2.Cursor;
import net.imglib2.Interval;
import net.imglib2.IterableInterval;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.numeric.RealType;
import net.imglib2.view.Views;

public class MedianFilter {
	
	public static < T extends RealType<  T > & Comparable<T> > void medianFilter(
			final RandomAccessibleInterval< T > src, final RandomAccessibleInterval< T > dst, final int[] kernelDim )
	{
		medianFilter(Views.extendMirrorSingle(src), src, dst, kernelDim);
	}

	/** 
	 * Apply median filter to the whole image
	 * @param infSrc -- extended image
	 * @param SrcInterval -- interval for the input
	 * @param dst -- output image
	 * @param kernelDim -- the size of the kernel in each dimension (should be *o
	 */
	
	public static < T extends RealType<  T > & Comparable<T> > void medianFilter(
			final RandomAccessible< T > infSrc, final Interval srcInterval, final RandomAccessibleInterval< T > dst, final int[] kernelDim){

		final RandomAccessibleInterval<T> src = Views.interval(infSrc, srcInterval);
		final Cursor<T> cSrc = Views.iterable(src).localizingCursor();
		final RandomAccess<T> rDst = dst.randomAccess();

		final int n = src.numDimensions();

		// store kernel boundaries
		final long[] min = new long[n];
		final long[] max = new long[n];

		final int[] kernelHalfDim = new int[n]; 
		for (int d = 0; d < n; ++d){
			// check that the dimension of the kernel is correct
			if ( kernelDim[d]%2 == 0 || kernelDim[d] < 3)	
				throw new RuntimeException("kernelDim[d] should be odd and (>= 3) for each d. For d = " + d +  " kernelDim[d] = " + kernelDim[d] + ".");
			kernelHalfDim[d] = kernelDim[d]/2; // store dim/2
		}

		// store previous/current position of cursor
		final long[] pPos = new long[n]; 
		final long[] cPos = new long[n]; 
		cSrc.localize(cPos);
		updateKernelMinMax(min, max, cPos, kernelHalfDim, n);

		// contains all elements of the kernel
		IterableInterval<T> histogram = Views.interval(infSrc, min, max);
		List<T> histogramList = new ArrayList<T>(); 
		addAll(histogram, histogramList);

		final long[] localMin  = new long[n - 1];
		final long[] localMax  = new long[n - 1];

		while(cSrc.hasNext()){
			cSrc.localize(pPos);
			cSrc.fwd();
			cSrc.localize(cPos);

			final long checkDist = checkDist(pPos, cPos, n);	 // check if the cursor moved only by one step						
			if (checkDist == 0){ // moved too far
				// define new boundaries of the new kernel-window
				updateKernelMinMax(min, max, cPos, kernelHalfDim, n);
				// adjust the histogram window
				histogram = Views.interval(infSrc, min, max);
				histogramList.clear();					
				addAll(histogram, histogramList);
			}
			else{ // moved by one
				final int dim = (int)Math.abs(checkDist) - 1; 	 //
				final long step = (long) Math.signum(checkDist); // shows the direction of movement (+1/-1)

				min[ dim ] += step;
				max[ dim ] += step;

				final RandomAccessible<T> removeSlice = Views.hyperSlice(infSrc, dim, step < 0 ? max[dim] + 1: min[dim] - 1);
				final RandomAccessible<T> addSlice  = Views.hyperSlice(infSrc, dim, step < 0 ? min[dim] : max[dim]);		

				// remove one dimension
				for (int i = 0; i < n; ++i)
					if (i != dim){
						localMin[i > dim ? i - 1 : i] = min[i];
						localMax[i > dim ? i - 1 : i] = max[i];
					}

				final RandomAccessibleInterval<T> removeHistogram = Views.interval(removeSlice, localMin, localMax);
				final RandomAccessibleInterval<T> addHistogram  = Views.interval(addSlice, localMin, localMax);

				removeElements(removeHistogram, histogramList);
				addElements(addHistogram, histogramList);
			}

			// get/set the median value 
			rDst.setPosition(cSrc);
			rDst.get().set(histogramList.get(histogramList.size()/2));
		}
	}

	/**
	 * copy elements from the histogram to list 
	 * */
	public final static <T extends RealType<T> & Comparable<T>> void addAll(final IterableInterval<T> histogram, final List<T> list){
		for (final T h : histogram) 
			list.add(h.copy()); 			
		Collections.sort(list);	
	}

	/**
	 * adjust the min/max values for the kernel
	 * */
	public final static void updateKernelMinMax(final long[] min, final long[] max, final long[] position, final int[] kernelHalfDim, int n){
		for (int d = 0; d < n; ++d){
			min[d] = position[d] - kernelHalfDim[d];
			max[d] = position[d] + kernelHalfDim[d];
		}
	}

	/**
	 *  add every element from histogram to the list
	 */
	public final static <T extends RealType<T> & Comparable<T>> void addElements(final RandomAccessibleInterval<T> histogram, final List<T> list){

		final List<T> histogramArray = new ArrayList<T>(); 
		final List<T> returnArray = new ArrayList<T>(); 

		addAll(Views.iterable(histogram), histogramArray);

		int i = 0; // histogramArray index
		int j = 0; // list index

		final int histogramArraySize = histogramArray.size();
		final int listSize = list.size();

		while((i != histogramArraySize)  && (j != listSize)){
			if (histogramArray.get(i).compareTo(list.get(j)) > 0)
				returnArray.add(list.get(j++));
			else{
				if (histogramArray.get(i).compareTo(list.get(j)) < 0)
					returnArray.add(histogramArray.get(i++));
				else{ // equality case
					returnArray.add(list.get(j++));
					returnArray.add(histogramArray.get(i++));
				}
			}

			// flush the rest
			if (i == histogramArraySize)
				while(j < listSize)
					returnArray.add(list.get(j++));
			if (j == listSize)
				while(i < histogramArraySize)
					returnArray.add(histogramArray.get(i++));		
		}

		// is this copying cheaper than shifting elements in the ArrayList
		list.clear();	
		for(final T h : returnArray) 
			list.add(h.copy());
	}

	/**
	 *  remove elements from the list
	 */
	public final static <T extends RealType<T> & Comparable<T>> void removeElements(final RandomAccessibleInterval<T> histogram, final List<T> list){		
		final List<T> histogramArray = new ArrayList<T>(); 
		addAll(Views.iterable(histogram), histogramArray);

		int i = 0; // histogramArray index
		int j = 0; // list index

		int histogramArraySize = histogramArray.size();
		while(i != histogramArraySize){ // iterate while we do not go through all elements in histogram
			if (histogramArray.get(i).compareTo(list.get(j)) == 0){
				list.remove(j);
				i++; 
			}
			else
				j++; 
		}			
	}	



	/**
	 * This one checks if the cursor moved only by one pixel away
	 * @param pPos - position of the previous element
	 * @param cPos - position of the current element
	 * @param n - num of dimensions
	 * @return d+1 in which it moved, or zero if jumped to far
	 */
	public static long checkDist(final long[] pPos, final long[] cPos, final int n){
		long dim = -1;
		long dir = 0;

		for (int d = 0; d < n; ++d){
			final long dist = cPos[d] - pPos[d]; // dist > 0 if we moved forward, dist < 0 otherwise

			if (dist != 0){ // ?moved
				if((Math.abs(dist) != 1) || (dim != -1)){ //?too far or ?more than once
					return 0;
				}
				else{
					dim = d; 	// set the direction of movement
					dir = dist;  // set the step 
				}
			}
		}

		return (dim + 1)*dir;
	}
}

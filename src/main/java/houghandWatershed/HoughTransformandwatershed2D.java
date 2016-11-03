package houghandWatershed;

import java.util.ArrayList;
import com.sun.tools.javac.util.Pair;

import drawandOverlay.HoughPushCurves;
import drawandOverlay.OverlayLines;
import labeledObjects.Lineobjects;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.Point;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.BenchmarkAlgorithm;
import net.imglib2.algorithm.OutputAlgorithm;
import net.imglib2.algorithm.localextrema.RefinedPeak;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Util;
import net.imglib2.view.Views;
import preProcessing.GetLocalmaxmin;
import preProcessing.GlobalThresholding;

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

public class HoughTransformandwatershed2D extends BenchmarkAlgorithm
		implements OutputAlgorithm<Pair<RandomAccessibleInterval<IntType>, ArrayList<Lineobjects>>> {

	private static final String BASE_ERROR_MSG = "[HoughTransform2D] ";
	private final RandomAccessibleInterval<FloatType> source;
	private final RandomAccessibleInterval<BitType> bitimg;
	private final double minlength;
	private RandomAccessibleInterval<IntType> watershedimage;
	private ArrayList<Lineobjects> linelist;

	/**
	 * Instantiate a new Hough Transform object that does Hough Transform on 2D
	 * images or slices of a 3D image.
	 * 
	 * @param source
	 *            the source 2D image on which the Hough transform has to be
	 *            done
	 * @param bitimg
	 *            the bitimg of the source to be used for doing the distance
	 *            transform and watershedding, created using user defined
	 *            threshold
	 * @param minlength
	 *            the minimum length of the lines to be detected this is done to
	 *            avoid dots which appear in microscopy images.
	 * 
	 */

	public HoughTransformandwatershed2D(final RandomAccessibleInterval<FloatType> source,
			final RandomAccessibleInterval<BitType> bitimg, final double minlength) {

		this.source = source;
		this.minlength = minlength;
		this.bitimg = bitimg;

	}

	public RandomAccessibleInterval<IntType> getWatershedimage() {

		return watershedimage;
	}

	public ArrayList<Lineobjects> getLinelist() {

		return linelist;
	}

	@Override
	public boolean checkInput() {
		if (source.numDimensions() > 2) {
			errorMessage = BASE_ERROR_MSG + " Can only operate on 1D, 2D, make slices of your stack . Got "
					+ source.numDimensions() + "D.";
			return false;
		}
		if (minlength < 2) {
			errorMessage = BASE_ERROR_MSG + "minimum length of line to be detected cannot be smaller than 2. Got "
					+ minlength + ".";
			return false;
		}
		return true;
	}

	@Override
	public boolean process() {

		
		WatershedDistimg WaterafterDisttransform = new WatershedDistimg(source, bitimg);
		WaterafterDisttransform.checkInput();
		WaterafterDisttransform.process();
		watershedimage = WaterafterDisttransform.getResult();
		final int Maxlabel = WaterafterDisttransform.GetMaxlabelsseeded(watershedimage);
		
		

		System.out.println("Total labels: " + Maxlabel);

		final int ndims = source.numDimensions();
		ArrayList<RefinedPeak<Point>> ReducedMinlist = new ArrayList<RefinedPeak<Point>>(
				source.numDimensions());

		ImageJFunctions.show(watershedimage);
		final double[] sizes = new double[ndims];

		// Automatic threshold determination for doing the Hough transform
		Float val = GlobalThresholding.AutomaticThresholding(source);

		linelist = new ArrayList<Lineobjects>(source.numDimensions());

		for (int label = 1; label < Maxlabel - 1; label++) {

			System.out.println("Label Number:" + label);

			RandomAccessibleInterval<FloatType> outimg =  Boundingboxes.CurrentLabelImage(watershedimage, source, label);
			
			long[] minCorner = Boundingboxes.GetMincorners(watershedimage, label);
			long[] maxCorner = Boundingboxes.GetMaxcorners(watershedimage, label);

			FinalInterval intervalsmall = new FinalInterval(minCorner, maxCorner);

			 RandomAccessibleInterval<FloatType> outimgview = Views.interval(outimg, intervalsmall);
			//ImageJFunctions.show(outimg);
			
			final double area = Distance(minCorner, maxCorner);
			
			if (area < minlength * minlength){
			System.out.println("Skipping currentl label, not a real line here!");
				continue;
			}



			// Set size of pixels in Hough space
			int mintheta = 0;

			// Usually is 180 but to allow for detection of vertical
			// lines,allowing a few more degrees

			int maxtheta = 220;
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
			final RandomAccessibleInterval<FloatType> houghimage = new ArrayImgFactory<FloatType>().create(interval,
					new FloatType());

			HoughPushCurves.Houghspace(outimgview, houghimage, min, max, val);

			for (int d = 0; d < houghimage.numDimensions(); ++d)
				sizes[d] = houghimage.dimension(d);

			// Define Arraylist to get the slope and the intercept of the Hough
			// detected lines
			ArrayList<RefinedPeak<Point>> SubpixelMinlist = new ArrayList<RefinedPeak<Point>>(source.numDimensions());

			// Get the list of all the detections
			SubpixelMinlist = GetLocalmaxmin.HoughspaceMaxima(houghimage, interval, sizes, thetaPerPixel, rhoPerPixel);

			// Reduce the number of detections by picking One line per Label,
			// using the best detection for each label
			ReducedMinlist = OverlayLines.ReducedList(outimg, SubpixelMinlist, sizes, min, max);

			ArrayList<double[]> points = new ArrayList<double[]>();

			points = OverlayLines.GetRhoTheta(ReducedMinlist, sizes, min, max);

			/**
			 * This object has rho, theta, min dimensions, max dimensions of the
			 * label
			 * 
			 */
			final Lineobjects line = new Lineobjects(label, points, minCorner, maxCorner);

			linelist.add(line);
		}

		return true;
	}

	@Override
	public Pair<RandomAccessibleInterval<IntType>, ArrayList<Lineobjects>> getResult() {

		Pair<RandomAccessibleInterval<IntType>, ArrayList<Lineobjects>> linepair = new Pair<RandomAccessibleInterval<IntType>, ArrayList<Lineobjects>>(
				watershedimage, linelist);

		return linepair;
	}

	

	

	
	public static double Distance(final long[] minCorner, final long[] maxCorner) {

		double distance = 0;

		for (int d = 0; d < minCorner.length; ++d) {

			distance += Math.pow((minCorner[d] - maxCorner[d]), 2);

		}
		return Math.sqrt(distance);
	}
}

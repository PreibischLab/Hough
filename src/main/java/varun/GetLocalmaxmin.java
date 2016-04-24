package varun;

import java.awt.Color;
import java.util.ArrayList;

import com.sun.tools.javac.util.Pair;

import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Line;
import ij.gui.Overlay;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.Interval;
import net.imglib2.IterableInterval;
import net.imglib2.Localizable;
import net.imglib2.Point;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealPoint;
import net.imglib2.RealPointSampleList;
import net.imglib2.algorithm.dog.DogDetection;
import net.imglib2.algorithm.dog.DogDetection.ExtremaType;
import net.imglib2.algorithm.localextrema.RefinedPeak;
import net.imglib2.algorithm.neighborhood.Neighborhood;
import net.imglib2.algorithm.neighborhood.RectangleShape;
import net.imglib2.algorithm.region.hypersphere.HyperSphere;
import net.imglib2.algorithm.region.hypersphere.HyperSphereCursor;
import net.imglib2.algorithm.stats.Normalize;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.Cursor;
import net.imglib2.Point;
import net.imglib2.PointSampleList;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Intervals;
import net.imglib2.util.RealSum;
import net.imglib2.util.ValuePair;
import net.imglib2.view.Views;

public class GetLocalmaxmin {

	public static enum IntensityType {
		Gaussian, Original, One
	}

	protected IntensityType intensityType;

	public static ArrayList<RefinedPeak<Point>> Removesimilar(ArrayList<RefinedPeak<Point>> SubpixelMinlist,
			double thetatolerance, double rhotolerance) {
		/********
		 * The part below removes the close values in theta and rho coordinate
		 * (keeps only one of the multiple values)
		 ********/

		int j = 0;

		for (int i = 0; i < SubpixelMinlist.size(); ++i) {

			j = i + 1;
			while (j < SubpixelMinlist.size()) {

				if (Math.abs(SubpixelMinlist.get(i).getDoublePosition(0)
						- SubpixelMinlist.get(j).getDoublePosition(0)) <= thetatolerance
						&& Math.abs(SubpixelMinlist.get(i).getDoublePosition(1)
								- SubpixelMinlist.get(j).getDoublePosition(1)) <= rhotolerance) {

					SubpixelMinlist.remove(j);

				}

				else {
					++j;
				}

			}

		}

		return SubpixelMinlist;
	}

	public static ArrayList<RefinedPeak<Point>> RejectLines(RandomAccessibleInterval<FloatType> inputimg,
			ArrayList<RefinedPeak<Point>> SubpixelMinlist, double[] sizes, double[] min, double[] max, int pixeljump) {

		int n = inputimg.numDimensions();
		ArrayList<RefinedPeak<Point>> ReducedMinlist = new ArrayList<RefinedPeak<Point>>(n);

		for (int index = 0; index < SubpixelMinlist.size(); ++index) {
			double[] rhothetapoints = new double[n];

			rhothetapoints = TransformCordinates.transformfwd(new double[] {
					SubpixelMinlist.get(index).getDoublePosition(0), SubpixelMinlist.get(index).getDoublePosition(1) },
					sizes, min, max);
			double[] location = new double[n];
			double slope, intercept;

			slope = -1.0 / Math.tan(Math.toRadians(rhothetapoints[0]));
			intercept = rhothetapoints[1] / Math.sin(Math.toRadians(rhothetapoints[0]));
			final RandomAccessibleInterval<FloatType> imgout = new ArrayImgFactory<FloatType>().create(inputimg,
					new FloatType());

			// Draw the exact line for the detected Hough space parameters
			PushCurves.Drawexactline(imgout, slope, intercept, IntensityType.Gaussian);

			FloatType minval = new FloatType(0);
			FloatType maxval = new FloatType(255);
			Normalize.normalize(Views.iterable(imgout), minval, maxval);
			// ImageJFunctions.show(imgout);

			FloatType Pixelfwd = new FloatType(0);
			FloatType Pixelini = new FloatType(0);

			FloatType val = new FloatType(250);
			float Pixeldiff, Pixelsum;

			final RandomAccess<FloatType> outbound = inputimg.randomAccess();
			Cursor<FloatType> cursor = Views.iterable(imgout).localizingCursor();
			int count = 0;

			// Compare each detected line with the ROI in the input image to
			// reject or accept the line

			while (cursor.hasNext()) {

				cursor.fwd();
				cursor.localize(location);

				if (cursor.get().compareTo(val) >= 0) {
					outbound.setPosition(cursor);

					Pixelini = outbound.get();
				}
				for (int d = 0; d < n; ++d) {
					if (cursor.getDoublePosition(d) < inputimg.dimension(d) - pixeljump)
						cursor.jumpFwd(pixeljump);
					else
						break;
				}

				if (cursor.get().compareTo(val) >= 0) {
					outbound.setPosition(cursor);

					Pixelfwd = outbound.get();
				}

				Pixeldiff = Pixelfwd.get() - Pixelini.get();
				Pixelsum = Pixelfwd.get() + Pixelini.get();

				if (Pixeldiff >= 0 && Pixelsum >= 255)
					count++;
			}

			if (count == 0) {
				SubpixelMinlist.remove(index);
				System.out.println(" Removed Peak at :" + "Theta: " + rhothetapoints[0] + " Rho: " + rhothetapoints[1]);
			} else {
				ReducedMinlist.add(SubpixelMinlist.get(index));
			}

		}

		return ReducedMinlist;
	}

	public static void Thresholding(RandomAccessibleInterval<FloatType> img, RandomAccessibleInterval<FloatType> imgout,
			Float ThresholdValue, final IntensityType setintensity, double[] sigma) {

		final double[] backpos = new double[imgout.numDimensions()];
		final Cursor<FloatType> bound = Views.iterable(img).localizingCursor();

		final RandomAccess<FloatType> outbound = imgout.randomAccess();

		while (bound.hasNext()) {

			bound.fwd();

			outbound.setPosition(bound);

			if (bound.get().get()>(ThresholdValue)) {

				bound.localize(backpos);
				switch (setintensity) {

				case Original:
					outbound.get().set(bound.get());
					break;

				case Gaussian:
					AddGaussian.addGaussian(imgout, backpos, sigma, false);
					break;
				case One:
					outbound.get().setReal(1);
					break;
				default:
					outbound.get().setReal(1);
					break;

				}

			}

			else {

				outbound.get().setZero();

			}

		}
	}

	public static void ThresholdingBit(RandomAccessibleInterval<FloatType> img,
			RandomAccessibleInterval<BitType> imgout, FloatType ThresholdValue) {

		final double[] backpos = new double[imgout.numDimensions()];
		final Cursor<FloatType> bound = Views.iterable(img).localizingCursor();

		final RandomAccess<BitType> outbound = imgout.randomAccess();

		while (bound.hasNext()) {

			bound.fwd();

			outbound.setPosition(bound);

			if (bound.get().compareTo(ThresholdValue) > 0) {

				bound.localize(backpos);

				outbound.get().setReal(1);

			}

			else {

				outbound.get().setZero();

			}

		}
	}

	public static RandomAccessibleInterval<FloatType> FindandDisplayLocalMaxima(RandomAccessibleInterval<FloatType> img,
			ImgFactory<FloatType> imageFactory, final IntensityType setintensity, double[] sigma) {

		// Create a new image for the output
		RandomAccessibleInterval<FloatType> output = imageFactory.create(img, new FloatType());

		// define an interval that is span number of pixel smaller on each side
		// in each dimension
		int span = 1;

		Interval interval = Intervals.expand(img, -span);

		// create a view on the source with this interval
		img = Views.interval(img, interval);

		// create a Cursor that iterates over the source and checks in a
		// 8-neighborhood
		// if it is a maxima
		final Cursor<FloatType> center = Views.iterable(img).cursor();

		// instantiate a RectangleShape to access rectangular local
		// neighborhoods

		final RectangleShape shape = new RectangleShape(span, true);

		// iterate over the set of neighborhoods in the image
		for (final Neighborhood<FloatType> localNeighborhood : shape.neighborhoods(img)) {
			final FloatType centerValue = center.next();

			// keep this boolean true as long as no other value in the local
			// neighborhood
			// is smaller
			boolean isMaximum = true;

			// check if all pixels in the local neighborhood that are smaller
			for (final FloatType value : localNeighborhood) {
				// test if the center is smaller than the current pixel value
				if (centerValue.compareTo(value) < 0) {
					isMaximum = false;
					break;
				}
			}
			int n = img.numDimensions();
			double[] position = new double[n];
			if (isMaximum) {
				final RandomAccess<FloatType> outbound = output.randomAccess();
				outbound.setPosition(center);

				center.localize(position);
				switch (setintensity) {

				case Original:
					outbound.get().set(center.get());
					break;

				case Gaussian:
					AddGaussian.addGaussian(output, position, sigma, false);
					break;
				default:
					AddGaussian.addGaussian(output, position, sigma, false);
					break;

				}

			}
		}

		return output;
	}

	public static RandomAccessibleInterval<FloatType> FindandDisplayLocalMinima(RandomAccessibleInterval<FloatType> img,
			ImgFactory<FloatType> imageFactory, final IntensityType setintensity, double[] sigma) {

		// Create a new image for the output
		Img<FloatType> output = imageFactory.create(img, new FloatType());

		// define an interval that is span number of pixel smaller on each side
		// in each dimension
		int span = 1;

		Interval interval = Intervals.expand(img, -span);

		// create a view on the source with this interval
		img = Views.interval(img, interval);

		// create a Cursor that iterates over the source and checks in a
		// 8-neighborhood
		// if it is a maxima
		final Cursor<FloatType> center = Views.iterable(img).cursor();

		// instantiate a RectangleShape to access rectangular local
		// neighborhoods

		final RectangleShape shape = new RectangleShape(span, true);

		// iterate over the set of neighborhoods in the image
		for (final Neighborhood<FloatType> localNeighborhood : shape.neighborhoods(img)) {
			final FloatType centerValue = center.next();

			// keep this boolean true as long as no other value in the local
			// neighborhood
			// is smaller
			boolean isMinimum = true;

			// check if all pixels in the local neighborhood that are smaller
			for (final FloatType value : localNeighborhood) {
				// test if the center is smaller than the current pixel value
				if (centerValue.compareTo(value) >= 0) {
					isMinimum = false;
					break;
				}
			}
			double[] position = new double[img.numDimensions()];
			if (isMinimum) {

				final RandomAccess<FloatType> outbound = output.randomAccess();
				outbound.setPosition(center);
				center.localize(position);
				switch (setintensity) {

				case Original:
					outbound.get().set(center.get());
					break;

				case Gaussian:
					AddGaussian.addGaussian(output, position, sigma, false);
					break;
				default:
					AddGaussian.addGaussian(output, position, sigma, false);
					break;

				}

			}
		}

		return output;
	}

	public static ArrayList<RealPoint> FindLocalMinima(RandomAccessibleInterval<FloatType> img) {

		int n = img.numDimensions();

		ArrayList<RealPoint> Minlist = new ArrayList<RealPoint>(n);

		
		int span = 1;

		Interval interval = Intervals.expand(img, -span);

		img = Views.interval(img, interval);

		final Cursor<FloatType> center = Views.iterable(img).cursor();

		final RectangleShape shape = new RectangleShape(span, true);

		for (final Neighborhood<FloatType> localNeighborhood : shape.neighborhoods(img)) {
			final FloatType centerValue = center.next();

			boolean isMinimum = true;

			for (final FloatType value : localNeighborhood) {
				if (centerValue.compareTo(value) >= 0) {
					isMinimum = false;
					break;
				}
			}

			if (isMinimum) {
				RealPoint Minpoints = new RealPoint(center);
				Minlist.add(Minpoints);
			}
		}

		return Minlist;
	}

	public static ArrayList<RealPoint> FindLocalMaxima(RandomAccessibleInterval<FloatType> img) {

		int n = img.numDimensions();

		ArrayList<RealPoint> Maxlist = new ArrayList<RealPoint>(n);

		int span = 1;

		Interval interval = Intervals.expand(img, -span);

		img = Views.interval(img, interval);

		final Cursor<FloatType> center = Views.iterable(img).cursor();

		final RectangleShape shape = new RectangleShape(span, true);

		for (final Neighborhood<FloatType> localNeighborhood : shape.neighborhoods(img)) {
			final FloatType centerValue = center.next();

			boolean isMaximum = true;

			for (final FloatType value : localNeighborhood) {

				if (centerValue.compareTo(value) < 0) {
					isMaximum = false;
					break;
				}
			}

			if (isMaximum) {
				RealPoint Minpoints = new RealPoint(center);
				Maxlist.add(Minpoints);
			}
		}

		return Maxlist;
	}

	// Detect minima in Scale space
	public static ArrayList<RefinedPeak<Point>> ScalespaceMinima(RandomAccessibleInterval<FloatType> inputimg,
			FinalInterval interval, double thetaPerPixel, double rhoPerPixel, double minPeakValue, double smallsigma,
			double bigsigma) {
		
		ArrayList<RefinedPeak<Point>> SubpixelMinlist = new ArrayList<RefinedPeak<Point>>(inputimg.numDimensions());
		// Create a Dog Detection object in Hough space
		DogDetection<FloatType> newdog = new DogDetection<FloatType>(Views.extendMirrorSingle(inputimg), interval,
				new double[] { thetaPerPixel, rhoPerPixel }, smallsigma, bigsigma, DogDetection.ExtremaType.MINIMA,
				minPeakValue, true);

		// Detect minima in Scale space
		SubpixelMinlist = newdog.getSubpixelPeaks();

		return SubpixelMinlist;
	}

	public static Pair<FloatType, FloatType> computeMinMaxIntensity(final IterableInterval<FloatType> inputimg) {
		// create a cursor for the image (the order does not matter)
		final Cursor<FloatType> cursor = inputimg.cursor();

		// initialize min and max with the first image value
		FloatType type = cursor.next();
		FloatType min = type.copy();
		FloatType max = type.copy();

		// loop over the rest of the data and determine min and max value
		while (cursor.hasNext()) {
			// we need this type more than once
			type = cursor.next();

			if (type.compareTo(min) < 0) {
				min.set(type);

			}

			if (type.compareTo(max) > 0) {
				max.set(type);

			}
		}
		Pair<FloatType, FloatType> pair = new Pair<FloatType, FloatType>(min, max);
		return pair;
	}

	// Automatic thresholding done on the Normalized input image
	public static Float AutomaticThresholding(IterableInterval<FloatType> inputimg) {
		int n = inputimg.numDimensions();
		FloatType min = new FloatType();
		FloatType max = new FloatType();

		Float Threshold, ThresholdNew, Thresholdupdate;

		Pair<FloatType, FloatType> pair = new Pair<FloatType, FloatType>(min, max);
		pair = computeMinMaxIntensity(inputimg);

		Threshold = (pair.snd.get() - pair.fst.get()) / 2;
		PointSampleList<FloatType> listA = new PointSampleList<FloatType>(n);
		PointSampleList<FloatType> listB = new PointSampleList<FloatType>(n);
		Cursor<FloatType> cursor = inputimg.localizingCursor();

		while (cursor.hasNext()) {
			cursor.fwd();

			if (cursor.get().get() <= Threshold) {
				Point newpointA = new Point(n);
				newpointA.setPosition(cursor);
				listA.add(newpointA, cursor.get().copy());
			} else {
				Point newpointB = new Point(n);
				newpointB.setPosition(cursor);
				listB.add(newpointB, cursor.get().copy());
			}
		}
		final RealSum realSumA = new RealSum();
		long countA = 0;

		for (final FloatType type : listA) {
			realSumA.add(type.getRealDouble());
			++countA;
		}

		final double sumA = realSumA.getSum() / countA;

		final RealSum realSumB = new RealSum();
		long countB = 0;

		for (final FloatType type : listB) {
			realSumB.add(type.getRealDouble());
			++countB;
		}

		final double sumB = realSumB.getSum() / countB;

		ThresholdNew = (float) (sumA + sumB) / 2;
		Thresholdupdate = SegmentbyThresholding(inputimg, ThresholdNew);

		while (true) {

			ThresholdNew = SegmentbyThresholding(inputimg, Thresholdupdate);

			if (Math.abs(Thresholdupdate - ThresholdNew) < 1.0E-3)
				break;
			Thresholdupdate = ThresholdNew;

		}

		return ThresholdNew;

	}

	// Segment image by thresholding, used to determine automatic thresholding
	// level
	public static Float SegmentbyThresholding(IterableInterval<FloatType> inputimg, Float Threshold) {

		int n = inputimg.numDimensions();
		Float ThresholdNew;
		PointSampleList<FloatType> listA = new PointSampleList<FloatType>(n);
		PointSampleList<FloatType> listB = new PointSampleList<FloatType>(n);
		Cursor<FloatType> cursor = inputimg.localizingCursor();
		while (cursor.hasNext()) {
			cursor.fwd();

			if (cursor.get().get() < Threshold) {
				Point newpointA = new Point(n);
				newpointA.setPosition(cursor);
				listA.add(newpointA, cursor.get().copy());
			} else {
				Point newpointB = new Point(n);
				newpointB.setPosition(cursor);
				listB.add(newpointB, cursor.get().copy());
			}
		}
		final RealSum realSumA = new RealSum();
		long countA = 0;

		for (final FloatType type : listA) {
			realSumA.add(type.getRealDouble());
			++countA;
		}

		final double sumA = realSumA.getSum() / countA;

		final RealSum realSumB = new RealSum();
		long countB = 0;

		for (final FloatType type : listB) {
			realSumB.add(type.getRealDouble());
			++countB;
		}

		final double sumB = realSumB.getSum() / countB;

		ThresholdNew = (float) (sumA + sumB) / 2;

		return ThresholdNew;

	}
	// OverlayLines

	public static void Overlaylines(RandomAccessibleInterval<FloatType> inputimg,
			ArrayList<RefinedPeak<Point>> SubpixelMinlist, double[] sizes, double[] min, double[] max) {

		double[] points = new double[inputimg.numDimensions()];

		ImageStack stack = new ImageStack((int) inputimg.dimension(0), (int) inputimg.dimension(1));

		stack.addSlice(ImageJFunctions.wrap(inputimg, "").getProcessor());
		new ImageJ();

		ImagePlus imp = new ImagePlus("scale space hough", stack);
		imp.show();

		Overlay o = imp.getOverlay();

		if (o == null) {
			o = new Overlay();
			imp.setOverlay(o);
		}

		o.clear();

		for (int index = 0; index < SubpixelMinlist.size(); ++index) {
			points = TransformCordinates.transformfwd(new double[] { SubpixelMinlist.get(index).getDoublePosition(0),
					SubpixelMinlist.get(index).getDoublePosition(1) }, sizes, min, max);
			System.out.println(" Found Peaks at :" + "Theta: " + points[0] + " Rho: " + points[1]);

			Line newline = new Line(0, points[1] / Math.sin(Math.toRadians(points[0])), inputimg.dimension(0),
					points[1] / Math.sin(Math.toRadians(points[0]))
							- inputimg.dimension(0) / Math.tan(Math.toRadians(points[0])));

			newline.setStrokeColor(Color.RED);
			newline.setStrokeWidth(0.8);

			o.add(newline);
		}
		imp.updateAndDraw();
	}

	public static RandomAccessibleInterval<FloatType> GradientofImage(RandomAccessibleInterval<FloatType> inputimg) {

		final RandomAccessibleInterval<FloatType> gradientimg = new ArrayImgFactory<FloatType>().create(inputimg,
				new FloatType());
		Cursor<FloatType> cursor = Views.iterable(gradientimg).localizingCursor();
		RandomAccessible<FloatType> view = Views.extendMirrorSingle(inputimg);
		RandomAccess<FloatType> randomAccess = view.randomAccess();
		// iterate over all pixels
		while (cursor.hasNext()) {
			// move the cursor to the next pixel
			cursor.fwd();

			// compute gradient in each dimension
			double gradient = 0;

			for (int d = 0; d < inputimg.numDimensions(); ++d) {
				// set the randomaccess to the location of the cursor
				randomAccess.setPosition(cursor);

				// move one pixel back in dimension d
				randomAccess.bck(d);

				// get the value
				double v1 = randomAccess.get().getRealDouble();

				// move twice forward in dimension d, i.e.
				// one pixel above the location of the cursor
				randomAccess.fwd(d);
				randomAccess.fwd(d);

				// get the value
				double v2 = randomAccess.get().getRealDouble();

				
				gradient += ((v2 - v1) * (v2 - v1) ) / 4;
			}

			cursor.get().setReal(Math.sqrt(gradient));
		}

		return gradientimg;
	}
	
	public static RandomAccessibleInterval<FloatType> FindDirectionalLocalMaxima(RandomAccessibleInterval<FloatType> img,
			ImgFactory<FloatType> imageFactory, final IntensityType setintensity, double[] sigma) {

		RandomAccessibleInterval<FloatType> output = imageFactory.create(img, new FloatType());

		int span = 1;

		Interval interval = Intervals.expand(img, -span);

		img = Views.interval(img, interval);

		final Cursor<FloatType> center = Views.iterable(img).cursor();

		final RectangleShape shape = new RectangleShape(span, true);
		
		for (final Neighborhood<FloatType> localNeighborhood : shape.neighborhoods(img)) {
			final FloatType centerValue = center.next();

			boolean isMaximum = true;
        final double direction = Gradientdirection(center, img);
Cursor<FloatType>cursor =  localNeighborhood.cursor();
//while(cursor.hasNext()){
	// Only consider these points
	//(y-yp)/(x-xp)=direction;
//}
			for (final FloatType value : localNeighborhood) {
				if (centerValue.compareTo(value) < 0) {
					isMaximum = false;
					break;
				}
			}
			int n = img.numDimensions();
			double[] position = new double[n];
			if (isMaximum) {
				final RandomAccess<FloatType> outbound = output.randomAccess();
				outbound.setPosition(center);

				center.localize(position);
				switch (setintensity) {

				case Original:
					outbound.get().set(center.get());
					break;

				case Gaussian:
					AddGaussian.addGaussian(output, position, sigma, false);
					break;
				default:
					AddGaussian.addGaussian(output, position, sigma, false);
					break;

				}

			}
		}

		return output;
	}

	public static double  Gradientdirection(Localizable cursor, RandomAccessibleInterval<FloatType> inputimg) {
		
		RandomAccessible<FloatType> view = Views.extendMirrorSingle(inputimg);
		RandomAccess<FloatType> randomAccess = view.randomAccess();
	// iterate over all pixels
	double direction = 0;
	

		// compute gradient in each dimension
		double gradientx = 0;
		double gradienty = 0;
		
      // Works only for 2 dimensions
		for (int d = 0; d < inputimg.numDimensions(); ++d) {
			// set the randomaccess to the location of the cursor
			randomAccess.setPosition(cursor);

			// move one pixel back in dimension d
			randomAccess.bck(d);

			// get the value
			double v1 = randomAccess.get().getRealDouble();

			// move twice forward in dimension d, i.e.
			// one pixel above the location of the cursor
			randomAccess.fwd(d);
			randomAccess.fwd(d);

			// get the value
			double v2 = randomAccess.get().getRealDouble();

			if (d == 0)
				gradientx = ((v2 - v1)) / 2;
			if (d == 1)
				gradienty = ((v2 - v1)) / 2;
			
		}

		direction = gradienty/(gradientx+1.0E-200);
		
	return direction ;
	}
	
	
}

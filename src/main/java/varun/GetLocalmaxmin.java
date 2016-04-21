package varun;

import java.awt.Color;
import java.util.ArrayList;

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
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Intervals;
import net.imglib2.view.Views;

public class GetLocalmaxmin {

	public static enum IntensityType {
		Gaussian, Original
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
		//	 ImageJFunctions.show(imgout);

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
				
				
				if (cursor.get().compareTo(val)>=0){
					outbound.setPosition(cursor);
					Pixelini = outbound.get();
		}
				for (int d = 0; d<n; ++d ){
				if (cursor.getDoublePosition(d)<inputimg.dimension(d)-pixeljump)
					cursor.jumpFwd(pixeljump);
				else
					break;
				}
			
					if (cursor.get().compareTo(val)>=0){
						outbound.setPosition(cursor);
						Pixelfwd = outbound.get();
					}
					for (int d = 0; d<n; ++d ){
					if (cursor.getDoublePosition(d)<inputimg.dimension(d)-pixeljump)
						cursor.jumpFwd(-pixeljump);
					else
						break;
					}
					Pixeldiff = Pixelfwd.get() - Pixelini.get();
					Pixelsum = Pixelfwd.get() + Pixelini.get();
					if (Pixeldiff > 0 || Pixelsum > 255) 
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
			FloatType ThresholdValue, final IntensityType setintensity, double[] sigma) {

		final double[] backpos = new double[imgout.numDimensions()];
		final Cursor<FloatType> bound = Views.iterable(img).localizingCursor();

		final RandomAccess<FloatType> outbound = imgout.randomAccess();

		while (bound.hasNext()) {

			bound.fwd();

			outbound.setPosition(bound);

			if (bound.get().compareTo(ThresholdValue) > 0) {

				bound.localize(backpos);
				switch (setintensity) {

				case Original:
					outbound.get().set(bound.get());
					break;

				case Gaussian:
					AddGaussian.addGaussian(imgout, backpos, sigma, false);
					break;

				}

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

				}

			}
		}

		return output;
	}

	public static ArrayList<RealPoint> FindLocalMinima(RandomAccessibleInterval<FloatType> img) {

		int n = img.numDimensions();

		ArrayList<RealPoint> Minlist = new ArrayList<RealPoint>(n);

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
			newline.setStrokeWidth(0.3);

			o.add(newline);
		}
		imp.updateAndDraw();

	}

}

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

	public static ArrayList<RefinedPeak<Point>> Removesimilar(ArrayList<RefinedPeak<Point>> SubpixelMinlist) {
		/********
		 * The part below removes the close values in theta and rho coordinate
		 * (keeps only one of the multiple values)
		 ********/

		int j = 0;

		for (int i = 0; i < SubpixelMinlist.size(); ++i) {

			j = i + 1;
			while (j < SubpixelMinlist.size()) {

				if (Math.abs(SubpixelMinlist.get(i).getDoublePosition(0)
						- SubpixelMinlist.get(j).getDoublePosition(0)) == 0 
						&& Math.abs(SubpixelMinlist.get(i).getDoublePosition(1)
								- SubpixelMinlist.get(j).getDoublePosition(1)) == 0 ) {

					SubpixelMinlist.remove(j);

				}

				else {
					++j;
				}

			}

		}

		return SubpixelMinlist;
	}

	public static void Thresholding(RandomAccessibleInterval<FloatType> img, RandomAccessibleInterval<FloatType> imgout,
			FloatType ThresholdValue) {

		final Cursor<FloatType> bound = Views.iterable(img).localizingCursor();

		final RandomAccess<FloatType> outbound = imgout.randomAccess();

		while (bound.hasNext()) {

			bound.fwd();

			outbound.setPosition(bound);

			if (bound.get().compareTo(ThresholdValue) > 0) {

				//outbound.get().set(bound.get());
                  outbound.get().set(1);
			}

			else {

				outbound.get().setZero();

			}

		}
	}

	public static RandomAccessibleInterval<FloatType> FindandDisplayLocalMaxima(RandomAccessibleInterval<FloatType> img,
			ImgFactory<FloatType> imageFactory) {

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
			double[] sigma = new double[n];
			if (isMaximum) {
				 final RandomAccess<FloatType> outbound =
				 output.randomAccess();
				 outbound.setPosition(center);
				 outbound.get().set(center.get());
				 for (int d = 0; d<n; ++d){
					 position[d]= outbound.getDoublePosition(d);
					 sigma[d] = 1.0;
				 }
				 AddGaussian.addGaussian(output,position, sigma, false);
			}
		}

		return output;
	}

	public static RandomAccessibleInterval<FloatType> FindandDisplayLocalMinima(RandomAccessibleInterval<FloatType> img,
			ImgFactory<FloatType> imageFactory) {

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

			if (isMinimum) {
				 final RandomAccess<FloatType> outbound =
				 output.randomAccess();
				 outbound.setPosition(center);
				 outbound.get().set(center.get());

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

	// Do mean filtering on the inputimage
	public static void MeanFilter(RandomAccessibleInterval<FloatType> inputimage,RandomAccessibleInterval<FloatType> outimage,
			double sigma, FloatType ThresholdValue) {
		
		// Normalize the input image
		
		FloatType minval = new FloatType(0);
		FloatType maxval = new FloatType(255);
		Normalize.normalize(Views.iterable(inputimage), minval, maxval);
		
		
		
		// Mean filtering for a given sigma
		Cursor<FloatType> cursorInput = Views.iterable(inputimage).cursor();
		Cursor<FloatType> cursorOutput = Views.iterable(outimage).cursor();
		FloatType mean = Views.iterable(inputimage).firstElement().createVariable();
		while (cursorInput.hasNext()) {
			cursorInput.fwd();
			cursorOutput.fwd();
			HyperSphere<FloatType> hyperSphere = new HyperSphere<FloatType>(Views.extendMirrorSingle(inputimage),
					cursorInput, (long) sigma);
			HyperSphereCursor<FloatType> cursorsphere = hyperSphere.cursor();
			cursorsphere.fwd();
			mean.set(cursorsphere.get());
			int n = 1;
			while (cursorsphere.hasNext()) {
				cursorsphere.fwd();
				n++;
				mean.add(cursorsphere.get());
			}
			mean.div(new FloatType(n));
			cursorOutput.get().set(mean);
		}

	}
	
	// Detect minima in Scale space
	public static ArrayList<RefinedPeak<Point>> ScalespaceMinima(RandomAccessibleInterval<FloatType> inputimg, FinalInterval interval,
			double thetaPerPixel, double rhoPerPixel, double minPeakValue, double smallsigma, double bigsigma){
		ArrayList<RefinedPeak<Point>> SubpixelMinlist = new ArrayList<RefinedPeak<Point>>(inputimg.numDimensions());
		// Create a Dog Detection object in Hough space
				DogDetection<FloatType> newdog = new DogDetection<FloatType>(Views.extendMirrorSingle(inputimg),
						interval, new double[] { thetaPerPixel, rhoPerPixel }, smallsigma, bigsigma, DogDetection.ExtremaType.MINIMA,
						minPeakValue, false);

				// Detect minima in Scale space
				SubpixelMinlist = newdog.getSubpixelPeaks();
				
				// Remove duplicate  values in theta and rho
				SubpixelMinlist = GetLocalmaxmin.Removesimilar(SubpixelMinlist);
		
		
		return SubpixelMinlist;
	}
	
	// OverlayLines
	
	public static void Overlaylines(RandomAccessibleInterval<FloatType> inputimg, ArrayList<RefinedPeak<Point>> SubpixelMinlist,
			double[] sizes, double [] min, double[] max){
		
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
	

}

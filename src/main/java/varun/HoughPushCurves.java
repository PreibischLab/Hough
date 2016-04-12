package varun;

import java.awt.Color;
import java.io.File;
import java.util.ArrayList;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import net.imglib2.algorithm.dog.DifferenceOfGaussian;
import net.imglib2.algorithm.dog.DogDetection;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Line;
import ij.gui.Overlay;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.Point;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealCursor;
import net.imglib2.RealPoint;
import net.imglib2.RealPointSampleList;
import net.imglib2.algorithm.stats.Normalize;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import util.ImgLib2Util;

public class HoughPushCurves {

	public static void Houghspace(Img<FloatType> inputimage, Img<FloatType> imgout, double[] min, double[] max,
			FloatType threshold) {

		int n = inputimage.numDimensions();

		final long[] position = new long[n];
		double Amplitude, Phase;

		final Cursor<FloatType> inputcursor = inputimage.localizingCursor();

		// for every function (as defined by an individual pixel)
		while (inputcursor.hasNext()) {

			inputcursor.fwd();
			inputcursor.localize(position);
			if (inputcursor.get().compareTo(threshold) > 0) {
				Amplitude = Math.sqrt(Math.pow(position[0], 2) + Math.pow(position[1], 2));
				Phase = Math.toDegrees(Math.atan2(position[0], position[1]));
				// draw the function into the hough space

				PushCurves.DrawSine(imgout, min, max, Amplitude, Phase);
			}

		}
	}

	public static void main(String[] args) {

		final Img<FloatType> inputimg = ImgLib2Util.openAs32Bit(new File("src/main/resources/Horizontal_line.tif"));
		// Normalize the inputimg
		new Normalize();
		FloatType minval = new FloatType(0);
		FloatType maxval = new FloatType(255);
		Normalize.normalize(inputimg, minval, maxval);
		new ImageJ();
		ImageJFunctions.show(inputimg);
		// Set size of pixels in Hough space
		double thetaPerPixel = 0.5;
		double rhoPerPixel = 0.5;

		int mintheta = 0;
		int maxtheta = 180;
		double size = Math
				.sqrt((inputimg.dimension(0) * inputimg.dimension(0) + inputimg.dimension(1) * inputimg.dimension(1)));
		int minRho = (int) -Math.round(size);
		int maxRho = -minRho;

		double[] min = { mintheta, minRho };
		double[] max = { maxtheta, maxRho };

		int pixelsTheta = (int) Math.round((maxtheta - mintheta) / thetaPerPixel);
		int pixelsRho = (int) Math.round((maxRho - minRho) / rhoPerPixel);

		final double ratio = ((max[0] - min[0])) / ((max[1] - min[1]));
		// Size of Hough space
		FinalInterval interval = new FinalInterval(new long[] { pixelsTheta, (long) (pixelsRho * ratio) });
		final Img<FloatType> houghimage = new ArrayImgFactory<FloatType>().create(interval, new FloatType());
		
		ArrayList<Point> Minlist = new ArrayList<Point>(inputimg.numDimensions());
		FloatType val = new FloatType(100);
		
		final double[] sizes = new double[inputimg.numDimensions()];
		
		for (int d = 0; d < houghimage.numDimensions(); ++d)
			sizes[d] = houghimage.dimension(d);
		
		// Do the Hough transform
		Houghspace(inputimg, houghimage, min, max, val);

		ImageJFunctions.show(houghimage);

		// Create a Dog Detection object in Hough space
		DogDetection<FloatType> newdog = new DogDetection<FloatType>(Views.extendMirrorSingle(houghimage), interval,
				new double[] {pixelsTheta, pixelsRho},
				1.1, 1.1*1.1, DogDetection.ExtremaType.MINIMA, 10, false);
		
		// Detect minima in Scale space
		Minlist = newdog.getPeaks();
		
		double[] points = new double[inputimg.numDimensions()];
		for (int index = 0; index < Minlist.size(); ++index) {

			points = TransformCordinates.transformfwd(new double[] { Minlist.get(index).getDoublePosition(0),
					Minlist.get(index).getDoublePosition(1) }, sizes, min, max);

			System.out.println(" Found Peaks at :" +points[0] + " " + points[1]);
			
		}
		
		
		
	/*
	 * // Gettinglocal maximas
		final int numthreads = Runtime.getRuntime().availableProcessors();
		final ExecutorService service = Executors.newFixedThreadPool(numthreads);

		ImageStack stack = new ImageStack((int) houghimage.dimension(0), (int) houghimage.dimension(1));
		
	 
		for (double s = 1.1; s <= 15; s = s + 0.5) {
			System.out.println(s);
			final double[] sigmaA = { s, s };
			final double[] sigmaB = { sigmaA[0] * sigmaA[0], sigmaA[1] * sigmaA[1] };
			final Img<FloatType> dogimage = new ArrayImgFactory<FloatType>().create(interval, new FloatType());
			DifferenceOfGaussian.DoG(sigmaA, sigmaB, Views.extendMirrorSingle(houghimage), tmpdogimage, dogimage,
					service);
			stack.addSlice(ImageJFunctions.wrap(dogimage, "sigma=" + s).getProcessor());
			Minlist = GetLocalmaxmin.FindLocalMinima(dogimage);

			double[] points = new double[inputimg.numDimensions()];
			for (int index = 0; index < Minlist.size(); ++index) {

				points = TransformCordinates.transformfwd(new double[] { Minlist.get(index).getDoublePosition(0),
						Minlist.get(index).getDoublePosition(1) }, sizes, min, max);

				System.out.println(points[0] + " " + points[1]);
				System.out.println(
						Minlist.get(index).getDoublePosition(0) + " " + Minlist.get(index).getDoublePosition(1));
			}
		}

		new ImageJ();

		ImagePlus imp = new ImagePlus("scale space hough", stack);
		imp.show();
		/*
		 * // localminimage = GetLocalmaxmin.FindandDisplayLocalMinima(dogimage,
		 * new ArrayImgFactory<FloatType>());
		 * 
		 * new ImageJ(); ImageJFunctions.show(localminimage);
		 * 
		 * // Minlist = GetLocalmaxmin.FindLocalMinima(dogimage); final double[]
		 * sizes = new double[inputimg.numDimensions()]; RealCursor<FloatType>
		 * listcursor = Minlist.localizingCursor(); for (int d = 0; d <
		 * houghimage.numDimensions(); ++d) sizes[d] = houghimage.dimension(d);
		 * 
		 * double[] points = new double[inputimg.numDimensions()]; while
		 * (listcursor.hasNext()) { listcursor.fwd(); points =
		 * TransformCordinates.transformfwd( new double[] {
		 * listcursor.getDoublePosition(0), listcursor.getDoublePosition(1) },
		 * sizes, min, max);
		 * 
		 * System.out.println(points[0]+ " "+ points[1]);
		 * System.out.println(listcursor.getDoublePosition(0)+ " "+
		 * listcursor.getDoublePosition(1)); }
		 * 
		 * /* ImagePlus imp= new ImagePlus( "scale space hough", stack );
		 * imp.show();
		 * 
		 * Overlay o = imp.getOverlay();
		 * 
		 * if ( o == null ) { o = new Overlay(); imp.setOverlay( o ); }
		 * 
		 * o.clear();
		 * 
		 * for ( int i = 0; i < 2; ++i ) { Line l = new Line(10, 10, 100, 100 +
		 * i * 2); l.setStrokeColor( new Color( i, 128, 128 ) );
		 * l.setStrokeWidth(0.5); o.add( l ); }
		 * 
		 * 
		 * imp.updateAndDraw();
		 */
	}

}

package pSF;

import java.util.ArrayList;
import java.util.Random;

import ij.ImageJ;
import net.imglib2.FinalInterval;
import net.imglib2.Interval;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;


public class SimulateBeads {


	public static RandomAccessibleInterval<FloatType> renderPoints(final ArrayList<double[]> meanlist,
			final Interval interval, final double[] sigma) {
		final RandomAccessibleInterval<FloatType> imgout = new ArrayImgFactory<FloatType>().create(interval, new FloatType());

		for (int listindex = 0; listindex < meanlist.size(); ++listindex){
			
			drawandOverlay.AddGaussian.addGaussian(imgout, 1.0, meanlist.get(listindex), sigma);
			
		}
		

		return imgout;
	}



	public static ArrayList<double[]> randomPoints(final int numPoints, final Interval range, final Random rnd) {
		final ArrayList<double[]> points = new ArrayList<double[]>();

		for (int i = 0; i < numPoints; ++i) {
			final double[] p = new double[range.numDimensions()];

			for (int d = 0; d < range.numDimensions(); ++d)
				p[d] = rnd.nextDouble() * (range.max(d) - range.min(d)) + range.min(d);

			points.add(p);
		}

		return points;
	}

	
	public static void main(String[] args) {
		final int numPoints = 100;
		final double[] sigma = new double[] { 1, 1 };
		final Random rnd = new Random(500);
		final FinalInterval range = new FinalInterval(512, 512);
		 RandomAccessibleInterval<FloatType> imgout = new ArrayImgFactory<FloatType>().create(range, new FloatType());
		
		 ArrayList<double[]> meanlist = new ArrayList<double[]>();
		
		meanlist = randomPoints(numPoints, range, rnd);
		
		imgout = renderPoints(meanlist, range, sigma);
		new ImageJ();
		ImageJFunctions.show(imgout);

		
	}
}
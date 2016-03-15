package varun;

import java.io.File;

import ij.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import util.ImgLib2Util;

public class PlotFunctionPull {
	// function: (x-x0)^2 +(y-y0)^2 = R^2

	// for every pixel in the output image, compute the distance to the
	// function
// intensity = distance_function (circle)
	//%% in general intensity = sum_f gauss( distance_f ); f ... function

	public static <T extends RealType<T>> void pull(Img<T> inputimage, RandomAccessible<T> imgout, double[] min,
			double[] max) {

		int n = inputimage.numDimensions();
		final int[] position = new int[n];

		final double[] realpos = new double[n];

		double radius = 100, radiussecond =200,sigma =10, distance, functionone, functiontwo, distancesecond;

		double[] delta = new double[n];

		double[] center = new double[n];
		
		double[] centersecond = new double[n];

		for (int d = 0; d < n; ++d) {
			delta[d] = (max[d] - min[d]) / inputimage.dimension(d);
			center[d] = (inputimage.dimension(d) / 2) * delta[d] + min[d];
			centersecond[d] = (3*inputimage.dimension(d) / 4) * delta[d] + min[d];
			
		}

		final Cursor<T> inputcursor = inputimage.localizingCursor();
		final RandomAccess<T> outbound = imgout.randomAccess();

		while (inputcursor.hasNext()) {
			inputcursor.fwd();
			inputcursor.localize(position);

			// Forward transformation

			for (int d = 0; d < n; ++d) {

				realpos[d] = inputcursor.getDoublePosition(d) * delta[d] + min[0];

			}

			// Half circle function one

			functionone = Math.sqrt(radius * radius - (realpos[0] - center[0]) * (realpos[0] - center[0]))
					+ center[1];
			
			// Half circle function two 
			
			functiontwo = Math.sqrt(radiussecond * radiussecond - (realpos[0] - centersecond[0]) * (realpos[0] - centersecond[0]))
					+ centersecond[1];

			// Now write the closed form distance between the transformed points

			distance = Math.abs(Math.sqrt(((realpos[1] - center[1]) * (realpos[1] - center[1]))
					+ (realpos[0] - center[0]) * (realpos[0] - center[0])) - radius);
			distancesecond = Math.abs(Math.sqrt(((realpos[1] - centersecond[1]) * (realpos[1] - centersecond[1]))
					+ (realpos[0] - centersecond[0]) * (realpos[0] - centersecond[0])) - radiussecond);

			outbound.setPosition(inputcursor);
			
			
			outbound.get().setReal(Math.exp(-distance*distance/sigma)+Math.exp(-distancesecond*distancesecond/sigma));

		}
	}

	public static void main(String[] args) {

		final Img<FloatType> inputimage = ImgLib2Util.openAs32Bit(new File("src/main/resources/1line.tif"));

	//	ImageJFunctions.show(inputimage);

		double[] min = { -150, -150 };
		double[] max = { 150, 150 };

		int pixelsY = (int) Math.round((max[1] - min[1]));
		int pixelsX = (int) Math.round((max[0] - min[0]));

		FinalInterval intervalquad = new FinalInterval(new long[] { pixelsX, pixelsY });

		final Img<FloatType> houghquadimage = new ArrayImgFactory<FloatType>().create(intervalquad, new FloatType());

		pull(inputimage, houghquadimage, min, max);
//new ImageJ();
		ImageJFunctions.show(houghquadimage).setTitle("Pull-Circle function");

	}
}

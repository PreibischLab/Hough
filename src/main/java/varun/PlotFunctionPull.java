package varun;

import java.io.File;
import java.io.FileNotFoundException;

import javax.annotation.Generated;

import ij.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.Localizable;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.neighborsearch.NearestNeighborSearch;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import util.ImgLib2Util;





public class PlotFunctionPull {
	// function: (x-x0)^2 +(y-y0)^2 = R^2

	// for every pixel in the output image, compute the distance to the
	// function
// intensity = distance_function (circle)
	//%% in general intensity = sum_f gauss( distance_f ); f ... function

	public static <T extends RealType<T>> void pull(RandomAccessibleInterval<T> imgout, double[] min,
			double[] max) {

		int n = imgout.numDimensions();
		final int[] position = new int[n];

		final double[] realpos = new double[n];

		double radius = 20 , radiussecond = 20, sigma =1.5, cutoff = 5*sigma, distance, distancesecond;

		double[] delta = new double[n];

		double[] center = new double[n];
		
		double[] centersecond = new double[n];

		
		for (int d = 0; d < n; ++d) {
			delta[d] = (max[d] - min[d]) / (imgout.dimension(d));
			center[d] = 0; //(imgout.dimension(d) / 2)* delta[d] + min[d]; // This centers the circle at half the image size
			centersecond[d] = 0; // (3*imgout.dimension(d) / 4) * delta[d] + min[d];
			
		}

		final Cursor<T> inputcursor = Views.iterable( imgout ).localizingCursor();
		final RandomAccess<T> outbound = imgout.randomAccess();

		while (inputcursor.hasNext()) {
			inputcursor.fwd();
			inputcursor.localize(position);

			// Forward transformation

			for (int d = 0; d < n; ++d) {

				realpos[d] = position[d] * delta[d] + min[d];

			}
			
			// for the circle function (center and radius describes a circle)
			Finaldistance function = new Finalfunction(realpos, center, radius, 0);
			
			distance = function.Circlefunctiondist(realpos, center, radius);

			
			outbound.setPosition(inputcursor);
			
			
			if (Math.abs(distance) < cutoff)
				
		    outbound.get().setReal(Math.exp(-distance*distance/sigma));
			
			else 

			outbound.get().setReal(0);	
		}
		
		}
	
	

	public static void main(String[] args) throws FileNotFoundException {
	//	ImageJFunctions.show(inputimage);

		double[] min = { -50, -50 };
		double[] max = { 50, 50 };

		final double ratio = (max[1]-min[1]) / (max[0]-min[0]);
		final int sizeX = 100;
		final int sizeY =  (int)Math.round( sizeX * ratio ); 

		
		
		final Img<FloatType> houghquadimage = new ArrayImgFactory<FloatType>().create(new long[]{sizeX, sizeY}, new FloatType());
		
		
		
		pull(houghquadimage, min, max);
        
		//new ImageJ();
		ImageJFunctions.show(houghquadimage).setTitle("Pull-Circle function");

		
		
		
	}
}

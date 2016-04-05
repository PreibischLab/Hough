package varun;

import java.io.File;

import ij.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.Point;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealPoint;
import net.imglib2.RealRandomAccess;
import net.imglib2.RealRandomAccessible;
import net.imglib2.RealRandomAccessibleRealInterval;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.cell.CellImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.interpolation.randomaccess.NearestNeighborInterpolatorFactory;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import util.ImgLib2Util;

public class HoughTransform {

	public static <T extends RealType<T>> void Houghspace(Img<T> inputimage, RandomAccessible<T> imgout,
			double [] min, double [] max, T threshold) {

	
		int n = inputimage.numDimensions();

		final long[] position = new long[n];
		
		int[] point = new int[n];
		double [] delta = new double[n];

		final RandomAccess<T> outbound = imgout.randomAccess();
		final double[] realpos = new double[n];

		final Cursor<T> inputcursor = inputimage.localizingCursor();

		for (int d = 0; d < n; ++d) {
			delta[d] = (max[d] - min[d]) / (inputimage.dimension(d));
			
			
		}
		// for every function (as defined by an individual pixel)
		while (inputcursor.hasNext()) {

			inputcursor.fwd();
			inputcursor.localize(position);

			
			// Forward transformation

						for (int d = 0; d < n; ++d) {

							realpos[d] = position[d] * delta[d] + min[d];

						}
			
			// draw the function into the hough space
			for (int angle = 0; angle < max[1]; ++angle) {
				if (inputcursor.get().compareTo(threshold) > 0){
 
				
				double rho = Math.cos( Math.toRadians( angle ) ) * position[0] + Math.sin( Math.toRadians( angle ) ) * position[1];

			
				point[0] = angle;
				point[1] = (int)Math.round( (rho - min[1])/delta[1] );
				
				

				outbound.setPosition(point);
				outbound.get().set(inputcursor.get());
				
			}
		}
		}
	}

	public static void main(String[] args) {
		
		
		
		final Img<FloatType> inputimg = ImgLib2Util.openAs32Bit(new File("src/main/resources/vertical_line.tif"));
		ImageJFunctions.show(inputimg);
		double thetaPerPixel = 1;
		double rhoPerPixel = 1;
		int maxtheta = 360;
		
		double size = Math
				.sqrt((inputimg.dimension(0) * inputimg.dimension(0) + inputimg.dimension(1) * inputimg.dimension(1)));
		int minRho = (int) -Math.round( size ); 
		int maxRho = -minRho;
		
		double [] min = {0,minRho};
		double [] max = {maxtheta, maxRho};
		
		int pixelsTheta = (int)Math.round( maxtheta / thetaPerPixel );
		int pixelsRho = (int)Math.round( (maxRho - minRho) / rhoPerPixel );
		
		FinalInterval interval = new FinalInterval(new long[] { pixelsTheta, pixelsRho });

		final Img<FloatType> houghimage = new ArrayImgFactory<FloatType>().create(interval, new FloatType());

		FloatType val = new FloatType(200); 
		
		Houghspace(inputimg, houghimage, min, max, val);
 
		new ImageJ();
		
		ImageJFunctions.show(houghimage);

	}

}

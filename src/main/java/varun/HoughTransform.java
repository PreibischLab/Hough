package varun;

import java.io.File;

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
			double mintheta, int maxtheta, T threshold) {

		// mintheta = 0;
		// maxtheta = 180;
		int n = inputimage.numDimensions();

		final long[] position = new long[n];
		
		int[] point = new int[n];

		final RandomAccess<T> outbound = imgout.randomAccess();
		//for (int angle = 0; angle < maxtheta; ++angle) {
		//	theta[angle] = angle;
		//}

		final Cursor<T> inputcursor = inputimage.localizingCursor();

		
		// for every function (as defined by an individual pixel)
		while (inputcursor.hasNext()) {

			inputcursor.fwd();
			inputcursor.localize(position);

			// draw the function into the hough space
			for (int angle = 0; angle < maxtheta; ++angle) {
				if (inputcursor.get().compareTo(threshold) > 0){
 
				
				double rho = Math.cos( Math.toRadians( angle ) ) * position[0] + Math.sin( Math.toRadians( angle ) ) * position[1];

			
				
				point[0] = (int)Math.round( rho );
				point[1] = angle;
				

				outbound.setPosition(point);
				outbound.get().set(inputcursor.get());
				
			}
		}
		}
	}

	public static void main(String[] args) {
		
		
		
		final Img<FloatType> inputimg = ImgLib2Util.openAs32Bit(new File("src/main/resources/1line_short.tif"));

		double thetaPerPixel = 1;
		double rhoPerPixel = 1;
		int maxtheta = 180;
		
		double size = Math
				.sqrt((inputimg.dimension(0) * inputimg.dimension(0) + inputimg.dimension(1) * inputimg.dimension(1)));
		int minRho = (int) -Math.round( size ); 
		int maxrho = -minRho;
		
		int pixelsTheta = (int)Math.round( maxtheta / thetaPerPixel );
		int pixelsRho = (int)Math.round( (maxrho - minRho) / rhoPerPixel );
		
		FinalInterval interval = new FinalInterval(new long[] { pixelsRho, pixelsTheta });

		final Img<FloatType> houghimage = new ArrayImgFactory<FloatType>().create(interval, new FloatType());

		FloatType val = new FloatType(200); 
		
		Houghspace(inputimg, houghimage, 0, maxtheta, val);

		ImageJFunctions.show(houghimage);

	}

}

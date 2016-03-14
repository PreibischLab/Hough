package varun;

import java.io.File;

import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.Point;
import net.imglib2.RandomAccess;
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

	public static <T extends RealType<T>> void Houghspace(Img<T> inputimage, RealRandomAccessible<T> imgout,
			double mintheta, int maxtheta) {

		// mintheta = 0;
		// maxtheta = 180;
		int n = inputimage.numDimensions();

		final long[] position = new long[n];
		double deltatheta = Math.PI / maxtheta;

		double[] theta = new double[maxtheta];
		double[] rho = new double[maxtheta];
		double[] point = new double[maxtheta];

		final RealRandomAccess<T> outbound = imgout.realRandomAccess();
		for (int angle = 0; angle < maxtheta; ++angle) {
			theta[angle] = deltatheta * angle + mintheta;
		}

		final Cursor<T> inputcursor = inputimage.localizingCursor();

		T val = inputimage.firstElement().createVariable();
		while (inputcursor.hasNext()) {

			inputcursor.fwd();
			inputcursor.localize(position);

			for (int angle = 0; angle < maxtheta; ++angle) {
				if (inputcursor.get().compareTo(val) > 0){

					rho[angle] = Math.cos(theta[angle]) * position[0] + Math.sin(theta[angle]) * position[1];

				
				point[0] = angle;
				point[1] = rho[angle];
				

				outbound.setPosition(point);
//				outbound.get().set(inputcursor.get());
				outbound.get().setReal(position[0]);
			}
		}
		}
	}

	public static void main(String[] args) {
		final Img<FloatType> inputimg = ImgLib2Util.openAs32Bit(new File("src/main/resources/dt.png"));

		
		int maxtheta = 180;
		int maxrho = (int) ( Math
				.sqrt((inputimg.dimension(0) * inputimg.dimension(0) + inputimg.dimension(1) * inputimg.dimension(1))));
		
		FinalInterval interval = new FinalInterval(new long[] { maxrho, maxtheta });

		final Img<FloatType> houghimage = new ArrayImgFactory<FloatType>().create(interval, new FloatType());

		NearestNeighborInterpolatorFactory<FloatType> factory1 = new NearestNeighborInterpolatorFactory<FloatType>();
		RealRandomAccessible<FloatType> houghinterpolation = Views.interpolate(Views.extendMirrorSingle(houghimage),
				factory1);

		Houghspace(inputimg, houghinterpolation, 0, maxtheta);

		ImageJFunctions.show(houghimage);

	}

}

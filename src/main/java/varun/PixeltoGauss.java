package varun;

import java.io.File;
import java.io.FileNotFoundException;
import ij.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.IterableInterval;
import net.imglib2.Localizable;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealCursor;
import net.imglib2.RealLocalizable;
import net.imglib2.RealPointSampleList;
import net.imglib2.RealRandomAccess;
import net.imglib2.RealRandomAccessible;
import net.imglib2.RealRandomAccessibleRealInterval;
import net.imglib2.algorithm.region.hypersphere.HyperSphereCursor;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import util.ImgLib2Util;

public class PixeltoGauss {

	// Uses AddGaussian to add a Gaussian to a pixel position

	public static void main(String[] args) throws FileNotFoundException {
		final Img<FloatType> inputimg = ImgLib2Util.openAs32Bit(new File("src/main/resources/exact_circle_intensity.tif"));
		ImageJFunctions.show(inputimg);
		
        int n = inputimg.numDimensions();
		double[] position = new double[n];
		double[] sigma = {0.8,0.5};
		Cursor<FloatType> cursor = inputimg.localizingCursor();
		while(cursor.hasNext()){
			cursor.fwd();
			cursor.localize(position);
			
			AddGaussian.addGaussian(inputimg, position, sigma);
			
		}
		
		ImageJFunctions.show(inputimg).setTitle("Add Gauss");

	}
}

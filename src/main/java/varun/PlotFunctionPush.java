package varun;

import java.io.File;

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

public class PlotFunctionPush {
	// how you started now, draw one function after the other into the image

	// for (x-img.dimensions(0)/2)^2+(y-image.dimensions(1)/2)^2 = radius^2,
	// radius = image.dimensions(0)/2

	
		
		public static  <T extends RealType<T>> void circlefunction(Img<T> inputimage, RandomAccessible<T> imgout,
				T threshold, double radius) {

			int n = inputimage.numDimensions();

			final long[] position = new long[n];

			int[] point = new int[n];

			double[] center = new double [n];
			final RandomAccess<T> outbound = imgout.randomAccess();

			final Cursor<T> inputcursor = inputimage.localizingCursor();

			// for every function (as defined by an individual pixel)
			while (inputcursor.hasNext()) {

				inputcursor.fwd();
				inputcursor.localize(position);

				// draw the function space

				if (inputcursor.get().compareTo(threshold) > 0) {

					
                    
                    center[0] = inputimage.dimension(0)/2;
                    center[1] = inputimage.dimension(1)/2;
					double rho = Math.sqrt(radius*radius-(position[0]-center[0])*(position[0]-center[0]))+center[1] ;

					point[0] = (int) position[0];
					point[1] = (int) Math.round(rho);

					outbound.setPosition(point);
					outbound.get().set(inputcursor.get());

				}

			}
		}
	
	
	public static <T extends RealType<T>> void sqrtfunction(Img<T> inputimage, RandomAccessible<T> imgout,
			T threshold, double center) {

		int n = inputimage.numDimensions();

		final long[] position = new long[n];

		int[] point = new int[n];

		final RandomAccess<T> outbound = imgout.randomAccess();

		final Cursor<T> inputcursor = inputimage.localizingCursor();

		// for every function (as defined by an individual pixel)
		while (inputcursor.hasNext()) {

			inputcursor.fwd();
			inputcursor.localize(position);

			// draw the function space

			if (inputcursor.get().compareTo(threshold) > 0) {

				

				double rho =  Math.sqrt((position[0]))+center ;

				point[0] = (int) position[0];
				point[1] = (int) Math.round(rho);

				outbound.setPosition(point);
				outbound.get().set(inputcursor.get());

			}

		}
	}

	public static void main(String[] args) {

		final Img<FloatType> inputimg = ImgLib2Util.openAs32Bit(new File("src/main/resources/2lines.tif"));

		ImageJFunctions.show(inputimg);
		int maxcircleX = (int) inputimg.dimension(0);
 
		int maxsqrtX = maxcircleX;
		
		double radius = 100;
		double center = 10;
		
		int maxcircleY =   (int) (((inputimg.dimension(1)) )+radius);
		
		int maxsqrtY =  (int) (Math.sqrt(inputimg.dimension(0))+2*center);
		
		double YPerPixel = 1;
		double XPerPixel = 1;
		int pixelscircleY = (int) Math.round((maxcircleY) / YPerPixel);
		int pixelscircleX = (int) Math.round((maxcircleX) / XPerPixel);
		
		int pixelssqrtY = (int) Math.round((maxsqrtY) / YPerPixel);
		int pixelssqrtX = (int) Math.round((maxsqrtX) / XPerPixel);

		FinalInterval intervalcircle = new FinalInterval(new long[] { pixelscircleX, pixelscircleY });
		FinalInterval intervalsqrt = new FinalInterval(new long[] { pixelssqrtX, pixelssqrtY });

		final Img<FloatType> houghcircleimage = new ArrayImgFactory<FloatType>().create(intervalcircle, new FloatType());

		final Img<FloatType> houghsqrtimage = new ArrayImgFactory<FloatType>().create(intervalsqrt, new FloatType());
		
		FloatType val = new FloatType(200);

		circlefunction(inputimg, houghcircleimage, val, radius);

		ImageJFunctions.show(houghcircleimage).setTitle("Circle function");
		
		
		sqrtfunction(inputimg, houghsqrtimage, val,center );

		ImageJFunctions.show(houghsqrtimage).setTitle("Sqrt function");
		
		

	}
}

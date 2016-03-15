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

	public static <T extends RealType<T>> void push (Img<T> inputimage,RandomAccessible<T> imgout, double[] min, double[] max){
		
		int n = inputimage.numDimensions();
		final int[] position = new int[n];
		final int[] newpositionfunctone = new int[n];
		final int[] newpositionfuncttwo = new int[n];
		final int[] newpositionfunctthree = new int[n];
		final double[] realposfunctone = new double[n];
		final double[] realposfuncttwo = new double [n];
		final double[] realposfunctthree = new double [n];
		
	
		
		double[] delta = new double[n];
		
		for (int d = 0; d<n;++d){
		delta[d] = (max[d] - min[d])/inputimage.dimension(d);
		
		}
		
		double center =inputimage.dimension(0)/2*delta[0]+min[0]; // For centering function at the middle of the grid.
		
		final Cursor<T> inputcursor = inputimage.localizingCursor();
		final RandomAccess<T> outbound = imgout.randomAccess();
		
		while(inputcursor.hasNext()){
		inputcursor.fwd();
		inputcursor.localize(position);
		
		
			realposfunctone[0] = position[0]*delta[0]+min[0]; // Transforming from pixels to function space
			
			realposfunctone[1] = -100*Math.exp(-Math.pow((realposfunctone[0]-center),2)/100); // Gaussian function
			
			realposfuncttwo[0] = realposfunctone[0];
			realposfuncttwo[1] = 10*Math.sin(Math.toRadians(realposfunctone[0]*10)); // Sin function
			
			realposfunctthree[0] = realposfunctone[0];
			realposfunctthree[1] = -100*(-Math.pow((realposfunctone[0]-center),2)/100); // Quadratic function
			
		newpositionfunctone[1] = (int) Math.round((realposfunctone[1]-min[1])/delta[1]); // Transforming back the y coordinate for Gaussian
		newpositionfunctone[0] = position[0];
		
		
		newpositionfuncttwo[1] = (int) Math.round((realposfuncttwo[1]-min[1])/delta[1]); // Transforming back the y coordinate for Sin
		newpositionfuncttwo[0] = position[0];
		
		
		newpositionfunctthree[1] = (int) Math.round((realposfunctthree[1]-min[1])/delta[1]);// Transforming back the y coordinate for Quadratic
		newpositionfunctthree[0] = position[0];
		

		if (newpositionfunctone[1] < inputimage.dimension(1)){
		
	
		
		outbound.setPosition(newpositionfunctone);
		
		
	outbound.get().setReal(1);
		
		
		}
	
		if (newpositionfuncttwo[1] < inputimage.dimension(1)){
			
			
			
		outbound.setPosition(newpositionfuncttwo);
			
			
		outbound.get().setReal(1);
			
			
			}	
		if (newpositionfunctthree[1] < inputimage.dimension(1)){
			
			
			
			outbound.setPosition(newpositionfunctthree);
				
				
		outbound.get().setReal(1);
			
				
				}	
		
		
		
		}
	}

	

	public static void main(String[] args) {

		final Img<FloatType> inputimage = ImgLib2Util.openAs32Bit(new File("src/main/resources/1line.tif"));

		
		
		ImageJFunctions.show(inputimage);
		
		double[] min ={-150,-150};
		double[] max ={150,150};
		
		
		
		
		
		
		
		
		int pixelsY = (int) Math.round((max[1]-min[1]));
		int pixelsX = (int) Math.round((max[0]-min[0]));
		
		

		FinalInterval intervalquad = new FinalInterval(new long[] { pixelsX, pixelsY });
		
		

	
		
		final Img<FloatType> houghquadimage = new ArrayImgFactory<FloatType>().create(intervalquad, new FloatType());
		
		
		
		
		
		

		
		
	push(inputimage, houghquadimage, min, max);

		
		
		ImageJFunctions.show(houghquadimage).setTitle("Test function");

	}
}

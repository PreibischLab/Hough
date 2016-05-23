package varun;
import java.util.ArrayList;

import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import varun.PerformWatershedding.Lineobjects;
public class LengthDetection {

	
	public static void DetectLength(
			RandomAccessibleInterval<FloatType> inputimg,
			final RandomAccessibleInterval<FloatType> processedimg,
			final Img< FloatType > ransacWeight,
			Img<IntType> Intimg,
			final double[] sigma)
	{
		final int n = inputimg.numDimensions();
		
		final int Maxlabel = PerformWatershedding.GetMaxlabels(inputimg);
		
		RandomAccessibleInterval<FloatType> outimg = new ArrayImgFactory<FloatType>()
				.create(inputimg, new FloatType());
		
		//for (int label = 0; label < Maxlabel-1; ++label){
		int label = 5;
		Lineobjects line =	PerformWatershedding.Getlabelobject(inputimg, processedimg,label);
		final long[] startposition = new long[n];
		startposition[0] = line.boxXmin;
		startposition[1] = line.boxYmin;
		final double rho = line.Rho;
		final double theta = line.Theta;
		final double[] location = new double[n];
		outimg = PerformWatershedding.CurrentLabelImage(Intimg,inputimg, label);
		
		Cursor<FloatType> imgcursor = Views.iterable(outimg).localizingCursor();
		RandomAccess<FloatType> ranac = outimg.randomAccess();
		while(imgcursor.hasNext()){
			imgcursor.fwd();
			
			ranac.setPosition(startposition);
			ranac.move((long)(-imgcursor.getDoublePosition(0)/Math.tan(Math.toRadians(theta))+rho/Math.sin(Math.toRadians(theta))), 1);
			ranac.localize(location);
			GaussianMaskFit.gaussianMaskFit( inputimg,location,sigma, ransacWeight );
		}
		
		System.out.println(location[0]+" "+ location[1]);
	//	}
		
		
	}
	
}

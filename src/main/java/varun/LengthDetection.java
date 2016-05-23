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
			Img<IntType> Intimg, 
			ArrayList<Lineobjects> linelist, 
			final double[] sigma)
	{
		
		final int Maxlabel = PerformWatershedding.GetMaxlabels(inputimg);
		
		RandomAccessibleInterval<FloatType> outimg = new ArrayImgFactory<FloatType>()
				.create(inputimg, new FloatType());
		
		for (int label = 0; label < Maxlabel-1; ++label){
		outimg = PerformWatershedding.CurrentLabelImage(Intimg,inputimg, label);
		
		Cursor<FloatType> imgcursor = Views.iterable(outimg).localizingCursor();
		RandomAccess<FloatType> ranac = outimg.randomAccess();
		
		
		}
		
	//	GaussianMaskFit.gaussianMaskFit( inputimg,location,sigma, null );
		
	}
	
}

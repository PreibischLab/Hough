package pSF;

import java.util.ArrayList;

import drawandOverlay.PushCurves;
import houghandWatershed.PerformWatershedding;
import net.imglib2.Cursor;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import peakFitter.LengthDetection;
import preProcessing.GetLocalmaxmin;
import preProcessing.Kernels;
import preProcessing.GetLocalmaxmin.IntensityType;
import preProcessing.Kernels.ProcessingType;

public class Extractpsfinfo {

	private final RandomAccessibleInterval<FloatType> inputimg;
	private final int ndims;
	
	public Extractpsfinfo(RandomAccessibleInterval<FloatType> inputimg){
		
		this.inputimg = inputimg;
		this.ndims = inputimg.numDimensions();
		
	}
	
	public void Extractparams(ArrayList<double[]> totalgausslist) throws Exception{
		
		RandomAccessibleInterval<IntType> intimg = PerformWatershedding.Dowatersheddingonly(inputimg);
		
		RandomAccessibleInterval<FloatType> localmaximgout = new ArrayImgFactory<FloatType>().create(inputimg,
				new FloatType());
		
	
		
		localmaximgout = GetLocalmaxmin.FindandDisplayLocalMaxima(inputimg, IntensityType.Original, new double []{1,1});
		
		final double noiseLevel = 0.05;
		
		 Cursor<FloatType> localmaxcursor = Views.iterable(localmaximgout).localizingCursor();
			while(localmaxcursor.hasNext()){
				localmaxcursor.fwd();
				if (localmaxcursor.get().get()< noiseLevel)
					localmaxcursor.get().set(0);
			}
	
		ImageJFunctions.show(localmaximgout);
					
		PointSampleList<FloatType> centroidlist = new PointSampleList<FloatType>(ndims);
		
		PushCurves.Getcentroids(localmaximgout, centroidlist);
		
		 LengthDetection MTlength = new LengthDetection(inputimg,intimg);
		 
		 double[] final_param= new double[2*ndims+2];
		 final double [] typical_sigma = new double[ndims];
		 // Input the psf-sigma here to be used as a replacment for very large sigma in the solver
		 for (int d = 0; d < ndims; ++d)
			 typical_sigma[d] = 10; // To choose the data size around a point
		 
		 Cursor<FloatType> listcursor = centroidlist.localizingCursor();
		 while(listcursor.hasNext()){
			 listcursor.fwd();
			 final_param = MTlength.Getfinalpointsparam(listcursor, typical_sigma);
			 
			 if (final_param[3] > 0  && final_param[4] > 0 && 1/final_param[3]!=0 && 1/final_param[4]!=0 ){
			 totalgausslist.add(final_param);
			 
			 }
		 }
		 
		
			 
		 
			
		
	}
	
	
	
	
}

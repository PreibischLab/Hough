package pSF;

import java.util.ArrayList;

import houghandWatershed.WatershedDistimg;
import net.imglib2.Cursor;
import net.imglib2.Point;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;
import peakFitter.Gausswback;
import preProcessing.GetLocalmaxmin;

public class Extractpsfinfo {

	private final RandomAccessibleInterval<FloatType> inputimg;
	private final RandomAccessibleInterval<BitType> bitimg;
	private final int ndims;
	
	public Extractpsfinfo(RandomAccessibleInterval<FloatType> inputimg,
			RandomAccessibleInterval<BitType> bitimg){
		
		this.inputimg = inputimg;
		this.bitimg = bitimg;
		this.ndims = inputimg.numDimensions();
		
	}
	
	public void Extractparams(ArrayList<double[]> totalgausslist, final long radius,final boolean ignorebrightpeaks) throws Exception{
	
		WatershedDistimg Watershedobject = new WatershedDistimg(inputimg, bitimg);
		Watershedobject.checkInput();
		Watershedobject.process();
		RandomAccessibleInterval<IntType> intimg = Watershedobject.getResult();
	
		
		
		final int maxlabel = Watershedobject.GetMaxlabelsseeded(intimg); 
		PointSampleList<FloatType> centroidlist = new PointSampleList<FloatType>(ndims);
		 
		for (int label = 1; label < maxlabel - 1; ++label){
			
			long[] pos = GetLocalmaxmin.computeMaxinLabel(inputimg, intimg, label, ignorebrightpeaks);
			 
			
				// Ignore the boundary points
			
				if (pos != null && (pos[0] + radius) < inputimg.dimension(0) && (pos[0] - radius) > 0 &&
						(pos[1] + radius) < inputimg.dimension(1) && (pos[1] - radius) > 0){
					Point newpoint = new Point(pos);
         			centroidlist.add(newpoint, new FloatType(1));
				}
			
		}
		
			
		 
		 Gausswback  MTlength = new Gausswback(inputimg,intimg);
		 
		 double[] final_param= new double[2*ndims+2];
		
		 // Input the psf-sigma here to be used as a replacment for very large sigma in the solver
		 
		 
		 Cursor<FloatType> listcursor = centroidlist.localizingCursor();
		
		 while(listcursor.hasNext()){
			 listcursor.fwd();
			 final_param = MTlength.Getfinalpointsparam(listcursor, radius);
			 
			 totalgausslist.add(final_param);
			 
			 
		 }
		 
		 
			
		
	}
	
	
	
	
	
}

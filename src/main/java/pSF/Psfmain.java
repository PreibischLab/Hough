package pSF;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;

import com.sun.tools.javac.util.Pair;

import drawandOverlay.PushCurves;
import ij.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.gauss.Gauss;
import net.imglib2.algorithm.stats.Normalize;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import preProcessing.GetLocalmaxmin;
import preProcessing.Kernels;
import preProcessing.Kernels.ProcessingType;
import util.ImgLib2Util;

public class Psfmain {
	
	
	public static void main(String[] args) throws Exception {
		RandomAccessibleInterval<FloatType> biginputimg = ImgLib2Util
				.openAs32Bit(new File("src/main/resources/example_beads.tif"));
		
		new ImageJ();
      
		new Normalize();
		FloatType minval = new FloatType(0);
		FloatType maxval = new FloatType(1);
		Normalize.normalize(Views.iterable(biginputimg), minval, maxval);
	//	ImageJFunctions.show(biginputimg);
		RandomAccessibleInterval<FloatType> inputimg = new ArrayImgFactory<FloatType>().create(biginputimg,
				new FloatType());
		inputimg = biginputimg;
	
		// Initialize empty images to be used later
				RandomAccessibleInterval<FloatType> gaussimg = new ArrayImgFactory<FloatType>().create(inputimg,
						new FloatType());
				 final int ndims = inputimg.numDimensions();
				
				
		Extractpsfinfo getpsf = new Extractpsfinfo(inputimg);
		
		ArrayList<double[]> totalgausslist = new ArrayList<double[]>();
		

		final long radius = 8; //Raidus of the Hypersphere to choose data size around the point
		
		// Say true if you want to ignore the brightest beads, say false if you want to take all beads.
		getpsf.Extractparams(totalgausslist, radius, false);
		
		
		
		
		PushCurves.DrawDetectedGaussians(gaussimg, totalgausslist);	
		ImageJFunctions.show(gaussimg).setTitle("Iterated Result");
		
		
		
		for (int index = 0; index < totalgausslist.size(); ++index){
			
			
			System.out.println("Amp: " + 2 * totalgausslist.get(index)[0] + " " + "Mean X: "
					+ totalgausslist.get(index)[1] + " " + "Mean Y: " + totalgausslist.get(index)[2] + " " + "SigX: "
					+ Math.sqrt(1.0/totalgausslist.get(index)[3]) + " " + "SigY: "
					+ Math.sqrt(1.0/totalgausslist.get(index)[4])+  " "+ "Noise: "+ totalgausslist.get(index)[5]);
		
			
			
			
		}
		
		
		
	}
}

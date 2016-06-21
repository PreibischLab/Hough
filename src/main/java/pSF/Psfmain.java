package pSF;

import java.io.File;
import java.util.ArrayList;

import drawandOverlay.PushCurves;
import ij.ImageJ;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.gauss.Gauss;
import net.imglib2.algorithm.stats.Normalize;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import preProcessing.Kernels;
import preProcessing.Kernels.ProcessingType;
import util.ImgLib2Util;

public class Psfmain {
	
	
	public static void main(String[] args) throws Exception {
		RandomAccessibleInterval<FloatType> biginputimg = ImgLib2Util
				.openAs32Bit(new File("src/main/resources/Fresh_data/psf_488_11.tif"));
		// small_mt.tif image to be used for testing
		// 2015-01-14_Seeds-1.tiff for actual
		// mt_experiment.tif for big testing
		new ImageJ();
      
		new Normalize();
		FloatType minval = new FloatType(0);
		FloatType maxval = new FloatType(1);
		Normalize.normalize(Views.iterable(biginputimg), minval, maxval);
		ImageJFunctions.show(biginputimg);
		// Initialize empty images to be used later
				RandomAccessibleInterval<FloatType> inputimg = new ArrayImgFactory<FloatType>().create(biginputimg,
						new FloatType());
				RandomAccessibleInterval<FloatType> gaussimg = new ArrayImgFactory<FloatType>().create(inputimg,
						new FloatType());
				 final int n = inputimg.numDimensions();
				 
				
		Extractpsfinfo getpsf = new Extractpsfinfo(biginputimg);
		
		final double noiseLevel = 0.05;
		ArrayList<double[]> totalgausslist = new ArrayList<double[]>();
		getpsf.Extractparams(totalgausslist, noiseLevel);
		
		PushCurves.DrawDetectedGaussians(gaussimg, totalgausslist);	
		//ImageJFunctions.show(gaussimg).setTitle("Iterated Result");
		
		final double[] maxoneoversigma = { Double.MIN_VALUE, Double.MIN_VALUE };
		final double[] oneoversigma = new double[n];
		int maxindex = 0;
		for (int listindex = 0; listindex < totalgausslist.size(); ++listindex){
			
			for (int d = 0; d < n; ++d){
				
				oneoversigma[d] = totalgausslist.get(listindex)[n+d+1];
				
				if (oneoversigma[d] > maxoneoversigma[d]){
					
					maxoneoversigma[d] = oneoversigma[d];
					maxindex = listindex;
					
				}
				
			}
			
		}
		
		
		for (int index = 0; index < totalgausslist.size(); ++index){
			
			System.out.println("Amplitude: " + totalgausslist.get(index)[0] + " " + "Mean X: "
					+ totalgausslist.get(index)[1] + " " + "Mean Y: " + totalgausslist.get(index)[2] + " " + "SigmaX: "
					+ Math.sqrt(1.0/totalgausslist.get(index)[3]) + " " + "SigmaY: "
					+ Math.sqrt(1.0/totalgausslist.get(index)[4]));
			
			//System.out.println(Math.sqrt(1.0/totalgausslist.get(index)[3]));
		}
		
		 System.out.println("Printing the parameters for the PSF detected :");
		 
		 
		 System.out.println("Amplitude: " + totalgausslist.get(maxindex)[0] + " " + "Mean X: "
							+ totalgausslist.get(maxindex)[1] + " " + "Mean Y: " + totalgausslist.get(maxindex)[2] + " " + "SigmaX: "
							+ Math.sqrt(1.0/totalgausslist.get(maxindex)[3]) + " " + "SigmaY: "
							+ Math.sqrt(1.0/totalgausslist.get(maxindex)[4]));
		
	}
}

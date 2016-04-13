package varun;

import ij.ImageJ;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;

public class Drawpushcurves {
	
	public static void main(String[] args) {
	double[] min = { 0, -50};
	double[] max = { 540, 50 };

	final double ratio = (max[1]-min[1]) / (max[0]-min[0]);
	final int sizeX = 540;
	final int sizeY =  (int)Math.round( sizeX * ratio ); 
	final Img<FloatType> Sineimage = new ArrayImgFactory<FloatType>().create(new long[]{sizeX, sizeY}, new FloatType());
    double Amplitude = 10;
    double Phase = 0; // In Degrees
  	// Draw Sine
		PushCurves.DrawSine(Sineimage, min, max, Amplitude, Phase);
		
	     new ImageJ();
		ImageJFunctions.show(Sineimage).setTitle("Sine by push");
/*
	
    double [] center = {0,5};
    double radius =10;
    
    double[] mincircle = { -40, -10};
	double[] maxcircle = { 40, 40 };

	final double ratiocircle = (maxcircle[1]-mincircle[1]) / (maxcircle[0]-mincircle[0]);
	final int sizeXcircle = 800;
	final int sizeYcircle =  (int)Math.round( sizeXcircle * ratiocircle ); 
	
    final Img<FloatType> Circleimage = new ArrayImgFactory<FloatType>().create(new long[]{sizeXcircle, sizeYcircle}, new FloatType());

	
	// Draw Circle
	PushCurves.drawCircle(Circleimage, mincircle,
		 maxcircle, center,radius);
	new ImageJ();
	ImageJFunctions.show(Circleimage).setTitle("Circle by push");

	*/
}
}

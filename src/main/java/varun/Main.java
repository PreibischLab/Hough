package varun;

import ij.ImageJ;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;

public class Main {
	
	public static void main(String[] args) {
	double[] min = { -350, -20};
	double[] max = { 350, 20 };

	final double ratio = (max[1]-min[1]) / (max[0]-min[0]);
	final int sizeX = 700;
	final int sizeY =  (int)Math.round( sizeX * ratio ); 
double [] center = {0,5};
double radius =10;

	final Img<FloatType> houghquadimage = new ArrayImgFactory<FloatType>().create(new long[]{sizeX, sizeY}, new FloatType());

	//PushCurves.drawCircle(houghquadimage, min,
	//	 max, center,radius);
	
	PushCurves.DrawSine(houghquadimage, min, max, 10, 0.0);
	
     new ImageJ();
	ImageJFunctions.show(houghquadimage).setTitle("Moving along the curve");

}
}

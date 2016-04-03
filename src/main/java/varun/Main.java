package varun;

import ij.ImageJ;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;

public class Main {
	
	public static void main(String[] args) {
	double[] min = { -20, -10};
	double[] max = { 40, 40 };

	final double ratio = (max[1]-min[1]) / (max[0]-min[0]);
	final int sizeX = 800;
	final int sizeY =  (int)Math.round( sizeX * ratio ); 
double [] center = {0,5};
double radius =10;

	final Img<FloatType> houghquadimage = new ArrayImgFactory<FloatType>().create(new long[]{sizeX, sizeY}, new FloatType());

	PushCircle.drawCircle(houghquadimage, min,
		 max, center,radius);
new ImageJ();
	ImageJFunctions.show(houghquadimage).setTitle("Push-circle function");

}
}

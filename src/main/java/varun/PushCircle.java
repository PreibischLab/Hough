package varun;

import com.sun.tools.javah.Util.Exit;

import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;

public class PushCircle {

	public static void drawCircle(Img<FloatType> imgout, double[] min, double[] max, double[] center, double radius) {
		int n = imgout.numDimensions();
		double[] realpos = new double[n];
		double[] size = new double[n];
		double[] location = new double[n];
		double[] position = new double[n];
		double[] iniposition = new double[n];
		double[] backini = new double[n];
		double[] newpos = new double[n];
		double[] backpos = new double[n];
	
		for (int d = 0; d < n; ++d)
			size[d] = imgout.dimension(d);
		Cursor<FloatType> cursor = Views.iterable(imgout).localizingCursor();
		while (cursor.hasNext()) {
			cursor.fwd();
			cursor.localize(location);
			realpos = TransformCordinates.transformfwd(location, size, min, max);
			
			// To get a starting point on the circle
			if (Math.pow(realpos[0] - center[0], 2) + Math.pow(realpos[1] - center[1], 2)
					- radius * radius <=1.0E-50) {
				for (int d = 0; d < n; ++d)
					position[d] = realpos[d];
				break;
				
			}
		}

		final RandomAccess<FloatType> outbound = imgout.randomAccess();
		
		double stepsize = 0.1;
		
		for (int d = 0; d < n; ++d)
			iniposition[d] = position[d];
		double initheta =  Math.atan2(iniposition[1]-center[1],iniposition[0]-center[0]);
		double increment = Math.acos((2*radius*radius-stepsize)/(2*radius*radius));

		backini = TransformCordinates.transformback(iniposition, size, min, max);

		while (true) {
			
			// Move the current point along the curve
				
				newpos[0] = center[0]+radius* Math.cos((initheta-increment));
				newpos[1] = center[1]+radius* Math.sin((initheta-increment));
				initheta = Math.atan2(newpos[1]-center[1],newpos[0]-center[0]);
			
			// Transform the co-ordinates back as double[]
			backpos = TransformCordinates.transformback(newpos, size, min, max);
			int[] setpos = new int[n];
			setpos[0] = (int) Math.round(backpos[0]);
			setpos[1] = (int) Math.round(backpos[1]);

		//	System.out.println(Math.pow(newpos[0] - center[0], 2) + Math.pow(newpos[1] - center[1], 2)
		//	- radius * radius);
			
			// To set the pixel intensity
			AddGaussian.addGaussian(imgout, backpos, new double[] { 1, 1 });

			// To make sure that the values transformed back are not out of bounds
			if (backpos[0] < imgout.realMax(0) - imgout.realMin(0) || backpos[0] > imgout.realMin(0) || backpos[1] < imgout.realMax(1) - imgout.realMin(1)
					|| backpos[1] > imgout.realMin(1))

				outbound.setPosition(setpos);
			System.out.println(setpos[0]-(int) Math.round(backini[0])+" "+ (setpos[1]-(int) Math.round(backini[1])));
			// Stopping criteria of moving along the circular arc
			if ( Math.abs(setpos[0]-(int) Math.round(backini[0]))==0 &&  Math.abs(setpos[1]-(int) Math.round(backini[1]))==0   )
				break;
			
			
			 if (newpos[0] >= max[0] || newpos[0] <= min[0] || newpos[1] >= max[1] || newpos[1] <= min[1])

				break;
			

		}
	}

}

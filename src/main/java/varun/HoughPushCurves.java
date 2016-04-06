package varun;
import java.io.File;

import ij.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import util.ImgLib2Util;
public class HoughPushCurves {

		public static  void Houghspace(Img<FloatType> inputimage, Img<FloatType> imgout,
				double [] min, double [] max, FloatType threshold) {

		
			int n = inputimage.numDimensions();

			final long[] position = new long[n];
			double Amplitude, Phase;
			
			final Cursor<FloatType> inputcursor = inputimage.localizingCursor();

			// for every function (as defined by an individual pixel)
			while (inputcursor.hasNext()) {

				inputcursor.fwd();
				inputcursor.localize(position);
				if (inputcursor.get().compareTo(threshold) > 0){
				Amplitude = Math.sqrt(Math.pow(position[0],2) + Math.pow(position[1], 2));
				Phase = Math.toDegrees(Math.atan2(position[0], position[1]));
				// draw the function into the hough space
				
				PushCurves.DrawSine(imgout, min, max, Amplitude, Phase);
				}
				
			}
		}

		public static void main(String[] args) {
			
			
			
			final Img<FloatType> inputimg = ImgLib2Util.openAs32Bit(new File("src/main/resources/multiple_lines.tif"));
			ImageJFunctions.show(inputimg);
			double thetaPerPixel = 0.08;
			double rhoPerPixel = 1;
			int mintheta = 0;
			int maxtheta = 180;
			double size = Math
					.sqrt((inputimg.dimension(0) * inputimg.dimension(0) + inputimg.dimension(1) * inputimg.dimension(1)));
			int minRho = (int) -Math.round( size ); 
			int maxRho = -minRho;
			
			double [] min = {mintheta,minRho};
			double [] max = {maxtheta, maxRho};
			
			int pixelsTheta = (int)Math.round( (maxtheta-mintheta) / thetaPerPixel );
			int pixelsRho = (int)Math.round( (maxRho - minRho) / rhoPerPixel );
			
			FinalInterval interval = new FinalInterval(new long[] { pixelsTheta, pixelsRho });

			final Img<FloatType> houghimage = new ArrayImgFactory<FloatType>().create(interval, new FloatType());

			FloatType val = new FloatType(200); 
			
			Houghspace(inputimg, houghimage, min, max, val);
	 
			new ImageJ();
			
			ImageJFunctions.show(houghimage);

		}

	}

	


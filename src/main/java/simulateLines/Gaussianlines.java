package simulateLines;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;

import drawandOverlay.PushCurves;
import net.imglib2.Interval;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.numeric.real.FloatType;
import preProcessing.GetLocalmaxmin.IntensityType;

public class Gaussianlines {

	public static void Drawsimulatedlines(RandomAccessibleInterval<FloatType> outimg, final Interval range,
			final Random rnd,  final double[] sigma,  final int numlines) {

		final int n = outimg.numDimensions();

		// Declare the number of lines to be plotted
		ArrayList<Fakeline> linearray = new ArrayList<Fakeline>();

		for (int lineindex = 0; lineindex < numlines; ++lineindex) {

			double[] startpos = new double[n];
			double endpos[] = new double[n];

			for (int d = 0; d < range.numDimensions(); ++d) {
				startpos[d] = rnd.nextDouble() * (range.max(d) - range.min(d)) + range.min(d);

			}

			endpos[0] = startpos[0] + 20;
			endpos[1] = startpos[1] + rnd.nextDouble()*(endpos[0] - startpos[0]);
			
			double slope = (endpos[1] - startpos[1]) / (endpos[0] - startpos[0]);
			double intercept = startpos[1] - slope * startpos[0];

			double computedlength = 0;
			for (int d = 0; d < n ; ++d){
				
				computedlength += Math.pow((endpos[d] - startpos[d]),2);
				
			}
			
			PushCurves.Drawshortline(outimg, linearray, slope, intercept, startpos, endpos, sigma);
			
			try {
	            FileWriter writer = new FileWriter("initiallengths.txt", true);
	            writer.write( "StartX: "  + startpos[0]+  " " +
	           		 "StartY: "+ startpos[1] + " " + "EndposX: " + endpos[0] +  
	        		 " EndposY :" + endpos[1] + " Slope: " + slope 
	        				 +" Intercept: " + intercept + " Length " + Math.sqrt(computedlength) );
	            writer.write("\r\n"); 
	            writer.close();
	        } catch (IOException e) {
	            e.printStackTrace();
	        }

		}
		}
	

}

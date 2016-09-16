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

		ArrayList<Fakeline> linearray = new ArrayList<Fakeline>();

		for (int lineindex = 0; lineindex < numlines; ++lineindex) {

			double[] startpos = new double[n];
			double endpos[] = new double[n];

			for (int d = 0; d < range.numDimensions(); ++d) {
				startpos[d] = rnd.nextDouble() * (range.max(d) - range.min(d)) + range.min(d);

			}

			endpos[0] = startpos[0] + rnd.nextDouble()*2 + 20;
			endpos[1] = startpos[1] + 2*rnd.nextDouble()*(endpos[0] - startpos[0]) + 20;
			
			
			
			double slope = (endpos[1] - startpos[1]) / (endpos[0] - startpos[0]);
			double intercept = startpos[1] - slope * startpos[0];

			
			
			PushCurves.Drawshortline(outimg, linearray, slope, intercept, startpos, endpos, sigma);
			

		}
		
		for (int index = 0; index < linearray.size(); ++index){
		try {
	        FileWriter writer = new FileWriter("initiallengthsbig_fifth.txt", true);
	        writer.write( "StartX: "  + linearray.get(index).startpos[0]+  " " +
	       		 "StartY: "+ linearray.get(index).startpos[1] + " " + "EndposX: " + linearray.get(index).endpos[0] +  
	    		 " EndposY :" + linearray.get(index).endpos[1]+ "  Length " + linearray.get(index).length );
	        writer.write("\r\n"); 
	        writer.write("\r\n");
	        writer.close();
	    } catch (IOException e) {
	        e.printStackTrace();
	    }
		}
		}
	

}

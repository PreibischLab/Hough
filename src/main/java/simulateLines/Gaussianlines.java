package simulateLines;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;

import drawandOverlay.PushCurves;
import mISC.Tree.Distance;
import net.imglib2.Interval;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.exception.IncompatibleTypeException;
import net.imglib2.type.numeric.real.FloatType;
import preProcessing.GetLocalmaxmin.IntensityType;

public class Gaussianlines {

	public static void Drawsimulatedlines(RandomAccessibleInterval<FloatType> outimg, final Interval range,
			final Random rnd,  final double[] sigma,  final int numlines) throws IncompatibleTypeException {

		final int n = outimg.numDimensions();

		ArrayList<Fakeline> linearray = new ArrayList<Fakeline>();

		for (int lineindex = 0; lineindex < numlines; ++lineindex) {

			double[] startpos = new double[n];
			double endpos[] = new double[n];

			for (int d = 0; d < range.numDimensions(); ++d) {
				startpos[d] = rnd.nextDouble() * (range.max(d) - range.min(d)) + range.min(d);

			}

			endpos[0] = startpos[0] + rnd.nextDouble()*2.9 - 12;
			endpos[1] = startpos[1] + 4*rnd.nextDouble()*(endpos[0] - startpos[0]) + 20;
			
			
			
			double slope = (endpos[1] - startpos[1]) / (endpos[0] - startpos[0]);
			double intercept = startpos[1] - slope * startpos[0];

			
			
			//PushCurves.Drawshortline(outimg, linearray, slope, intercept, startpos, endpos, sigma);

			PushCurves.Drawshortline(outimg, linearray, slope, intercept, startpos, endpos, sigma);
			

		}
		
		for (int index = 0; index < linearray.size(); ++index){
		try {
	        FileWriter writer = new FileWriter("initiallengthsbig_eight.txt", true);
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
	
	public static ArrayList<Fakeline> Drawstartline(final Interval range, 
			final Random rnd,  final double[] sigma,  final int numlines){
		int n = range.numDimensions();

		ArrayList<Fakeline> linearrayini = new ArrayList<Fakeline>();
		
		for (int lineindex = 0; lineindex < numlines; ++lineindex) {

			
			
			double[] inistartpos = new double[n];
			double[] iniendpos = new double[n];

			
				
				for (int d = 0; d < range.numDimensions(); ++d) {
					inistartpos[d] = rnd.nextDouble() * (range.max(d) - range.min(d)) + range.min(d);

				}

				iniendpos[0] = inistartpos[0] + rnd.nextDouble()*2.9 - 12;
				iniendpos[1] = inistartpos[1] + 2*rnd.nextDouble()*(iniendpos[0] - inistartpos[0]) + 45 ;
				double inislope = (iniendpos[1] - inistartpos[1]) / (iniendpos[0] - inistartpos[0]);
				double iniintercept = inistartpos[1] - inislope * inistartpos[0];

				
				Fakeline line = new Fakeline(Distance(inistartpos, iniendpos), inislope, iniintercept, inistartpos, iniendpos);
				linearrayini.add(line);
		}
				for (int index = 0; index < linearrayini.size(); ++index){
					try {
				        FileWriter writer = new FileWriter("initial_length_move.txt", true);
				        writer.write(  "StartX: "  + linearrayini.get(index).startpos[0]+  " " +
				       		 "StartY: "+ linearrayini.get(index).startpos[1] + " " + "EndposX: " + linearrayini.get(index).endpos[0] +  
				    		 " EndposY :" + linearrayini.get(index).endpos[1]+ "  Length " + linearrayini.get(index).length );
				        writer.write("\r\n"); 
				        writer.write("\r\n");
				        writer.close();
				    } catch (IOException e) {
				        e.printStackTrace();
				    }
			
		}
		return linearrayini;
	}
	
	public static void Drawmovingsimulatedlines(RandomAccessibleInterval<FloatType> outimg, final Interval range, final int rate,
			final ArrayList<Fakeline> linearrayini,
			 final double[] sigma) throws IncompatibleTypeException {
		final Random velocity = new Random(10);
		final int n = outimg.numDimensions();

		ArrayList<Fakeline> linearray = new ArrayList<Fakeline>();
		
		for (int index = 0; index < linearrayini.size(); ++index) {

				
				
				double[] startpos = new double[n];
				double[] endpos = new double[n];
				
				
			for (int d = 0; d < range.numDimensions(); ++d) {
				startpos[d] = linearrayini.get(index).startpos[d];
				
			}
			endpos[0] = linearrayini.get(index).endpos[0]+ 4 * Math.sin(velocity.nextDouble()* rate);
			double inislope = (linearrayini.get(index).endpos[1] - linearrayini.get(index).startpos[1]) 
					/ (linearrayini.get(index).endpos[0] - linearrayini.get(index).startpos[0]);
			double iniintercept = linearrayini.get(index).startpos[1] - inislope * linearrayini.get(index).startpos[0];
			endpos[1] = inislope * endpos[0] +  iniintercept;
			
			
			double slope = (endpos[1] - startpos[1]) / (endpos[0] - startpos[0]);
			double intercept = startpos[1] - slope * startpos[0];

			
			
			PushCurves.Drawshortline(outimg, linearray, slope, intercept, startpos, endpos, sigma);
			

		}
	/*	
		for (int index = 0; index < linearray.size(); ++index){
		try {
	        FileWriter writer = new FileWriter("moving_length.txt", true);
	        writer.write("Frame:" + " " + rate + " " + "StartX: "  + linearray.get(index).startpos[0]+  " " +
	       		 "StartY: "+ linearray.get(index).startpos[1] + " " + "EndposX: " + linearray.get(index).endpos[0] +  
	    		 " EndposY :" + linearray.get(index).endpos[1]+ "  Length " + linearray.get(index).length );
	        writer.write("\r\n"); 
	        writer.write("\r\n");
	        writer.close();
	    } catch (IOException e) {
	        e.printStackTrace();
	    }
		}*/
		}
	public static double Distance(final double[] cordone, final double[] cordtwo) {

		double distance = 0;

		for (int d = 0; d < cordone.length; ++d) {

			distance += Math.pow((cordone[d] - cordtwo[d]), 2);

		}
		return Math.sqrt(distance);
	}
	
}

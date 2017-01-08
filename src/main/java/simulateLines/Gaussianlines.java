package simulateLines;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;

import com.sun.tools.javac.util.Pair;

import drawandOverlay.AddGaussian;
import drawandOverlay.PushCurves;
import mISC.Tree.Distance;
import net.imglib2.Interval;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.gauss3.Gauss3;
import net.imglib2.exception.IncompatibleTypeException;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import preProcessing.GetLocalmaxmin.IntensityType;

public class Gaussianlines {

	public static void Drawsimulatedlines(RandomAccessibleInterval<FloatType> outimg, final Interval range,
			 final double[] sigma,  final int numlines) throws IncompatibleTypeException {

		final int n = outimg.numDimensions();
		 // Pnoise1: (10, -10) Pnoise2: (20 , 10) Pnoise3:  (30, 18) Pnoise4: (15, 56) Pnoise5: (24, 44)
		ArrayList<Fakeline> linearray = new ArrayList<Fakeline>();
		final Random rnd = new Random(30);
		final Random rndsec = new Random(18);
		
		for (int lineindex = 0; lineindex < numlines; ++lineindex) {

			double startpos[] = new double[n];
			double endpos[] = new double[n];
			double MaxLength = 29.82;

			for (int d = 0; d < range.numDimensions(); ++d) {
				startpos[d] = (rnd.nextDouble() * (range.max(d) - range.min(d)) + range.min(d)) ;
				endpos[d] = ((rndsec.nextDouble() * (range.max(d) - range.min(d)) + range.min(d)))  ;
			}

			while (true){
			if (Distance(startpos, endpos) > MaxLength){
				
				for (int d = 0; d < range.numDimensions(); ++d) {
					
					endpos[d] = (startpos[d] + endpos[d]) / 2;
				}
				
			}
			if (Distance(startpos, endpos) <= MaxLength)
				break;
			}
			
			
			double slope = (endpos[1] - startpos[1]) / (endpos[0] - startpos[0]);
			double intercept = startpos[1] - slope * startpos[0];

		
			

			PushCurves.Drawshortline(outimg, slope, intercept, startpos, endpos, sigma);
			
			
		
			

			
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
	
	
	public static void GetSeeds(RandomAccessibleInterval<FloatType> outimg, ArrayList<double[]> startseeds,ArrayList<double[]> endseeds,
			
			final Interval range, final int numlines, final double[] sigma) throws IncompatibleTypeException{
		
		final Random rnd = new Random(-40);
		final Random rndsec = new Random(80);
		final int n = outimg.numDimensions();
		
	
		
			for (int index = 0; index < numlines; ++index) {
				
				double startpos[] = new double[n];
				double endpos[] = new double[n];
				double[] startline = new double[n];
				double[] endline = new double[n];
				double MaxLength = 55.82;
				
				for (int d = 0; d < range.numDimensions(); ++d) {
					startpos[d] = 150 + (rnd.nextDouble() * (range.max(d) - range.min(d)) + range.min(d)) ;
					endpos[d] = ((rndsec.nextDouble() * (range.max(d) - range.min(d)) + range.min(d)))  ;
				}
				final double[] tmppos = new double[n];
				
				final double[] minVal = new double[n];
				final double[] maxVal = new double[n];
				while (true){
					if (Distance(startpos, endpos) > MaxLength){
						
						for (int d = 0; d < range.numDimensions(); ++d) {
							
							endpos[d] = (startpos[d] + endpos[d]) / 2;
						}
						
					}
					if (Distance(startpos, endpos) <= MaxLength)
						break;
					}
					
					
					double slope = (endpos[1] - startpos[1]) / (endpos[0] - startpos[0]);
					double intercept = startpos[1] - slope * startpos[0];
				for (int d = 0; d < n; ++d) {

					final double locationdiff = startpos[d] - endpos[d];
					final boolean minsearch = locationdiff > 0;
					tmppos[d] = startpos[d];

					
						minVal[d] = minsearch ? endpos[d] : startpos[d];
						maxVal[d] = minsearch ? tmppos[d] : endpos[d];
					
					}

				
				
				if (slope >= 0) {
					for (int d = 0; d < n; ++d) {

						startline[d] = minVal[d];
						endline[d] = maxVal[d];
					}
					
				
				}

				if (slope < 0) {

					startline[0] = minVal[0];
					startline[1] = maxVal[1];
					endline[0] = maxVal[0];
					endline[1] = minVal[1];
				

				}
			
				double stepsize = 1.0;
				double steppos[] = {startline[0], startline[1]};
				double dx = stepsize / Math.sqrt(1 + slope * slope);
				double dy = slope * dx;
				while (true) {
					
					AddGaussian.addGaussian(outimg, steppos, sigma);
					
					if (steppos[0] > endline[0] || steppos[1] > endline[1] && slope >= 0)
						break;
					if (steppos[0] > endline[0] || steppos[1] < endline[1] && slope < 0)
						break;
					steppos[0] += dx;
					steppos[1] += dy;
					
					
				}
				for (int d  = 0; d < n ; ++d)
					endline[d] = steppos[d];
			
				
				final double[] startinfo = {startline[0], startline[1], slope, intercept};
				final double[] endinfo = {endline[0], endline[1], slope, intercept};
				startseeds.add(startinfo);
				endseeds.add(endinfo);
				
			
	}
				
			
			
			
		
		
		
		
	}
	
	
	public static Pair<ArrayList<double[]>, ArrayList<double[]>> Growseeds (RandomAccessibleInterval<FloatType> outimg, 
			ArrayList<double[]> startseeds, ArrayList<double[]> endseeds, int frame, double[] sigma) throws IncompatibleTypeException{
		
	
		final int n = outimg.numDimensions();
		
        double growrate = 12* Math.sin(0.2 * frame) ;

		
		 ArrayList<double[]> newcords = new ArrayList<double[]>();
		 ArrayList<double[]> newcordsend = new ArrayList<double[]>();
		
		 for (int index = 0; index < startseeds.size(); ++index){
			 newcords.add(startseeds.get(index));
			 double[] startpos = new double[n];
			 double[] endpos = new double[n];
			 double slope = startseeds.get(index)[n];
			 double intercept = startseeds.get(index)[n + 1];
			 
			 
			 for (int d = 0; d < n; ++d){
				 
				 endpos[d] = startseeds.get(index)[d];
				 
			 }
				 startpos[0] = startseeds.get(index)[0] -  1.5* Math.abs((growrate)) - 10.5 ;
				 startpos[1] = slope * startpos[0] + intercept;
				
			 
			    newcords.add(startpos);
			
			 
					final double stepsize =  1 ;
					double steppos[] = {startpos[0], startpos[1]};
					double dx = stepsize / Math.sqrt(1 + slope * slope);
					double dy = slope * dx;
					
					while (true) {
						
						AddGaussian.addGaussian(outimg, steppos, sigma);
						
						if (steppos[0] > endpos[0] || steppos[1] > endpos[1] && slope >= 0)
							break;
						if (steppos[0] > endpos[0] || steppos[1] < endpos[1] && slope < 0)
							break;
						steppos[0] += dx;
						steppos[1] += dy;
						
					}
					
				
			 
		 }
		 
		
       
		 for (int index = 0; index < endseeds.size(); ++index){
			 newcordsend.add(endseeds.get(index));
			 double[] startpos = new double[n];
			 double[] endpos = new double[n];
			 double slope = endseeds.get(index)[n];
			 double intercept = endseeds.get(index)[n + 1];
			 
			 for (int d = 0; d < n; ++d){
				 
				 startpos[d] = endseeds.get(index)[d];
			 }
				 endpos[0] = endseeds.get(index)[0] +  5.5* Math.abs((growrate)) ;
				endpos[1] = slope * endpos[0] + intercept;
			 
			 newcordsend.add(endpos);
			
			 PushCurves.Drawshortstrictcurve(outimg,  slope, intercept, startpos, endpos, sigma, frame);
			 
			 
		 }
	
	
		 Pair<ArrayList<double[]>, ArrayList<double[]>> pair = new Pair<ArrayList<double[]>, ArrayList<double[]>>(newcords, newcordsend);
	
	
	return pair;
		
	}
	
	
	
	public static Pair<Pair<ArrayList<double[]>, ArrayList<double[]>>, Pair<ArrayList<Double>, ArrayList<Double>>> Getgrowlength(
			ArrayList<double[]> startseeds, ArrayList<double[]> endseeds, int frame, double[] sigma) throws IncompatibleTypeException{
		
		
		final int n = sigma.length;
		
        final ArrayList<double[]> newstartpoints = new ArrayList<double[]>();
        final ArrayList<double[]> newendpoints = new ArrayList<double[]>();
        final ArrayList<Double> startlength = new ArrayList<Double>();
        final ArrayList<Double> endlength = new ArrayList<Double>();
		double growrate = 12* Math.sin(0.2 * frame) ;

		
		 
		 
		 for (int index = 0; index < startseeds.size(); ++index){
			 double[] startpos = new double[n];
			 double[] endpos = new double[n];
			 double slope = startseeds.get(index)[n];
			 double intercept = startseeds.get(index)[n + 1];
			 double length = 0;
			 
			 for (int d = 0; d < n; ++d){
				 
				 endpos[d] = startseeds.get(index)[d];
			 }
				 
				 startpos[0] = startseeds.get(index)[0] -  1.5* Math.abs((growrate)) - 10.5 ;
				startpos[1] = slope * startpos[0] + intercept;
					
			 
			 length = Distance(startpos, endpos);
				startlength.add(length);
			 final double[] startinfo = { startpos[0], startpos[1], slope, intercept};
			 newstartpoints.add(startinfo);
			 
			 
			 
			 
		 }
		 
		
        
		 for (int index = 0; index < endseeds.size(); ++index){
			 double[] startpos = new double[n];
			 double[] endpos = new double[n];
			 double slope = endseeds.get(index)[n];
			 double intercept = endseeds.get(index)[n + 1];
			 double length = 0;
			 
			 for (int d = 0; d < n; ++d){
				 
				 startpos[d] = endseeds.get(index)[d];
				 
			 }
				 endpos[0] = endseeds.get(index)[0] +  5.5* Math.abs((growrate)) ;
				 endpos[1] = slope * endpos[0] + intercept;
				
			 
			 length = Distance(startpos, endpos);
				endlength.add(length);
			 final double[] endinfo = { endpos[0], endpos[1], slope, intercept};
			 newendpoints.add(endinfo);
			
			 
			 
		 }
	
	
	Pair<ArrayList<double[]>, ArrayList<double[]>> pair = new Pair<ArrayList<double[]>, ArrayList<double[]>>(newstartpoints, newendpoints);
	Pair<ArrayList<Double>, ArrayList<Double>> pairlength = new Pair<ArrayList<Double>, ArrayList<Double>>(startlength, endlength);
	
	Pair<Pair<ArrayList<double[]>, ArrayList<double[]>>, Pair<ArrayList<Double>, ArrayList<Double>>> totalpair = new
			Pair<Pair<ArrayList<double[]>, ArrayList<double[]>>, Pair<ArrayList<Double>, ArrayList<Double>>>(pair, pairlength);
	
	return totalpair;
		
	}
	
	
	public static Pair<ArrayList<Pair<Integer, double[]>>, ArrayList<Pair<Integer, double[]>>> Drawmovingsimulatedlines(RandomAccessibleInterval<FloatType> outimg, final Interval range, final int rate,final int numlines,
			 final double[] sigma) throws IncompatibleTypeException {
		final int n = outimg.numDimensions();

		ArrayList<Pair<Integer, double[]>> linearraystart = new ArrayList<Pair<Integer, double[]>>();
		ArrayList<Pair<Integer, double[]>> linearrayend = new ArrayList<Pair<Integer, double[]>>();
		ArrayList<Fakeline> linearray = new ArrayList<Fakeline>();
		final Random rnd = new Random(-40);
		final Random rndsec = new Random(80);
		for (int index = 0; index < numlines; ++index) {

				
			double startpos[] = new double[n];
			double endpos[] = new double[n];
			double newstartpos[] = new double[n];
			double newendpos[] = new double[n];
			double MaxLength = 55.82;
			
			for (int d = 0; d < range.numDimensions(); ++d) {
				startpos[d] = 150 + (rnd.nextDouble() * (range.max(d) - range.min(d)) + range.min(d)) ;
				endpos[d] = ((rndsec.nextDouble() * (range.max(d) - range.min(d)) + range.min(d)))  ;
			}

			
			while (true){
			if (Distance(startpos, endpos) > MaxLength){
				
				for (int d = 0; d < range.numDimensions(); ++d) {
					
					endpos[d] = (startpos[d] + endpos[d]) / 2;
				}
				
			}
			if (Distance(startpos, endpos) <= MaxLength)
				break;
			}
			
			
			double slope = (endpos[1] - startpos[1]) / (endpos[0] - startpos[0]);
			double intercept = startpos[1] - slope * startpos[0];

		if (Math.abs(slope) < 10){
			
			if (rate == 0)
			PushCurves.Drawshortline(outimg, slope, intercept, startpos, endpos, sigma);
			final double[] tmppos = new double[n];
			final double[] startline = new double[n];
			double[] endline = new double[n];
			final double[] minVal = new double[n];
			final double[] maxVal = new double[n];
			for (int d = 0; d < n; ++d) {

				final double locationdiff = startpos[d] - endpos[d];
				final boolean minsearch = locationdiff > 0;
				tmppos[d] = startpos[d];

				
					minVal[d] = minsearch ? endpos[d] : startpos[d];
					maxVal[d] = minsearch ? tmppos[d] : endpos[d];
				
				}

			
			
			if (slope >= 0) {
				for (int d = 0; d < n; ++d) {

					startline[d] = minVal[d];
					endline[d] = maxVal[d];
				}
				
			
			}

			if (slope < 0) {

				startline[0] = minVal[0];
				startline[1] = maxVal[1];
				endline[0] = maxVal[0];
				endline[1] = minVal[1];
			

			}
			
			
			double growrate = 0;
		
					growrate = 12* Math.sin(0.2 * rate) ;
			
				
				
					
				newstartpos[0] = startline[0] -  1.5* Math.abs((growrate)) - 10.5 ;
				newstartpos[1] = slope * newstartpos[0] + intercept;
				
				newendpos[0] = endline[0] +  5.5* Math.abs((growrate)) ;
				newendpos[1] = slope * newendpos[0] + intercept;
				
				
				Pair<Integer, double[]> framepairstart = new Pair<Integer, double[]>(rate,newstartpos);
				Pair<Integer, double[]> framepairend = new Pair<Integer, double[]>(rate,newendpos);
				
				
				
				linearraystart.add(framepairstart);
				linearrayend.add(framepairend);
				
			if (rate > 0){
				
				
				PushCurves.Drawshortstrictline(outimg,  slope, intercept, newstartpos, startline, sigma, rate);
				PushCurves.Drawshortstrictcurve(outimg,  slope, intercept, endline, newendpos, sigma, rate);
				
			}
			
			
		}
		}
		Pair<ArrayList<Pair<Integer, double[]>>, ArrayList<Pair<Integer, double[]>> > pair = new Pair<ArrayList<Pair<Integer, double[]>>, ArrayList<Pair<Integer, double[]>> >(linearraystart, linearrayend);
		
		return pair;
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
	    }*/
		//}
		}
	public static double Distance(final double[] cordone, final double[] cordtwo) {

		double distance = 0;

		for (int d = 0; d < cordone.length; ++d) {

			distance += Math.pow((cordone[d] - cordtwo[d]), 2);

		}
		return Math.sqrt(distance);
	}
	
}

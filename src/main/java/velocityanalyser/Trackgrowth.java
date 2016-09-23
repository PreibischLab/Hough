package velocityanalyser;

import java.util.ArrayList;

import com.sun.tools.javac.util.Pair;

import graphconstructs.Staticproperties;
import houghandWatershed.PerformWatershedding;
import net.imglib2.Point;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.IntervalView;
import peakFitter.Linefitter;

public class Trackgrowth {

	private final IntervalView<FloatType> currentframe;
	private final double minlength;
	private final ArrayList<double[]> PrevFrameparam;
	private final int framenumber;
	private final double[] psf;
	
	public Trackgrowth(final IntervalView<FloatType> currentframe, final double minlength, final ArrayList<double[]> PrevFrameparam,
			final int framenumber, final double[] psf){
		
		this.currentframe = currentframe;
		this.minlength = minlength;
		this.PrevFrameparam = PrevFrameparam;
		this.framenumber = framenumber;
		this.psf = psf;
		
	}
	
	public Pair<ArrayList<double[]>, ArrayList<Staticproperties>> Updatetrackpoints() throws Exception{
		
		PerformWatershedding Watershedobject = new PerformWatershedding(currentframe, minlength);
		RandomAccessibleInterval<IntType> currentlabelledimg = Watershedobject.Dowatersheddingonly();
		
		Linefitter currentline = new Linefitter(currentframe, currentlabelledimg);
		int ndims = currentframe.numDimensions();
		
		ArrayList<double[]> NextFrameparam = new ArrayList<double[]>();
		
		ArrayList<Staticproperties> startandendinframe = new ArrayList<Staticproperties>();
		
		
		for (int index = 0; index < PrevFrameparam.size(); ++index) {

			Point linepoint = new Point(ndims);
			linepoint.setPosition(
					new long[] { (long) PrevFrameparam.get(index)[0], (long) PrevFrameparam.get(index)[1] });
			
			
			
			 int currentlabel = currentline.Getlabel(linepoint);
			
			 double[] paramnextframe =
					currentline.Getfinaltrackparam(PrevFrameparam.get(index),
							currentlabel, psf, framenumber, true);
			 NextFrameparam.add(paramnextframe);
			 
			 
			 final double[] oldstartpoint = {PrevFrameparam.get(index)[0], PrevFrameparam.get(index)[1]};
			 
			 final double[] oldendpoint = {PrevFrameparam.get(index)[2], PrevFrameparam.get(index)[3]};
			 
			 final double[] newstartpoint = {paramnextframe[0], paramnextframe[1]};
			 
			 final double[] newendpoint = {paramnextframe[2], paramnextframe[3]};
			 
			 final double[] directionstart = {newstartpoint[0] - oldstartpoint[0] , newstartpoint[1] - oldstartpoint[1] };
			 
			 final double[] directionend = {newendpoint[0] - oldendpoint[0] , newendpoint[1] - oldendpoint[1] };
			 
			System.out.println("Frame:" + framenumber + " " +  "Fits :" + currentlabel + " "+ "StartX:" + paramnextframe[0] 
					+ " StartY:" + paramnextframe[1] + " " + "EndX:"
					+ paramnextframe[2] + "EndY: " + paramnextframe[3]);
			
			final Staticproperties edge = 
		   new Staticproperties(currentlabel, oldstartpoint, oldendpoint, newstartpoint, newendpoint, directionstart , directionend );
			

					startandendinframe.add(edge);	
		}
		
		Pair<ArrayList<double[]> , ArrayList<Staticproperties>> pair = 
		new Pair<ArrayList<double[]> , ArrayList<Staticproperties>>(NextFrameparam, startandendinframe);
		
	
		return pair;
		
	}
	
	
	
}

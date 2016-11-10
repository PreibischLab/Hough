package peakFitter;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import graphconstructs.Staticproperties;
import houghandWatershed.Boundingboxes;
import ij.gui.EllipseRoi;
import labeledObjects.LabelledImg;
import net.imglib2.Cursor;
import net.imglib2.Point;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.BenchmarkAlgorithm;
import net.imglib2.algorithm.OutputAlgorithm;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import peakFitter.GaussianMaskFit.Endfit;
import preProcessing.GetLocalmaxmin;

public class SubpixelVelocityMSER extends BenchmarkAlgorithm
implements OutputAlgorithm<ArrayList<double[]>> {

	private static final String BASE_ERROR_MSG = "[SubpixelVelocity] ";
	private final RandomAccessibleInterval<FloatType> source;
	private final ArrayList<LabelledImg> imgs;
	private final ArrayList<double[]> PrevFrameparam;
	private final int ndims;
	private final int framenumber;
	private ArrayList<double[]> final_paramlist;
	private ArrayList<Staticproperties> startandendinframe;
	private final double[] psf;
	
	// LM solver iteration params
	final int maxiter = 500;
	final double lambda = 1e-3;
	final double termepsilon = 1e-1;
	//Mask fits iteration param
	final int iterations = 500;
	final double cutoffdistance = 10;
	final boolean halfgaussian = false;
	final double Intensityratio = 0.5;
	
	public  SubpixelVelocityMSER(final RandomAccessibleInterval<FloatType> source, 
			                      final ArrayList<LabelledImg> imgs,
			                       final ArrayList<double[]> PrevFrameparam,
			                       final double[] psf,
			                       final int framenumber) {
		
		this.source = source;
		this.imgs = imgs;
		this.PrevFrameparam = PrevFrameparam;
		this.psf = psf;
		this.framenumber = framenumber;
		this.ndims = source.numDimensions();
		
	}
	
	@Override
	public boolean checkInput() {
		if (source.numDimensions() > 2) {
			errorMessage = BASE_ERROR_MSG + " Can only operate on 1D, 2D, make slices of your stack . Got "
					+ source.numDimensions() + "D.";
			return false;
		}
		return true;
	}

	@Override
	public boolean process() {
		
		final_paramlist = new ArrayList<double[]>();
		startandendinframe = new ArrayList<Staticproperties>();
		
		for (int index = 0; index < PrevFrameparam.size(); ++index) {

			Point linepoint = new Point(ndims);
			linepoint.setPosition(
					new long[] { (long) PrevFrameparam.get(index)[0], (long) PrevFrameparam.get(index)[1] });
			
			
			
			 ArrayList<Integer> label = Getlabel(linepoint);
			 Set<Integer> multilabel = new HashSet<Integer>(label);
			 
			
			 Point secondlinepoint = new Point(ndims);
				secondlinepoint.setPosition(
						new long[] { (long) PrevFrameparam.get(index)[2], (long) PrevFrameparam.get(index)[3] });
				
				
				 
				 ArrayList<Integer> seclabel = Getlabel(secondlinepoint);
				 Set<Integer> secmultilabel = new HashSet<Integer>(seclabel);
				 
				secmultilabel.retainAll(multilabel);
				Iterator<Integer> iter = secmultilabel.iterator();
				while(iter.hasNext()){
				int currentlabel = iter.next();
				 
			
				 
				 
				 if (currentlabel != Integer.MIN_VALUE){
					 System.out.println(currentlabel);
			 double[] paramnextframe =Getfinaltrackparam(PrevFrameparam.get(index),
							currentlabel, psf, framenumber);
			 if (paramnextframe!=null){
			 final_paramlist.add(paramnextframe);
			 
			 
			 final double[] oldstartpoint = {PrevFrameparam.get(index)[0], PrevFrameparam.get(index)[1]};
			 
			 final double[] oldendpoint = {PrevFrameparam.get(index)[2], PrevFrameparam.get(index)[3]};
			 
			 final double[] newstartpoint = oldstartpoint;
			 
			 final double[] newendpoint = oldendpoint;
			
			 if (paramnextframe!=null){
			 for (int d = 0; d < ndims; ++d){
				 
				 newstartpoint[d] = paramnextframe[d];
				 newendpoint[d] = paramnextframe[d + ndims];
			 }
			 }
				
			 final double[] directionstart = {newstartpoint[0] - oldstartpoint[0] , newstartpoint[1] - oldstartpoint[1] };
			 
			 final double[] directionend = {newendpoint[0] - oldendpoint[0] , newendpoint[1] - oldendpoint[1] };
			 
			System.out.println("Frame:" + framenumber + " " +  "Fits :" + currentlabel + " "+ "StartX:" + paramnextframe[0] 
					+ " StartY:" + paramnextframe[1] + " " + "EndX:"
					+ paramnextframe[2] + "EndY: " + paramnextframe[3]);
			
			final Staticproperties edge = 
		   new Staticproperties(currentlabel, oldstartpoint, oldendpoint, newstartpoint, newendpoint, directionstart , directionend );
			

					startandendinframe.add(edge);	
				 }
				 }
				 // The line is missing in the current frame
				 else{
					 final double[] oldstartpoint = {PrevFrameparam.get(index)[0], PrevFrameparam.get(index)[1]};
					 
					 final double[] oldendpoint = {PrevFrameparam.get(index)[2], PrevFrameparam.get(index)[3]};
					 final double[] newstartpoint = oldstartpoint;
					 
					 final double[] newendpoint = oldendpoint;
					 final double[] directionstart = {newstartpoint[0] - oldstartpoint[0] , newstartpoint[1] - oldstartpoint[1] };
					 
					 final double[] directionend = {newendpoint[0] - oldendpoint[0] , newendpoint[1] - oldendpoint[1] };

						final Staticproperties edge = 
					   new Staticproperties(currentlabel, oldstartpoint, oldendpoint, newstartpoint, newendpoint, directionstart , directionend );
						
					 startandendinframe.add(edge);
				 }
		}
		}
		return false;
	}

	@Override
	public ArrayList<double[]> getResult() {
		return final_paramlist;
	} 
	
	public ArrayList<Staticproperties> getStateVectors() {
		return startandendinframe;
	} 
	
	private final double[] MakerepeatedLineguess(double[] iniparam, int label)  {
		long[] newposition = new long[ndims];
		double[] minVal = { Double.MAX_VALUE, Double.MAX_VALUE };
		double[] maxVal = { -Double.MIN_VALUE, -Double.MIN_VALUE };

		RandomAccessibleInterval<FloatType> currentimg = Boundingboxes.CurrentLabelImage(imgs, label);
				//imgs.get(label).Actualroiimg;
		final double[] cordone = { iniparam[0], iniparam[1] };
		final double[] cordtwo = { iniparam[2], iniparam[3] };

		double slope = (cordone[1] - cordtwo[1]) / (cordone[0] - cordtwo[0]);
		double intercept = cordone[1] - slope * cordone[0];
		double newintercept = intercept;


		final Cursor<FloatType> outcursor = Views.iterable(currentimg).localizingCursor();

		final double maxintensityline = GetLocalmaxmin.computeMaxIntensity(currentimg);

		while (outcursor.hasNext()) {

			outcursor.fwd();

			if (outcursor.get().get() / maxintensityline > Intensityratio) {
				outcursor.localize(newposition);

				long pointonline = (long) (outcursor.getLongPosition(1) - slope * outcursor.getLongPosition(0) - newintercept);

				// To get the min and max co-rodinates along the line so we
				// have starting points to
				// move on the line smoothly

				if (pointonline == 0) {
					for (int d = 0; d < ndims; ++d) {
						if (outcursor.getDoublePosition(d) <= minVal[d])
							minVal[d] = outcursor.getDoublePosition(d);

						if (outcursor.getDoublePosition(d) >= maxVal[d])
							maxVal[d] = outcursor.getDoublePosition(d);

					}

				}
			}
		}

		final double[] MinandMax = new double[2 * ndims + 3];

		if (slope >= 0) {
			for (int d = 0; d < ndims; ++d) {

				MinandMax[d] = minVal[d];
				MinandMax[d + ndims] = maxVal[d];
			}

		}

		if (slope < 0) {

			MinandMax[0] = minVal[0];
			MinandMax[1] = maxVal[1];
			MinandMax[2] = maxVal[0];
			MinandMax[3] = minVal[1];

		}

		MinandMax[2 * ndims] = iniparam[2 * ndims];
		MinandMax[2 * ndims + 1] = iniparam[2 * ndims + 1];
		MinandMax[2 * ndims + 2] = iniparam[2 * ndims + 2];
		System.out.println("Label: " + label + " " + "Initial guess: " + " StartX: " + MinandMax[0] + " StartY: "
				+ MinandMax[1] + " EndX: " + MinandMax[2] + " EndY: " + MinandMax[3]);

		
			for (int d = 0; d < ndims; ++d) {

				if (MinandMax[d] == Double.MAX_VALUE || MinandMax[d + ndims] == -Double.MIN_VALUE)
					return null;
				if (MinandMax[d] >= source.dimension(d) || MinandMax[d + ndims] >= source.dimension(d))
					return null;
				if (MinandMax[d] <= 0 || MinandMax[d + ndims] <= 0)
					return null;

			}
		

		
		return MinandMax;

	}
	public double[] Getfinaltrackparam(final double[] iniparam, final int label, final double[] psf, final int rate)  {

		if (iniparam == null || label == Integer.MIN_VALUE)
			return null;

		else {

			PointSampleList<FloatType> datalist = gatherfullData(label);
			final Cursor<FloatType> listcursor = datalist.localizingCursor();
			double[][] X = new double[(int) datalist.size()][ndims];
			double[] I = new double[(int) datalist.size()];
			int index = 0;
			while (listcursor.hasNext()) {
				listcursor.fwd();

				for (int d = 0; d < ndims; d++) {
					X[index][d] = listcursor.getDoublePosition(d);
				}

				I[index] = listcursor.get().getRealDouble();

				index++;
			}

			final double[] finalparamstart = MakerepeatedLineguess(iniparam, label);
			if (finalparamstart == null)
				return null;

			else {
				RandomAccessibleInterval<FloatType> currentimg = Boundingboxes.CurrentLabelImage(imgs, label);
						//imgs.get(label).Actualroiimg;

				final double[] fixed_param = new double[ndims];

				for (int d = 0; d < ndims; ++d) {

					fixed_param[d] = 1.0 / Math.pow(psf[d], 2);
				}

				final double[] inistartpos = { finalparamstart[0], finalparamstart[1] };
				final double[] iniendpos = { finalparamstart[2], finalparamstart[3] };

				double inicutoffdistance = Distance(inistartpos, iniendpos);

				// LM solver part
				if (inicutoffdistance > 2) {
					try {
						LevenbergMarquardtSolverLine.solve(X, finalparamstart, fixed_param, I, new GaussianLineds(), lambda,
								termepsilon, maxiter);
					} catch (Exception e) {
						e.printStackTrace();
					}

					final double[] startpos = { finalparamstart[0], finalparamstart[1] };
					final double[] endpos = { finalparamstart[2], finalparamstart[3] };
					// NaN protection: we prefer returning the crude estimate
					// than
					// NaN
					for (int j = 0; j < finalparamstart.length; j++) {
						if (Double.isNaN(finalparamstart[j]))
							finalparamstart[j] = iniparam[j];
					}

					final double LMdist = sqDistance(startpos, endpos);

					double[] returnparam = new double[2 * ndims + 3];

					final double maxintensityline = GetLocalmaxmin.computeMaxIntensity(currentimg);

					

					double newslope = (endpos[1] - startpos[1]) / (endpos[0] - startpos[0]);
					double newintercept = (endpos[1] - newslope * endpos[0]);
					double dx = finalparamstart[4] / Math.sqrt(1 + newslope * newslope);
					double dy = newslope * dx;
					double[] dxvector = { dx,  dy };

					double[] startfit = new double[ndims];
					double[] endfit = new double[ndims];

					
					

					try {
						startfit = peakFitter.GaussianMaskFit.sumofgaussianMaskFit(currentimg,  startpos.clone(),
								psf, iterations, dxvector, newslope, newintercept, maxintensityline,  halfgaussian,
								Endfit.Start, label);
					} catch (Exception e) {
						e.printStackTrace();
					}

					try {
						endfit = peakFitter.GaussianMaskFit.sumofgaussianMaskFit(currentimg,  endpos.clone(), psf,
								iterations, dxvector, newslope, newintercept, maxintensityline,  halfgaussian,
								Endfit.End, label);
					} catch (Exception e) {
						e.printStackTrace();
					}

					final double Maskdist = sqDistance(startfit, endfit);
					// If mask fits fail, return LM solver results


					for (int d = 0; d < ndims; ++d) {
						returnparam[d] = startfit[d];
						returnparam[ndims + d] = endfit[d];
					}

					
					if (Math.abs(Math.sqrt(Maskdist)) - Math.sqrt(LMdist) > cutoffdistance){
						if (Math.abs(startpos[0] - startfit[0]) >= cutoffdistance / 2 && Math.abs(startpos[1] - startfit[1]) >= cutoffdistance / 2
								|| Math.abs(endpos[0] - endfit[0]) >= cutoffdistance / 2 && Math.abs(endpos[1] - endfit[1]) >= cutoffdistance / 2 ){
							System.out.println("Mask fits fail, returning LM solver results!");
						
							for (int d = 0; d < ndims; ++d) {
							returnparam[d] = startpos[d];
							returnparam[ndims + d] = endpos[d];
						}
						}
					
						if (Math.abs(startpos[0] - startfit[0]) >= cutoffdistance || Math.abs(startpos[1] - startfit[1]) >= cutoffdistance 
								|| Math.abs(endpos[0] - endfit[0]) >= cutoffdistance  || Math.abs(endpos[1] - endfit[1]) >= cutoffdistance  ){
							System.out.println("Mask fits fail, returning LM solver results!");
							for (int d = 0; d < ndims; ++d) {
								returnparam[d] = startpos[d];
								returnparam[ndims + d] = endpos[d];
							}
							
						}
					
					
					}
					
					
					
					for (int d = 0; d < ndims; ++d) {
						if (Double.isNaN(startfit[d]) || Double.isNaN(endfit[d])) {
							System.out.println("Mask fits fail, returning LM solver results!");
							returnparam[d] = startpos[d];
							returnparam[ndims + d] = endpos[d];

						}

					}

				

					returnparam[2 * ndims] = finalparamstart[4];
					returnparam[2 * ndims + 1] = finalparamstart[5];
					returnparam[2 * ndims + 2] = finalparamstart[6];
					return returnparam;
				} else
					return null;
			}

		}

	}
	
	

	private PointSampleList<FloatType> gatherfullData(final int label) {
		final PointSampleList<FloatType> datalist = new PointSampleList<FloatType>(ndims);

		RandomAccessibleInterval<FloatType> currentimg = Boundingboxes.CurrentLabelImage(imgs, label);

		
		Cursor<FloatType> localcursor = Views.iterable(currentimg).localizingCursor();

		while (localcursor.hasNext()) {
			localcursor.fwd();
			
			Point newpoint = new Point(localcursor);
			datalist.add(newpoint, localcursor.get().copy());

		}

		return datalist;
	}
	public ArrayList<Integer> Getlabel(final Point linepoint) {

		final int x = linepoint.getIntPosition(0);
		final int y = linepoint.getIntPosition(1);
		ArrayList<Integer> currentlabel = new ArrayList<Integer>();
		for (int index = 0; index < imgs.size(); ++index){
			
			EllipseRoi ellipse = imgs.get(index).roi;
			
			if (ellipse.contains(x, y)){
				
				currentlabel.add(index);
			}
			
		}
		
	

		return currentlabel;
	}
	public double Distance(final double[] cordone, final double[] cordtwo) {

		double distance = 0;

		for (int d = 0; d < ndims; ++d) {

			distance += Math.pow((cordone[d] - cordtwo[d]), 2);

		}
		return Math.sqrt(distance);
	}

	public double sqDistance(final double[] cordone, final double[] cordtwo) {

		double distance = 0;

		for (int d = 0; d < ndims; ++d) {

			distance += Math.pow((cordone[d] - cordtwo[d]), 2);

		}
		return (distance);
	}
}

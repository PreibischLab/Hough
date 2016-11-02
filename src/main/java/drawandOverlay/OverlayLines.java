package drawandOverlay;

import java.util.ArrayList;
import houghandWatershed.TransformCordinates;
import ij.gui.EllipseRoi;
import labeledObjects.LabelledImg;
import labeledObjects.Lineobjects;
import labeledObjects.Simpleobject;
import net.imglib2.Cursor;
import net.imglib2.Point;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.localextrema.RefinedPeak;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import preProcessing.GetLocalmaxmin;
import preProcessing.GetLocalmaxmin.IntensityType;

public class OverlayLines {
	


	

	public static ArrayList<RefinedPeak<Point>> ReducedList(RandomAccessibleInterval<FloatType> inputimg,
			ArrayList<RefinedPeak<Point>> SubpixelMinlist, double[] sizes, double[] min, double[] max) {

		RandomAccessibleInterval<FloatType> imgout = new ArrayImgFactory<FloatType>().create(inputimg, new FloatType());
		double[] points = new double[imgout.numDimensions()];
		int maxcount = 0;
		int maxindex = 0;
		
		ArrayList<RefinedPeak<Point>> MainMinlist = new ArrayList<RefinedPeak<Point>>(inputimg.numDimensions());
		for (int index = 0; index < SubpixelMinlist.size(); ++index) {
			points = TransformCordinates.transformfwd(new double[] { SubpixelMinlist.get(index).getDoublePosition(0),
					SubpixelMinlist.get(index).getDoublePosition(1) }, sizes, min, max);

			double slope = -1.0 / Math.tan(Math.toRadians(points[0]));
			double intercept = points[1] / Math.sin(Math.toRadians(points[0]));
			
			PushCurves.Drawexactline(imgout, slope, intercept, IntensityType.Gaussian);

			RandomAccess<FloatType> inran = inputimg.randomAccess();
			Cursor<FloatType> outcursor = Views.iterable(imgout).localizingCursor();
			

			int count = 0;
			while (outcursor.hasNext()) {
				outcursor.fwd();

				if (outcursor.get().get() > 0) {
					inran.setPosition(outcursor);

					if (inran.get().get() > 0)
						count++;

					
					if (count > maxcount) {
						maxcount = count;
						maxindex = index;

					}

				}

			}
		}
		
		
		if (maxcount > 0) {
			MainMinlist.add(SubpixelMinlist.get(maxindex));
			
		}
		
		
		
		
		return MainMinlist;
	}

	public static ArrayList<double[]> GetRhoTheta(ArrayList<RefinedPeak<Point>> MainMinlist, double[] sizes, double[] min,
			double[] max) {

		ArrayList<double[]> points = new ArrayList<double[]>(); //[sizes.length];
		for (int index = 0; index < MainMinlist.size(); ++index) {

		final double[]	point = TransformCordinates.transformfwd(new double[] { MainMinlist.get(index).getDoublePosition(0),
					MainMinlist.get(index).getDoublePosition(1) }, sizes, min, max);
		
		points.add(point);
		}
		return points;
	}

	
	public static void Getmserlines(RandomAccessibleInterval<FloatType> imgout,
			ArrayList<LabelledImg> imgslist,ArrayList<Simpleobject> lineobject){
		
		for (int index = 0; index < imgslist.size(); ++index){
			
			final int label = imgslist.get(index).label;
			final double slope = imgslist.get(index).slopeandintercept[0];
			final double intercept = imgslist.get(index).slopeandintercept[1];
			final EllipseRoi ellipse = imgslist.get(index).roi;
			final Simpleobject simpleobj = new Simpleobject(label, slope, intercept);
			lineobject.add(simpleobj);
			PushCurves.DrawRoiline(imgout, ellipse, slope, intercept);
			
		}
		
		
	}
	public static void GetAlllines(
			RandomAccessibleInterval<FloatType> imgout,
			
			RandomAccessibleInterval<IntType> intimg, 
			ArrayList<Lineobjects> linelist,
			ArrayList<Simpleobject> lineobject,
			final long radius) {

		for (int index = 0; index < linelist.size(); ++index) {

			final int label = linelist.get(index).Label;
			
			ArrayList<double[]> rhothetalist = linelist.get(index).rhotheta;
			
			for (int arrayindex = 0; arrayindex < rhothetalist.size(); ++arrayindex){
			final double rho = rhothetalist.get(arrayindex)[1];
			final double theta = rhothetalist.get(arrayindex)[0];
			
			double slope = -1.0 / (Math.tan(Math.toRadians(theta)));
			double intercept = rho / Math.sin(Math.toRadians(theta));

			final Simpleobject simpleobj = new Simpleobject(label, slope, intercept);
			lineobject.add(simpleobj);
			
			PushCurves.DrawTruncatedline(imgout, intimg, slope, intercept, label);

			}
		}
		
	}
	public static void GetCurrentlines(
			RandomAccessibleInterval<FloatType> imgout,
			RandomAccessibleInterval<FloatType> maximgout,
			Img<IntType> intimg, 
			ArrayList<Lineobjects> linelist,
			int currentlabel) {

		for (int index = 0; index < linelist.size(); ++index) {

			final int label = linelist.get(index).Label;
			
		ArrayList<double[]> rhothetalist = linelist.get(index).rhotheta;
		
		for (int arrayindex = 0; arrayindex < rhothetalist.size(); ++arrayindex){
			final double rho = linelist.get(index).rhotheta.get(arrayindex)[1];
			final double theta = linelist.get(index).rhotheta.get(arrayindex)[0];
		
		if (label == currentlabel){
			
			double slope = -1.0 / Math.tan(Math.toRadians(theta));
			double intercept = rho / Math.sin(Math.toRadians(theta));

			PushCurves.Drawexactline(imgout,intimg, slope,intercept, label);
			maximgout = GetLocalmaxmin.FindandDisplayLocalMaxima(imgout,
					IntensityType.Original, new double[]{1,1});
			
		}
		}
		}
	}

}
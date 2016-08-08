package peakFitter;

import java.util.ArrayList;

import com.sun.tools.javac.util.Pair;

import drawandOverlay.AddGaussian;
import drawandOverlay.PushCurves;
import houghandWatershed.Finalfunction;
import ij.plugin.HyperStackReducer;
import labeledObjects.Finalobject;
import labeledObjects.Indexedlength;
import labeledObjects.LabelMax;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.IterableInterval;
import net.imglib2.Localizable;
import net.imglib2.Point;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealCursor;
import net.imglib2.RealLocalizable;
import net.imglib2.RealPoint;
import net.imglib2.RealPointSampleList;
import net.imglib2.RealRandomAccess;
import net.imglib2.RealRandomAccessibleRealInterval;
import net.imglib2.algorithm.neighborhood.HyperSphereNeighborhood;
import net.imglib2.algorithm.neighborhood.RectangleShape;
import net.imglib2.algorithm.neighborhood.SquareStrelTest;
import net.imglib2.algorithm.region.hypersphere.HyperSphere;
import net.imglib2.algorithm.region.hypersphere.HyperSphereCursor;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.roi.RectangleRegionOfInterest;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Intervals;
import net.imglib2.view.Views;
import peakFitter.GaussianMaskFit.Endfit;
import preProcessing.GetLocalmaxmin;
import preProcessing.GetLocalmaxmin.IntensityType;

public class LengthDetection {

	private final RandomAccessibleInterval<FloatType> inputimg;
	private final RandomAccessibleInterval<IntType> intimg;
	private final int ndims;

	public LengthDetection(RandomAccessibleInterval<FloatType> inputimg, RandomAccessibleInterval<IntType> intimg) {

		this.inputimg = inputimg;
		this.intimg = intimg;
		this.ndims = inputimg.numDimensions();
		assert inputimg.numDimensions() == intimg.numDimensions();

	}
	
	

	
	public double Distance(final double[] cordone, final double[] cordtwo) {

		double distance = 0;

		for (int d = 0; d < ndims; ++d) {

			distance += Math.pow((cordone[d] - cordtwo[d]), 2);

		}
		return Math.sqrt(distance);
	}

	public ArrayList<Finalobject> Updateslopeandintercept(ArrayList<Finalobject> finalparam) {

		int Maxlabel = houghandWatershed.PerformWatershedding.GetMaxlabelsseeded(intimg);
		final ArrayList<Finalobject> updateparamlist = new ArrayList<Finalobject>();
		for (int label = 1; label < Maxlabel - 1; ++label) {

			ArrayList<RealPoint> centroidlist = new ArrayList<RealPoint>();

			int labelindex = 0;

			for (int index = 0; index < finalparam.size(); ++index) {

				if (finalparam.get(index).Label == label) {

					labelindex = index;
					centroidlist.add(finalparam.get(index).centroid);

				}

			}

			double newslope = 0;
			double newintercept = 0;

			final double[] pointone = new double[ndims];
			final double[] pointtwo = new double[ndims];

			if (centroidlist.size() > 0) {
				centroidlist.get(0).localize(pointone);
				centroidlist.get(centroidlist.size() - 1).localize(pointtwo);

				if (pointtwo[0] != pointone[0] ){
				newslope = (pointtwo[1] - pointone[1]) / (pointtwo[0] - pointone[0]);
				newintercept = pointtwo[1] - newslope * pointtwo[0];
				}
				
				else{
					centroidlist.get(centroidlist.size()/2).localize(pointone);
					centroidlist.get(3*centroidlist.size()/4).localize(pointtwo);
					newslope = (pointtwo[1] - pointone[1]) / (pointtwo[0] - pointone[0]);
					newintercept = pointtwo[1] - newslope * pointtwo[0];
				}
				

				if (Double.isNaN(newslope) ){
					newslope = finalparam.get(labelindex).slope ;
					newintercept = finalparam.get(labelindex).intercept;
				}
				
					
				
			}

			System.out.println("old: " + finalparam.get(labelindex).slope + " " + finalparam.get(labelindex).intercept);

			System.out.println("new: " + newslope + " " + newintercept);

			for (int index = 0; index < finalparam.size(); ++index) {

				if (finalparam.get(index).Label == label) {
					Finalobject update = new Finalobject(label, finalparam.get(index).centroid,
							finalparam.get(index).Intensity, finalparam.get(index).sigmaX, finalparam.get(index).sigmaY,
							newslope, newintercept);

					updateparamlist.add(update);

				}
			}

		}

		return updateparamlist;
	}

	public ArrayList<Finalobject> Removepoints(ArrayList<Finalobject> finalparam, ArrayList<LabelMax> labelmaxlist) {

		int Maxlabel = houghandWatershed.PerformWatershedding.GetMaxlabelsseeded(intimg);

		final ArrayList<Finalobject> correctparamlist = new ArrayList<Finalobject>();
		for (int label = 1; label < Maxlabel - 1; ++label) {
			
			float[] pos = new float[ndims];
			double maxintensity = Double.MIN_VALUE;

			for (int index = 0; index < finalparam.size(); ++index) {

				if (finalparam.get(index).Label == label) {

					
					if (finalparam.get(index).Intensity > maxintensity) {

						maxintensity = finalparam.get(index).Intensity;

					}

				}

			}

			final LabelMax labelmax = new LabelMax(label, maxintensity);
			labelmaxlist.add(labelmax);
			System.out.println("Label :" + label + " " + maxintensity);

			for (int listindex = 0; listindex < finalparam.size(); ++listindex) {

				if (finalparam.get(listindex).Label == label) {

					if (finalparam.get(listindex).Intensity / maxintensity > 0.5   )

						correctparamlist.add(finalparam.get(listindex));
				}

			}
		}

		return correctparamlist;

	}

	public void Returnlengths(ArrayList<Finalobject> finalparam, ArrayList<Indexedlength> finallength,
			ArrayList<LabelMax> labelmaxlist, double[] sigma) throws Exception {

		int Maxlabel = houghandWatershed.PerformWatershedding.GetMaxlabelsseeded(intimg);

		for (int label = 1; label < Maxlabel - 1; ++label) {
			final double[] minVal = { Double.MAX_VALUE, Double.MAX_VALUE };
			final double[] maxVal = { Double.MIN_VALUE, Double.MIN_VALUE };
			double length = 0;
			double lengthpre = 0;

			double[] minposition = new double[ndims];
			double[] maxposition = new double[ndims];

			double[] startpos = new double[ndims];
			double[] endpos = new double[ndims];
			double slope = 0;
			double intercept = 0;
			int labelindex = 0;
			for (int index = 0; index < finalparam.size(); ++index) {

				if (finalparam.get(index).Label == label) {
					labelindex = index;

					for (int d = 0; d < ndims; ++d) {
						minposition[d] = finalparam.get(index).centroid.getDoublePosition(d);
						maxposition[d] = finalparam.get(index).centroid.getDoublePosition(d);

						
						if (minposition[d] <= minVal[d]) {
							minVal[d] = minposition[d];
						}
						if (maxposition[d] >= maxVal[d]) {
							maxVal[d] = maxposition[d];
						}

					}

					if (finalparam.get(index).slope >= 0) {

						for (int d = 0; d < ndims; ++d) {

							startpos[d] = minVal[d];
							endpos[d] = maxVal[d];

						}
					}

					if (finalparam.get(index).slope < 0) {

						startpos[0] = minVal[0];
						startpos[1] = maxVal[1];

						endpos[0] = maxVal[0];
						endpos[1] = minVal[1];

					}

				}

			}

			slope = finalparam.get(labelindex).slope;
			intercept = finalparam.get(labelindex).intercept;
			lengthpre = Distance(startpos, endpos);

			System.out.println("Label :" + label + " " + " StartX: " + startpos[0] + " " + " StartY: " + startpos[1]
					+ "EndX :" + endpos[0] + " " + "EndY :" + endpos[1] + " " + "Length: " + lengthpre);

			final int iterations = 5000;
			double maxintensity = 0;

			for (int listlabelindex = 0; listlabelindex < labelmaxlist.size(); ++listlabelindex) {

				if (labelmaxlist.get(listlabelindex).Label == label)

					maxintensity = labelmaxlist.get(listlabelindex).maxIntensity;

			}

			final double[] newsigma = { sigma[0], sigma[1] };
			final double[] startfit = peakFitter.GaussianMaskFit.gaussianMaskFit(inputimg, intimg, startpos, newsigma,
					iterations, maxintensity,1.0,  slope, intercept, Endfit.Start);
			final double[] endfit = peakFitter.GaussianMaskFit.gaussianMaskFit(inputimg, intimg, endpos, newsigma,
					iterations, maxintensity, 1.0, slope, intercept, Endfit.End);

			
			

			length = Distance(startfit, endfit);
			final Indexedlength currentlength = new Indexedlength(label, length, startfit, endfit, slope, intercept);

			finallength.add(currentlength);

			System.out.println(
					"New:" + "Label :" + label + " " + " StartX: " + startfit[0] + " " + " StartY: " + startfit[1] + " "
							+ "EndX :" + endfit[0] + " " + "EndY :" + endfit[1] + "  " + "Length: " + length);

		}

	}

}

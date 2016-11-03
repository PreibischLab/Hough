package getRoi;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import houghandWatershed.Boundingboxes;
import houghandWatershed.HoughTransform2D;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.gui.EllipseRoi;
import ij.gui.Overlay;
import labeledObjects.LabelledImg;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.BenchmarkAlgorithm;
import net.imglib2.algorithm.OutputAlgorithm;
import net.imglib2.algorithm.componenttree.mser.Mser;
import net.imglib2.algorithm.componenttree.mser.MserTree;
import net.imglib2.img.ImagePlusAdapter;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.integer.UnsignedByteType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Util;
import net.imglib2.view.Views;
import peakFitter.SortListbyproperty;

public class RoiforMSER extends BenchmarkAlgorithm
		implements OutputAlgorithm <ArrayList<LabelledImg>> {

	private static final String BASE_ERROR_MSG = "[RoiforMSER] ";
	private final RandomAccessibleInterval<FloatType> source;
	private final RandomAccessibleInterval<FloatType> Actualsource;
	private final double delta;
	private final long minSize;
	private final long maxSize;
	private final double maxVar;
	private final double minDIversity;
	private final boolean darktoBright;
	private Overlay ov;
	
	private int Roiindex;
	private ArrayList<LabelledImg> imgs;

	public RoiforMSER(final RandomAccessibleInterval<FloatType> source, final RandomAccessibleInterval<FloatType> Actualsource,
			final double delta, final long minSize,
			final long maxSize, final double maxVar, final double minDiversity, final boolean darktoBright) {

		this.source = source;
		this.Actualsource = Actualsource;
		this.delta = delta;
		this.minSize = minSize;
		this.maxSize = maxSize;
		this.maxVar = maxVar;
		this.minDIversity = minDiversity;
		this.darktoBright = darktoBright;
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

		final FloatType type = source.randomAccess().get().createVariable();
		
		imgs = new ArrayList<LabelledImg>();

		ov = new Overlay();
		ArrayList<double[]> ellipselist = new ArrayList<double[]>();
		ArrayList<double[]> meanandcovlist = new ArrayList<double[]>();
		final Img<UnsignedByteType> newimg;

		try
		{
		ImageJFunctions.wrap(source, "curr");
		final ImagePlus currentimp = IJ.getImage();
		IJ.run("8-bit");

		newimg = ImagePlusAdapter.wrapByte(currentimp);

		}
		catch ( final Exception e )
		{
			e.printStackTrace();
			return false;
		}
		
		MserTree<UnsignedByteType> newtree = MserTree.buildMserTree(newimg, delta, minSize, maxSize, maxVar,
				minDIversity, darktoBright);
		final HashSet<Mser<UnsignedByteType>> rootset = newtree.roots();
		
		
		final Iterator<Mser<UnsignedByteType>> rootsetiterator = rootset.iterator();
		
		
		
		
		while (rootsetiterator.hasNext()) {

			Mser<UnsignedByteType> rootmser = rootsetiterator.next();

			if (rootmser.size() > 0) {

				final double[] meanandcov = { rootmser.mean()[0], rootmser.mean()[1], rootmser.cov()[0],
						rootmser.cov()[1], rootmser.cov()[2] };
				meanandcovlist.add(meanandcov);
				ellipselist.add(meanandcov);

			}
		}
		
		// We do this so the ROI remains attached the the same label and is not changed if the program is run again
	SortListbyproperty.sortpointList(ellipselist);
		
			for (int index = 0; index < ellipselist.size(); ++index) {
				Roiindex = index;
				
				final ImgFactory<FloatType> factory = Util.getArrayOrCellImgFactory(source, type);
				RandomAccessibleInterval<FloatType>  Roiimg = factory.create(source, type);
				RandomAccessibleInterval<FloatType>  ActualRoiimg = factory.create(Actualsource, type);
				
				final double[] mean = { ellipselist.get(index)[0], ellipselist.get(index)[1] };
				final double[] covar = { ellipselist.get(index)[2], ellipselist.get(index)[3],
						ellipselist.get(index)[4] };
				final EllipseRoi ellipseroi = createEllipse(mean, covar, 3);
				ellipseroi.setStrokeColor(Color.green);
				ov.add(ellipseroi);

				
				
				

				Cursor<FloatType> sourcecursor = Views.iterable(source).localizingCursor();
				RandomAccess<FloatType> ranac = Roiimg.randomAccess();
				while (sourcecursor.hasNext()) {

					sourcecursor.fwd();

					final int x = sourcecursor.getIntPosition(0);
					final int y = sourcecursor.getIntPosition(1);
					ranac.setPosition(sourcecursor);
					if (ellipseroi.contains(x, y)) {
						
						ranac.get().set(sourcecursor.get());

					}
					

				}
				Cursor<FloatType> Actualsourcecursor = Views.iterable(Actualsource).localizingCursor();
				RandomAccess<FloatType> Actualranac = ActualRoiimg.randomAccess();
				while (Actualsourcecursor.hasNext()) {

					Actualsourcecursor.fwd();

					final int x = Actualsourcecursor.getIntPosition(0);
					final int y = Actualsourcecursor.getIntPosition(1);
					Actualranac.setPosition(Actualsourcecursor);
					if (ellipseroi.contains(x, y)) {
						
						Actualranac.get().set(Actualsourcecursor.get());

					}
					

				}
				
				// Obtain the slope and intercept of the line by Hough Transform (slow but very accurate)
				//final double[] slopeandintercept = LargestEigenvector(mean, covar, Roiimg, ellipseroi, Roiindex);
				
				// Obtain the slope and intercept of the line by obtaining the major axis of the ellipse (super fast but could be inaccurate)
				final double[] slopeandintercept = LargestEigenvector(mean, covar);
				
				
				LabelledImg currentimg = new LabelledImg(Roiindex, Roiimg, ActualRoiimg, ellipseroi, slopeandintercept, mean, covar);
				if(Math.abs(slopeandintercept[0])!=Double.POSITIVE_INFINITY || Math.abs(slopeandintercept[1])!=Double.POSITIVE_INFINITY  )
				imgs.add(currentimg);
				
				
			}

		

		return true;
	}

	@Override
	public ArrayList<LabelledImg> getResult() {
		
		return imgs;
	}
	
    public Overlay getOverlay() {
		
		return ov;
	}
	

	/**
	 * 2D correlated Gaussian
	 * 
	 * @param mean
	 *            (x,y) components of mean vector
	 * @param cov
	 *            (xx, xy, yy) components of covariance matrix
	 * @return ImageJ roi
	 */
	public static EllipseRoi createEllipse(final double[] mean, final double[] cov, final double nsigmas) {
		final double a = cov[0];
		final double b = cov[1];
		final double c = cov[2];
		final double d = Math.sqrt(a * a + 4 * b * b - 2 * a * c + c * c);
		final double scale1 = Math.sqrt(0.5 * (a + c + d)) * nsigmas;
		final double scale2 = Math.sqrt(0.5 * (a + c - d)) * nsigmas;
		final double theta = 0.5 * Math.atan2((2 * b), (a - c));
		final double x = mean[0];
		final double y = mean[1];
		final double dx = scale1 * Math.cos(theta);
		final double dy = scale1 * Math.sin(theta);
		final EllipseRoi ellipse = new EllipseRoi(x - dx, y - dy, x + dx, y + dy, scale2 / scale1);
		return ellipse;
	}

	
	/**
	 * Returns the slope and the intercept of the line passing through the major axis of the ellipse
	 * 
	 * 
	 *@param mean
	 *            (x,y) components of mean vector
	 * @param cov
	 *            (xx, xy, yy) components of covariance matrix
	 * @return slope and intercept of the line along the major axis
	 */
	public  double[] LargestEigenvector( final double[] mean, final double[] cov){
		
		final double a = cov[0];
		final double b = cov[1];
		final double c = cov[2];
		final double d = Math.sqrt(a * a + 4 * b * b - 2 * a * c + c * c);
		final double[] eigenvector1 = {2 * b, c - a + d};
		double[] LargerVec = new double[eigenvector1.length];

		LargerVec =  eigenvector1;
		
        final double slope = LargerVec[1] / LargerVec[0];
        final double intercept = mean[1] - mean[0] * slope;
       
        if (LargerVec[0]!= 0){
        double[] pair = {slope, intercept};
        return pair;
        }
        else{
        	double[] pair = {Double.MAX_VALUE, Double.MAX_VALUE};
        return pair;
        }
		
	}
	
	/**
	 * Returns the slope and the intercept of the line by doing Hough Transform
	 * 
	 * 
	 *@param mean
	 *            (x,y) components of mean vector
	 * @param cov
	 *            (xx, xy, yy) components of covariance matrix
	 * @return slope and intercept of the line along the major axis
	 */
	public  double[] LargestEigenvector( final double[] mean, final double[] cov, RandomAccessibleInterval<FloatType> currentimg, final EllipseRoi roi, final int label){
		System.out.println("Doing Hough Transform in Label: " + label);
		int n = currentimg.numDimensions();
		long[] position = new long[n];
		long[] minVal = { Long.MAX_VALUE, Long.MAX_VALUE };
		long[] maxVal = { Long.MIN_VALUE, Long.MIN_VALUE };
		
		Cursor<FloatType> localcursor = Views.iterable(currentimg).localizingCursor();

		while (localcursor.hasNext()) {
			localcursor.fwd();
			int x = localcursor.getIntPosition(0);
			int y = localcursor.getIntPosition(1);
			if (roi.contains(x, y)){

				localcursor.localize(position);
				for (int d = 0; d < n; ++d) {
					if (position[d] < minVal[d]) {
						minVal[d] = position[d];
					}
					if (position[d] > maxVal[d]) {
						maxVal[d] = position[d];
					}

				}
				
			}
		}
		
		FinalInterval interval = new FinalInterval(minVal, maxVal);
		RandomAccessibleInterval<FloatType> currentimgsmall = Views.interval(currentimg, interval);
	
		
		HoughTransform2D Houghonly = new HoughTransform2D(currentimgsmall, source);
		Houghonly.checkInput();
		Houghonly.process();
		double[] pair = Houghonly.getResult();
		
		
        
        return pair;
		
	}
	
	
	
}
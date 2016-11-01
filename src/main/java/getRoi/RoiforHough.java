package getRoi;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;

import com.sun.tools.javac.util.Pair;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.gui.EllipseRoi;
import labeledObjects.LabelledImg;
import labeledObjects.LabelledRoi;
import net.imglib2.Cursor;
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

public class RoiforHough extends BenchmarkAlgorithm
		implements OutputAlgorithm <ArrayList<LabelledImg>> {

	private static final String BASE_ERROR_MSG = "[RoiforHough] ";
	private final RandomAccessibleInterval<FloatType> source;
	private final double delta;
	private final long minSize;
	private final long maxSize;
	private final double maxVar;
	private final double minDIversity;
	private final boolean darktoBright;
	
	private int Roiindex;
	private ArrayList<LabelledRoi> rois;
	private ArrayList<LabelledImg> imgs;

	public RoiforHough(final RandomAccessibleInterval<FloatType> source, final double delta, final long minSize,
			final long maxSize, final double maxVar, final double minDiversity, final boolean darktoBright) {

		this.source = source;
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
		
		rois = new ArrayList<LabelledRoi>();
		imgs = new ArrayList<LabelledImg>();

		ArrayList<double[]> ellipselist = new ArrayList<double[]>();
		final Img<UnsignedByteType> newimg;

		try
		{
		new ImageJ();
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

				ellipselist.add(meanandcov);

			}
		}

			for (int index = 0; index < ellipselist.size(); ++index) {

				
				final ImgFactory<FloatType> factory = Util.getArrayOrCellImgFactory(source, type);
				RandomAccessibleInterval<FloatType>  Roiimg = factory.create(source, type);
				
				
				final double[] mean = { ellipselist.get(index)[0], ellipselist.get(index)[1] };
				final double[] covar = { ellipselist.get(index)[2], ellipselist.get(index)[3],
						ellipselist.get(index)[4] };
				EllipseRoi ellipseroi = createEllipse(mean, covar, 3);
				LabelledRoi currentroi = new LabelledRoi(index, ellipseroi);

				
				final double[] meanandintercept = LargestEigenvector(mean, covar);
				Roiindex = index;
				rois.add(currentroi);

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
				
				LabelledImg currentimg = new LabelledImg(Roiindex, Roiimg, meanandintercept);
				imgs.add(currentimg);
				
				
			}

		

		return true;
	}

	@Override
	public ArrayList<LabelledImg> getResult() {
		
		return imgs;
	}
	
	public ArrayList<LabelledRoi> getLabelledRoilist() {

		return rois;

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
	public static double[] LargestEigenvector( final double[] mean, final double[] cov){
		
		final double a = cov[0];
		final double b = cov[1];
		final double c = cov[2];
		final double d = Math.sqrt(a * a + 4 * b * b - 2 * a * c + c * c);
		final double[] eigenvector1 = {2 * b, c - a + d};
		final double[] eigenvector2 = {2 * b, c - a - d};
		final double mageigenvec1 = eigenvector1[0] * eigenvector1[0] + eigenvector1[1] * eigenvector1[1];
		final double mageigenvec2 = eigenvector2[0] * eigenvector2[0] + eigenvector2[1] * eigenvector2[1];
		double[] LargerVec = new double[eigenvector1.length];
		final double locationdiff = mageigenvec2 - mageigenvec1;
		final boolean minVec = locationdiff > 0;
        LargerVec= minVec ? eigenvector2 : eigenvector1;
		
        final double slope = LargerVec[1] / LargerVec[0];
        final double intercept = mean[1] - mean[0] * slope;
        
        double[] pair = {slope, intercept};
        
        return pair;
		
	}
	
	
}
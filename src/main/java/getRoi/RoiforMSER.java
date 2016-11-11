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

import com.sun.tools.javac.util.Pair;

import houghandWatershed.Boundingboxes;
import houghandWatershed.HoughTransformandMser;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.gui.EllipseRoi;
import ij.gui.Overlay;
import labeledObjects.LabelledImg;
import mserMethods.GetDelta;
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
import preProcessing.GetLocalmaxmin;

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
	private final int minlength;
	private final boolean darktoBright;
	private final boolean doHough;
	private Overlay ov;
	private final int ndims;
	private int Roiindex;
	private ArrayList<LabelledImg> imgs;

	public RoiforMSER(final RandomAccessibleInterval<FloatType> source, final RandomAccessibleInterval<FloatType> Actualsource,
			final double delta, final long minSize, 
			final long maxSize, final double maxVar, final double minDiversity, final int minlength,  final boolean darktoBright, final boolean doHough) {

		this.source = source;
		this.Actualsource = Actualsource;
		this.delta = delta;
		this.minSize = minSize;
		this.maxSize = maxSize;
		this.maxVar = maxVar;
		this.minDIversity = minDiversity;
		this.darktoBright = darktoBright;
		this.minlength = minlength;
		this.doHough = doHough;
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
		int count = 0;
			for (int index = 0; index < ellipselist.size(); ++index) {
				
				
				final ImgFactory<FloatType> factory = Util.getArrayOrCellImgFactory(source, type);
				RandomAccessibleInterval<FloatType>  Roiimg = factory.create(source, type);
				RandomAccessibleInterval<FloatType>  ActualRoiimg = factory.create(Actualsource, type);
				
				final double[] mean = { ellipselist.get(index)[0], ellipselist.get(index)[1] };
				final double[] covar = { ellipselist.get(index)[2], ellipselist.get(index)[3],
						ellipselist.get(index)[4] };
				final EllipseRoi ellipseroi = GetDelta.createEllipse(mean, covar, 3);
				
	    		final double perimeter = ellipseroi.getLength();
	    		final double smalleigenvalue = SmallerEigenvalue(mean, covar);
	    		if (perimeter > 4 * Math.PI * minlength && smalleigenvalue < 30){
	    			
	    			Roiindex = count;
	    			count++;
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
				
				double[] slopeandintercept = new double[ndims + 1];
				
				
				
				
				// Obtain the slope and intercept of the line by obtaining the major axis of the ellipse (super fast and accurate)
				if (doHough){
					
					HoughTransformandMser viaHough = new HoughTransformandMser(Roiindex, Roiimg, ellipseroi, minlength);
					viaHough.checkInput();
					viaHough.process();
					slopeandintercept = viaHough.getResult();
					
				}
				
				else{
					
				
				slopeandintercept = LargestEigenvector(mean, covar);
				
				}
				
				LabelledImg currentimg = new LabelledImg(Roiindex, Roiimg, ActualRoiimg, ellipseroi, slopeandintercept, mean, covar);
				
				if(slopeandintercept!=null  )
				imgs.add(currentimg);
				
				}
				
			}

		

		return true;
	}

	@Override
	public ArrayList<LabelledImg> getResult() {
		
		
		//ArrayList<LabelledImg> reducedimgs = Overlappingregions(imgs);
		return imgs;
	}
	
	public ArrayList<LabelledImg> getAllResult() {
		
		
		return imgs;
	}
	
    public Overlay getOverlay() {
		
		return ov;
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
		
		// For inifinite slope lines support is provided
		final double a = cov[0];
		final double b = cov[1];
		final double c = cov[2];
		final double d = Math.sqrt(a * a + 4 * b * b - 2 * a * c + c * c);
		final double[] eigenvector1 = {2 * b, c - a + d};
		double[] LargerVec = new double[eigenvector1.length + 1];

		LargerVec =  eigenvector1;
		
        final double slope = LargerVec[1] / (LargerVec[0] );
        final double intercept = mean[1] - mean[0] * slope;
       
        if (Math.abs(slope) != Double.POSITIVE_INFINITY){
        double[] pair = {slope, intercept, Double.MAX_VALUE};
        return pair;
      
        }
        
        else {
        	
        	double[] prependicular = {Double.MAX_VALUE, Double.MAX_VALUE, mean[0]};
        	return prependicular;
        	}
        	 
       
		
	}
	
	/**
	 * Returns the smallest eigenvalue of the ellipse
	 * 
	 * 
	 *@param mean
	 *            (x,y) components of mean vector
	 * @param cov
	 *            (xx, xy, yy) components of covariance matrix
	 * @return slope and intercept of the line along the major axis
	 */
	public  double SmallerEigenvalue( final double[] mean, final double[] cov){
		
		// For inifinite slope lines support is provided
		final double a = cov[0];
		final double b = cov[1];
		final double c = cov[2];
		final double d = Math.sqrt(a * a + 4 * b * b - 2 * a * c + c * c);

		
        final double smalleigenvalue = (a + c - d) / 2;
       
        
        	
        	return smalleigenvalue;
        	
        	 
       
		
	}
	
	
	
	
	
}
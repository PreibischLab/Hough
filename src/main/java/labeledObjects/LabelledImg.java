package labeledObjects;

import ij.gui.EllipseRoi;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.numeric.real.FloatType;

public class LabelledImg {

	
	public final int label;
	public final RandomAccessibleInterval<FloatType> roiimg;
	public final RandomAccessibleInterval<FloatType> Actualroiimg;
	public final EllipseRoi roi;
	public final double[] slopeandintercept;
	public final double[] mean;
	public final double[] covar;
	
	public LabelledImg(final int label, final RandomAccessibleInterval<FloatType> roiimg,
			final RandomAccessibleInterval<FloatType> Actualroiimg,
			final EllipseRoi roi,
			final double[] slopeandintercept,
			final double[] mean,
			final double[] covar){
		
		this.label = label;
		this.roiimg = roiimg;
		this.Actualroiimg = Actualroiimg;
		this.roi = roi;
		this.slopeandintercept = slopeandintercept;
		this.mean = mean;
		this.covar = covar;
	}
	
}

package labeledObjects;

import ij.gui.EllipseRoi;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.numeric.real.FloatType;

public class LabelledImg {

	
	public final int label;
	public final RandomAccessibleInterval<FloatType> roiimg;
	public final EllipseRoi roi;
	public final double[] slopeandintercept;
	
	public LabelledImg(final int label, final RandomAccessibleInterval<FloatType> roiimg, final EllipseRoi roi,
			final double[] slopeandintercept){
		
		this.label = label;
		this.roiimg = roiimg;
		this.roi = roi;
		this.slopeandintercept = slopeandintercept;
	}
	
}

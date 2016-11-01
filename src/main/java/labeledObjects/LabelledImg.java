package labeledObjects;

import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.numeric.real.FloatType;

public class LabelledImg {

	
	public final int label;
	public final RandomAccessibleInterval<FloatType> roiimg;
	
	public LabelledImg(final int label, final RandomAccessibleInterval<FloatType> roiimg){
		
		this.label = label;
		this.roiimg = roiimg;
	}
	
}

package labeledObjects;

import ij.gui.EllipseRoi;

public class LabelledRoi {

	
	public int roiindex;
	
	public EllipseRoi roi;
	
	
	public LabelledRoi(final int roiindex, final EllipseRoi roi){
		
		this.roiindex = roiindex;
		this.roi = roi;
		
	}
}

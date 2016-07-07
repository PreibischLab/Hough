package labeledObjects;

import net.imglib2.Point;
import net.imglib2.type.numeric.real.FloatType;

    public final class PreFinalobject {

	public final int Label;
	public final Point centroid;
	public final FloatType Intensity;
	public final double slope;
	public final double intercept;
	

	public PreFinalobject(final int Label, final Point centroid, final FloatType Intensity, final double slope,
			final double intercept) {
		this.Label = Label;
		this.centroid = centroid;
		this.Intensity = Intensity;
		this.slope = slope;
		this.intercept = intercept;

	}
}

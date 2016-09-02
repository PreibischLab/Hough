package labeledObjects;

import net.imglib2.RealPoint;

public class Finalobject {
	public final int Label;
	public final RealPoint centroid;
	public final double Intensity;
	public final double sigmaX;
	public final double sigmaY;
	public final double slope;
	public final double intercept;
	

	public Finalobject(final int Label, final RealPoint centroid, final double Intensity, final double sigmaX,
			final double sigmaY, final double slope,
			final double intercept) {
		this.Label = Label;
		this.centroid = centroid;
		this.Intensity = Intensity;
		this.sigmaX = sigmaX;
		this.sigmaY = sigmaY;
		this.slope = slope;
		this.intercept = intercept;

	}
}

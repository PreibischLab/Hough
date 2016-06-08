package varun;

import net.imglib2.type.numeric.real.FloatType;

public  final class Labelparam {
	final int Label;
	final double[] point;
	final FloatType Value;
	final double slope;
	final double intercept;
	

	protected Labelparam(final int Label, final double[] point, final FloatType Value, final double slope,
			final double intercept) {
		this.Label = Label;
		this.point = point;
		this.Value = Value;
		this.slope = slope;
		this.intercept = intercept;

	}
}

package varun;

import net.imglib2.type.numeric.real.FloatType;

public final class Simulatedline {
	final int Label;
	final double[] point;
	final FloatType Value;

	protected Simulatedline(final int Label, final double[] point, final FloatType Value) {
		this.Label = Label;
		this.point = point;
		this.Value = Value;
		

	}
}

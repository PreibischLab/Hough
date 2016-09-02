package labeledObjects;

import net.imglib2.type.numeric.real.FloatType;

public final class Simulatedline {
	public final int Label;
	public final double[] point;
	public final FloatType Value;

	public Simulatedline(final int Label, final double[] point, final FloatType Value) {
		this.Label = Label;
		this.point = point;
		this.Value = Value;
		

	}
}

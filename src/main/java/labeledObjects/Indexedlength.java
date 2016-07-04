package labeledObjects;



public final class Indexedlength {

	public final int Label;
	public final double length;
	public final double[] startpos;
	public final double[] endpos;
	public final double slope;
	public final double intercept;
	

	public Indexedlength(final int Label, final double length, final double[] startpos, final double[] endpos, 
			final double slope, final double intercept) {
		this.Label = Label;
		this.length = length;
		this.startpos = startpos;
		this.endpos = endpos;
		this.slope = slope;
		this.intercept = intercept;

		
	}

}

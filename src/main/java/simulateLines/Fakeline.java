package simulateLines;

public class Fakeline {

	
	public final double length;
	public final double slope;
	public final double intercept;
	public final double[] startpos;
	public final double[] endpos;
	
	public Fakeline (final double length, final double slope, final double intercept, final double[] startpos, final double[] endpos){
		
		this.length = length;
		this.slope = slope;
		this.intercept = intercept;
		this.startpos = startpos;
		this.endpos = endpos;
		
		
	}
	
}

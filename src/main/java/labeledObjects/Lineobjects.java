package labeledObjects;

// Objects containing the label and the correspoing rho and theta information
	public  final class Lineobjects {
		public final int Label;
		public final double Rho;
		public final double Theta;
		public final double [] boxmin;
		public final double [] boxmax;
		
		

		public Lineobjects(
				final int Label,
				final double Rho, 
				final double Theta, 
				final double [] boxmin, 
				final double [] boxmax
				) {
			this.Label = Label;
			this.Rho = Rho;
			this.Theta = Theta;
			this.boxmin = boxmin;
			this.boxmax = boxmax;
			
			
		}
	}

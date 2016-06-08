package varun;

// Objects containing the label and the correspoing rho and theta information
	public  final class Lineobjects {
		final int Label;
		final double Rho;
		final double Theta;
		final long [] boxmin;
		final long [] boxmax;
		
		

		protected Lineobjects(
				final int Label,
				final double Rho, 
				final double Theta, 
				final long [] boxmin, 
				final long [] boxmax
				) {
			this.Label = Label;
			this.Rho = Rho;
			this.Theta = Theta;
			this.boxmin = boxmin;
			this.boxmax = boxmax;
			
			
		}
	}

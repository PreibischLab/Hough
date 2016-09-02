package labeledObjects;

import java.util.ArrayList;

// Objects containing the label and the correspoing rho and theta information
	public  final class Lineobjects {
		public final int Label;
		//public final double Rho;
		//public final double Theta;
        public final ArrayList<double[]> rhotheta;  		
		public final long [] boxmin;
		public final long [] boxmax;
		
		

		public Lineobjects(
				final int Label,
				final ArrayList<double[]> rhotheta,
				//final double Rho, 
				//final double Theta, 
				final long[] minCorner, 
				final long[] maxCorner
				) {
			this.Label = Label;
		    this.rhotheta = rhotheta;
			
			//	this.Rho = Rho;
		//	this.Theta = Theta;
			this.boxmin = minCorner;
			this.boxmax = maxCorner;
			
			
		}
	}

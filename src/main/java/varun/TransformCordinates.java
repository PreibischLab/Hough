package varun;

import net.imglib2.type.numeric.RealType;

public class TransformCordinates {

	public static <T extends RealType<T>> double[] transformfwd(double[] location, double[] size, double[] min,
			double[] max) {

		int n = location.length;

		double[] delta = new double[n];

		final double[] realpos = new double[n];

		for (int d = 0; d < n; ++d){
			
			delta[d] = (max[d] - min[d]) / size[d];

			realpos[d] = location[d] * delta[d] + min[d];
		}
		return realpos;

	}

	public static <T extends RealType<T>> double[] transformback(double[] location, double[] size, double[] min,
			double[] max) {

		int n = location.length;

		double[] delta = new double[n];

		final double[] realpos = new double[n];

		for (int d = 0; d < n; ++d){
			
			delta[d] = (max[d] - min[d]) / size[d];
		    
			realpos[d] = (location[d] - min[d]) / delta[d];
		}
		return realpos;

	}

}
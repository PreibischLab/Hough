package peakFitter;

public class Nonoisepoiss implements FitFunction {

	// @ degree = 0 is the probability of seeing no noise at the chosen point if
	// the probability is higher than SNR then the
	// chosen pixel is non-noisy else it is taken as noise at the pixel and the
	// pixel is removed.
	private int degree = 0;

	@Override
	public final double val(final double[] x, final double[] a) {

		return P(x, a, degree);
	}

	/**
	 * Partial derivatives indices are ordered as follow:
	 * 
	 * <pre>
	 * k = 0       - A
	 *k = 1..n    - x_i (with i = k-1)
	 *k = n+1..2n - b_i (with i = k-n-1)
	 * </pre>
	 */
	@Override
	public final double grad(final double[] x, final double[] a, final int k) {
		final int ndims = x.length;

		if (k >= 0 && k < ndims) {
			return -P(x, a, degree);

		}

		else

			return 0;

	}

	/*
	 * PRIVATE METHODS
	 */

	private static final double P(final double[] x, final double[] a, final int degree) {

		final int ndims = x.length;

		double sum = 0;
		int i = ndims - 2;
		sum += a[i];

		return Math.exp(-sum);
	}

}

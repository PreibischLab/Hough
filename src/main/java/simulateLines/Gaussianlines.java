package simulateLines;

import java.util.Random;

import drawandOverlay.PushCurves;
import net.imglib2.Interval;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.numeric.real.FloatType;
import preProcessing.GetLocalmaxmin.IntensityType;

public class Gaussianlines {

	public static void Drawsimulatedlines(RandomAccessibleInterval<FloatType> outimg, final Interval range,
			final Random rnd, final Random rndsec, final double[] sigma, final double length, final int numlines) {

		final int n = outimg.numDimensions();

		// Declare the number of lines to be plotted
		

		for (int lineindex = 0; lineindex < numlines; ++lineindex) {

			double[] startpos = new double[n];
			double endpos[] = new double[n];

			for (int d = 0; d < range.numDimensions(); ++d) {
				startpos[d] = rnd.nextDouble() * (range.max(d) - range.min(d)) + range.min(d);

				endpos[d] = rndsec.nextDouble()  * (range.max(d) - range.min(d)) + range.min(d);
				
				if (startpos[d] < length )
					startpos[d] = startpos[d] + 2 * length;
				if (endpos[d] < length )
					endpos[d] = endpos[d] + 2 * length;
				
				if (startpos[d] > outimg.dimension(d) - length)
					startpos[d] = startpos[d] - 2 * length;
				if (endpos[d] > outimg.dimension(d) - length)
					endpos[d] = endpos[d] - 2 * length;
				
				
				System.out.println("Start: " + startpos[d]);
				System.out.println("End: " + endpos[d]);
			}

			double slope = (endpos[1] - startpos[1]) / (endpos[0] - startpos[0]);
			double intercept = startpos[1] - slope * startpos[0];

			System.out.println(slope + "  " + intercept);

			PushCurves.Drawshortline(outimg, slope, intercept, startpos, endpos, sigma, length);

		}

	}

}

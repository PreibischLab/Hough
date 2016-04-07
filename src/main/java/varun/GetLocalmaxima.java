package varun;

import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;

public class GetLocalmaxima {

	public static void Thresholding(RandomAccessibleInterval<FloatType> img, RandomAccessibleInterval<FloatType> imgout,
			FloatType ThresholdValue) {

		final Cursor<FloatType> bound = Views.iterable(img).localizingCursor();

		final RandomAccess<FloatType> outbound = imgout.randomAccess();

		while (bound.hasNext()) {

			bound.fwd();

			outbound.setPosition(bound);

			if (bound.get().compareTo(ThresholdValue) > 0) {

				outbound.get().set(bound.get());

			}

			else {

				outbound.get().setZero();

			}

		}
	}

}

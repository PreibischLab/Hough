package varun;

import net.imglib2.Cursor;
import net.imglib2.Interval;
import net.imglib2.Localizable;
import net.imglib2.Point;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealPoint;
import net.imglib2.RealRandomAccess;
import net.imglib2.algorithm.fft2.FFTConvolution;
import net.imglib2.algorithm.gradient.PartialDerivative;
import net.imglib2.algorithm.neighborhood.Neighborhood;
import net.imglib2.algorithm.neighborhood.RectangleShape;
//import net.imglib2.algorithm.fft2.FFTConvolution;
import net.imglib2.algorithm.region.hypersphere.HyperSphere;
import net.imglib2.algorithm.region.hypersphere.HyperSphereCursor;
import net.imglib2.algorithm.stats.Normalize;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.interpolation.randomaccess.NLinearInterpolatorFactory;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.complex.ComplexFloatType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Intervals;
import net.imglib2.view.Views;
import varun.GetLocalmaxmin.IntensityType;

public class Kernels {

	public static void ButterflyKernel(final RandomAccessibleInterval<FloatType> inputimage) {

		final float[] butterflyKernel = new float[] { 0, -2, 0, 1, 2, 1, 0, -2, 0 };

		final Img<FloatType> Butterfly = ArrayImgs.floats(butterflyKernel, new long[] { 3, 3 });

		// apply convolution to convolve input data with kernels

		new FFTConvolution<FloatType>(inputimage, Butterfly, new ArrayImgFactory<ComplexFloatType>()).convolve();
	}

	public static void GeneralButterflyKernel(final RandomAccessibleInterval<FloatType> inputimage, double linewidth,
			double dtheta) {
		final double size = Math.sqrt((inputimage.dimension(0) * inputimage.dimension(0)
				+ inputimage.dimension(1) * inputimage.dimension(1)));
		final double gradientmaggth = linewidth / Math.sin(Math.toRadians(dtheta));
		final double mask = -(2 * gradientmaggth + size) / 2;
		final float[] butterflyKernel = new float[] { 0, (float) mask, 0, (float) gradientmaggth, (float) size,
				(float) gradientmaggth, 0, (float) mask, 0 };

		final Img<FloatType> Butterfly = ArrayImgs.floats(butterflyKernel, new long[] { 3, 3 });

		// apply convolution to convolve input data with kernels

		new FFTConvolution<FloatType>(inputimage, Butterfly, new ArrayImgFactory<ComplexFloatType>()).convolve();
		ImageJFunctions.show(Butterfly);
	}

	public static void BigButterflyKernel(final RandomAccessibleInterval<FloatType> inputimage) {

		final float[] butterflyKernel = new float[] { -10, -15, -22, -22, -22, -22, -22, -15, -10, -1, -6, -13, -22,
				-22, -22, -13, -6, -1, 3, 6, 4, -3, -22, -3, 4, 6, 3, 3, 11, 19, 28, 42, 28, 19, 11, 3, 3, 11, 27, 42,
				42, 42, 27, 11, 3, 3, 11, 19, 28, 42, 28, 19, 11, 3, 3, 6, 4, -3, -22, -3, 4, 6, 3, -1, -6, -13, -22,
				-22, -22, -13, -6, -1, -10, -15, -22, -22, -22, -22, -22, -15, -10

		};

		final Img<FloatType> Butterfly = ArrayImgs.floats(butterflyKernel, new long[] { 9, 9 });

		// apply convolution to convolve input data with kernels

		new FFTConvolution<FloatType>(inputimage, Butterfly, new ArrayImgFactory<ComplexFloatType>()).convolve();

	}

	public static void SobelXFilter(final RandomAccessibleInterval<FloatType> inputimage) {

		// create sobel edge filter kernels
		final float[] sX = new float[] { 1, 0, -1, 2, 0, -2, 1, 0, -1 };
		final Img<FloatType> sobelX = ArrayImgs.floats(sX, new long[] { 3, 3 });

		// apply convolution to convolve input data with kernels

		new FFTConvolution<FloatType>(inputimage, sobelX, new ArrayImgFactory<ComplexFloatType>()).convolve();

	}

	public static void SobelYFilter(final RandomAccessibleInterval<FloatType> inputimage) {

		// create sobel edge filter kernels
		final float[] sY = new float[] { 1, 2, -1, 0, 0, 0, -1, -2, -1 };
		final Img<FloatType> sobelY = ArrayImgs.floats(sY, new long[] { 3, 3 });

		// apply convolution to convolve input data with kernels

		new FFTConvolution<FloatType>(inputimage, sobelY, new ArrayImgFactory<ComplexFloatType>()).convolve();

	}

	public static void Edgedetector(final RandomAccessibleInterval<FloatType> inputimage) {
		final float[] HorizontalEdgeFilterKernel = new float[] { -1, 0, 1, -1, 0, 1, -1, 0, 1 };

		final float[] VerticalEdgeFilterKernel = new float[] { -1, -1, -1, 0, 0, 0, 1, 1, 1 };

		final Img<FloatType> HorizontalEdgeFilter = ArrayImgs.floats(HorizontalEdgeFilterKernel, new long[] { 3, 3 });
		final Img<FloatType> VerticalEdgeFilter = ArrayImgs.floats(VerticalEdgeFilterKernel, new long[] { 3, 3 });
		// apply convolution to convolve input data with kernels

		new FFTConvolution<FloatType>(inputimage, HorizontalEdgeFilter, new ArrayImgFactory<ComplexFloatType>())
				.convolve();
		new FFTConvolution<FloatType>(inputimage, VerticalEdgeFilter, new ArrayImgFactory<ComplexFloatType>())
				.convolve();

	}

	// Do mean filtering on the inputimage
	public static void MeanFilter(RandomAccessibleInterval<FloatType> inputimage,
			RandomAccessibleInterval<FloatType> outimage, double sigma) {
		// Mean filtering for a given sigma
		Cursor<FloatType> cursorInput = Views.iterable(inputimage).cursor();
		Cursor<FloatType> cursorOutput = Views.iterable(outimage).cursor();
		FloatType mean = Views.iterable(inputimage).firstElement().createVariable();
		while (cursorInput.hasNext()) {
			cursorInput.fwd();
			cursorOutput.fwd();
			HyperSphere<FloatType> hyperSphere = new HyperSphere<FloatType>(Views.extendMirrorSingle(inputimage),
					cursorInput, (long) sigma);
			HyperSphereCursor<FloatType> cursorsphere = hyperSphere.cursor();
			cursorsphere.fwd();
			mean.set(cursorsphere.get());
			int n = 1;
			while (cursorsphere.hasNext()) {
				cursorsphere.fwd();
				n++;
				mean.add(cursorsphere.get());
			}
			mean.div(new FloatType(n));
			cursorOutput.get().set(mean);
		}

	}

	// Naive Edge detector, first get the gradient of the image, then get local
	// maxima

	public static RandomAccessibleInterval<FloatType> NaiveEdge(RandomAccessibleInterval<FloatType> inputimg,
			FloatType minval, FloatType maxval, double[] sigma, boolean Thresholding) {
		RandomAccessibleInterval<FloatType> imgout = new ArrayImgFactory<FloatType>().create(inputimg, new FloatType());
		RandomAccessibleInterval<FloatType> premaximgout = new ArrayImgFactory<FloatType>().create(inputimg,
				new FloatType());
		RandomAccessibleInterval<FloatType> maximgout = new ArrayImgFactory<FloatType>().create(inputimg,
				new FloatType());
		// Compute gradient of the image
		imgout = GetLocalmaxmin.GradientmagnitudeImage(inputimg);

		Normalize.normalize(Views.iterable(imgout), minval, maxval);
		// Get the global threshold value for the gradient image
		final Float val = GlobalThresholding.AutomaticThresholding(Views.iterable(imgout));
		// Do conditional maxima on the pixels above the global threshold value
		premaximgout = GetLocalmaxmin.FindConditionalLocalMaxima(imgout, new ArrayImgFactory<FloatType>(),
				IntensityType.Gaussian, sigma, val);
		// Compute global threshold for the premaximgout
		final Float valsec = GlobalThresholding.AutomaticThresholding(Views.iterable(premaximgout));
		// Choose the se the pixels below valsec to 0 or not
		if (Thresholding) {
			GetLocalmaxmin.Thresholding(premaximgout, maximgout, valsec, IntensityType.Gaussian, sigma);
			return maximgout;
		} else {
			return premaximgout;
		}

	}

	

	public static RandomAccessibleInterval<FloatType> SimplifiedEdge(RandomAccessibleInterval<FloatType> inputimg,
			FloatType minval, FloatType maxval, double[] sigma, boolean Thresholding) {
		int n = inputimg.numDimensions();
		RandomAccessibleInterval<FloatType> gradientimg = new ArrayImgFactory<FloatType>().create(inputimg,
				new FloatType());
		gradientimg = GetLocalmaxmin.GradientmagnitudeImage(inputimg);
		Cursor<FloatType> cursor = Views.iterable(inputimg).localizingCursor();
		
		RandomAccess<FloatType> randomAccess = inputimg.randomAccess();
		
		final double[] direction = new double[n];
        
        double gradient=0;
		// iterate over all pixels
		while (cursor.hasNext()) {
			// Initialize a point
			cursor.fwd();
			// compute gradient and its direction in each dimension and move along the direction
			
			for (int d = 0; d < inputimg.numDimensions(); ++d) {
				// move one pixel back in dimension d
				randomAccess.bck(d);

				// get the value
				double Back = randomAccess.get().getRealDouble();

				// move twice forward in dimension d, i.e.
				// one pixel above the location of the cursor
				randomAccess.fwd(d);
				randomAccess.fwd(d);

				// get the value
				double Fwd = randomAccess.get().getRealDouble();

				gradient += ((Fwd - Back) * (Fwd - Back)) / 4;

				direction[d] = (Fwd - Back) / 2;
				
			}
			for (int d = 0; d < inputimg.numDimensions(); ++d) 
				
			direction[d] = direction[d] / gradient;
			
			
			
			
			
		}			
					System.out.println(direction[0] + "  "+ direction[1]);
					
					
			/*	HyperSphere<FloatType> hyperSphere = new HyperSphere<FloatType>( inputimg,randomAccess, span );
				
				final Float centerValue = randomAccess.getFloatPosition(d);
				
				boolean isMaximum = true;

				// check if all pixels in the local neighborhood that are smaller
				for (final FloatType value : hyperSphere) {
					// test if the center is smaller than the current pixel value
					if (centerValue < value.get()) {
						isMaximum = false;
						
					}
				}
				if (isMaximum) {
					final RandomAccess<FloatType> outbound = gradientimg.randomAccess();
					outbound.setPosition(randomAccess.getLongPosition(d),d);

					randomAccess.localize(position);
					
						AddGaussian.addGaussian(gradientimg, position, sigma, false);
					
					}
				randomAccess.move(Math.round(direction[d]),d);
				
				if (randomAccess.getDoublePosition(d) >= inputimg.dimension(d) || randomAccess.getDoublePosition(d) <= 0)
					break;
				
				
				}
				
			
			
			*/

		
		return inputimg;

	}



}
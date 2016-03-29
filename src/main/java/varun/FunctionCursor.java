package varun;


import net.imglib2.Cursor;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RandomAccess;
import net.imglib2.Sampler;
import net.imglib2.view.Views;

public class FunctionCursor < T > implements Cursor<T> {

	final RandomAccessibleInterval<T> source;

	final protected double[] funcparamone;

	final protected RandomAccess<T> randomAccess;
	
	final protected Cursor<T> cursor;
	
	final protected double funcparamtwo;
	
	final protected double funcparamthree;

	final int numDimensions, maxDim;

	// the remaining number of steps in each dimension we still have to go
	final long[] nextsteps;
	

	public FunctionCursor( final RandomAccessibleInterval< T > source, final double[] funcparamone, final double funcparamtwo,
			final double funcparamthree)
	{
		this.source = source;
		this.funcparamone = funcparamone.clone();
		this.funcparamtwo = funcparamtwo;
		this.funcparamthree = funcparamthree;
		this.numDimensions = source.numDimensions();
		this.maxDim = numDimensions-1;
		this.randomAccess = source.randomAccess();
		this.cursor = Views.iterable(source).localizingCursor();
		this.nextsteps = new long[ numDimensions ];

		
	}

	public FunctionCursor( final FunctionCursor< T > cursor )
	{
		this.source = cursor.source;
		this.funcparamone = cursor.funcparamone.clone();
		this.funcparamtwo = cursor.funcparamtwo;
		this.funcparamthree = cursor.funcparamthree;
		this.numDimensions = cursor.numDimensions();
		this.randomAccess = source.randomAccess();
		this.cursor = Views.iterable(source).localizingCursor();
		this.maxDim = cursor.maxDim;
		this.nextsteps = cursor.nextsteps;
		this.randomAccess.setPosition( cursor.randomAccess );
		
	}
	
	
	
	
	
	@Override
	public void localize(float[] position) {

		randomAccess.localize( position );		
		
	}

	@Override
	public void localize(double[] position) {

		randomAccess.localize( position );		
		
	}

	@Override
	public float getFloatPosition(int d) {
		
		return randomAccess.getFloatPosition( d );
	}

	@Override
	public double getDoublePosition(int d) {
		
		return randomAccess.getDoublePosition( d );
	}

	@Override
	public int numDimensions() {

		return numDimensions;
	}

	@Override
	public T get() {
		return randomAccess.get();
	}

	@Override
	public Sampler<T> copy() {
		
		return new FunctionCursor< T >( this );
	}

	@Override
	public void jumpFwd(long steps) {
		for ( long j = 0; j < steps; ++j )
			fwd();
		
	}

	@Override
	public void fwd() {
		
		
			
			
			
	}

	@Override
	public void reset() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public boolean hasNext() {
		
			return nextsteps[ maxDim ] > 0;
	}

	@Override
	public T next() {
		fwd();
		return get();
	}

	@Override
	public void remove() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void localize(int[] position) {
		
		randomAccess.localize( position );		
	}

	@Override
	public void localize(long[] position) {
		
		randomAccess.localize( position );		
	}

	@Override
	public int getIntPosition(int d) {
		
		return randomAccess.getIntPosition( d );
	}

	@Override
	public long getLongPosition(int d) {
	
		return randomAccess.getLongPosition( d );
		
	}

	@Override
	public Cursor<T> copyCursor() {
		return copyCursor();
	}

}

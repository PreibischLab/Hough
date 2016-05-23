package varun;

import java.io.DataOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.io.Writer;

import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.cell.CellImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.Type;
import net.imglib2.Cursor;
import net.imglib2.IterableInterval;
import net.imglib2.Point;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import util.ImgLib2Util;
import ij.ImageJ;
import ij.ImagePlus;
import ij.io.Opener;
import ij.process.ImageProcessor;

import java.io.File;
import java.util.ArrayList;

	


	public class Copy {
		
		
		
		
		
		public static void main(String[] args) throws IOException  {
			Img<FloatType> img1 = ImgLib2Util.openAs32Bit(new File("src/main/resources/box_canny_input.png"));
		System.out.println(img1.dimension(0) + " "+ img1.dimension(1));
			int n = img1.numDimensions();
			
			
			 File fac = new File("src/main/resources/image_intensity.dat");
	           if (!fac.exists())
	           {
	               fac.createNewFile();
	           }
	           System.out.println("The file has been created.");
	           
	           FileWriter write = new FileWriter(fac);

            
			
	        Cursor<FloatType> cursor = Views.iterable(img1).cursor();
	        final double[] direction = new double[n];
			// Extend the input image for gradient computation
			RandomAccessible<FloatType> view = Views.extendMirrorSingle(img1);
			RandomAccess<FloatType> randomAccess = view.randomAccess();
	        write.write("");
	        while(cursor.hasNext()){
	        	cursor.fwd();
	        	double gradient = 0;
	        	for (int d = 0; d < img1.numDimensions(); ++d){
	        		randomAccess.setPosition(cursor);
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

	        		
	        		
	        		
	        		
	        	//	write.append(out+ "\n");
	                 
	        	}
	        	// Normalize the gradient direction
				
				for (int d = 0; d < img1.numDimensions(); ++d) {
					if (gradient != 0)
						direction[d] = direction[d] / gradient;
					else
						direction[d] = Double.MAX_VALUE;
				
					write.append(direction[d]+ "\n");
				}
				
	        	
	        	
	        }
	        write.close();
		}
	}



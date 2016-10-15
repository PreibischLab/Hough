package drawandOverlay;


import java.awt.Color;
import java.util.ArrayList;

import org.jgrapht.graph.DefaultWeightedEdge;

import fiji.tool.SliceListener;
import fiji.tool.SliceObserver;
import graphconstructs.Staticproperties;
import ij.ImagePlus;
import ij.gui.Line;
import ij.gui.Overlay;
import labeledObjects.Subgraphs;

public class DisplaysubGraphend {
	
	private final ImagePlus imp;
	private final ArrayList<Subgraphs> subgraph;
    private final int ndims;
	
	public DisplaysubGraphend(final ImagePlus imp, final ArrayList<Subgraphs> subgraph){
		
		this.imp = imp;
		this.subgraph = subgraph;
		ndims = imp.getNDimensions();
		

		// add listener to the imageplus slice slider
				SliceObserver sliceObserver = new SliceObserver( imp, new ImagePlusListener() );
	}
	
public ImagePlus getImp() { return this.imp; } 
	
	
	protected  class ImagePlusListener implements SliceListener
	{
		@Override
		public void sliceChanged(ImagePlus arg0)
		{
			
			
			int maxSlice = subgraph.get(subgraph.size() - 1).Currentframe + 1;
			
			imp.show();
			
			
			Overlay o = imp.getOverlay();
			
			if( getImp().getOverlay() == null )
			{
				o = new Overlay();
				getImp().setOverlay( o ); 
			}

			o.clear();
			
			int currentSlice = getImp().getCurrentSlice();
			for (int index = 0; index < subgraph.size(); ++index){
			
				if (currentSlice == subgraph.get(index).Previousframe + 1 && currentSlice < maxSlice){
					for (DefaultWeightedEdge e : subgraph.get(index).subgraph.edgeSet()){
						double[] startedge = subgraph.get(index).subgraph.getEdgeSource(e);
						double[] targetedge = subgraph.get(index).subgraph.getEdgeTarget(e);
					        
					        
					        Line newline = new Line(startedge[0], startedge[1], targetedge[0], targetedge[1]);
							newline.setStrokeColor(Color.GREEN);
							newline.setStrokeWidth(0.8);

							o.add(newline);
							
						
					}
				}
				
				imp.updateAndDraw();
			}
			
			
			
			
			System.out.println( arg0.getCurrentSlice() );
		}		
	}
	
}


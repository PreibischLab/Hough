package velocityanalyser;
import java.util.ArrayList;
import java.util.Iterator;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.SimpleWeightedGraph;

import graphconstructs.Logger;
import graphconstructs.Staticproperties;

public class Trackstart implements Linetracker {

	

		private final ArrayList<ArrayList<Staticproperties>> Allstartandend;
		private final long maxframe;
		private SimpleWeightedGraph< double[], DefaultWeightedEdge > graph;
		protected Logger logger = Logger.DEFAULT_LOGGER;
		protected String errorMessage;

		public Trackstart(
				final ArrayList<ArrayList<Staticproperties>> Allstartandend,  
				final long maxframe){
			this.Allstartandend = Allstartandend;
			this.maxframe = maxframe;
			
			
		}
		
		
		
		

		@Override
		public boolean process() {

			reset();
			
			
			for (int frame = 0; frame < maxframe - 1  ; ++frame){
			
			
				ArrayList<Staticproperties> Baseframestartend = Allstartandend.get(frame);
				
				
				
				Iterator<Staticproperties> baseobjectiterator = Baseframestartend.iterator();
				
				
		      
				
				while(baseobjectiterator.hasNext()){
					
					final Staticproperties source = baseobjectiterator.next();
					
					
					double sqdist = Distance(source.oldstartpoint, source.newstartpoint);
					
					synchronized (graph) {
						
						graph.addVertex(source.oldstartpoint);
						graph.addVertex(source.newstartpoint);
						final DefaultWeightedEdge edge = graph.addEdge(source.oldstartpoint, source.newstartpoint);
						graph.setEdgeWeight(edge, sqdist);
						
						
					}
				
			       
				}
				
				System.out.println("Moving to next frame!");
			}
			
			
				return true;
				
			}
		

		@Override
		public void setLogger( final Logger logger) {
			this.logger = logger;
			
		}
		

		@Override
		public SimpleWeightedGraph< double[], DefaultWeightedEdge > getResult()
		{
			return graph;
		}
		
		@Override
		public boolean checkInput() {
			final StringBuilder errrorHolder = new StringBuilder();
			final boolean ok = checkInput();
			if (!ok) {
				errorMessage = errrorHolder.toString();
			}
			return ok;
		}
		
		public void reset() {
			graph = new SimpleWeightedGraph<double[], DefaultWeightedEdge>(DefaultWeightedEdge.class);
//			final Iterator<Staticproperties> it = Allstartandend.get(0).iterator();
				graph.addVertex(Allstartandend.get(0).get(0).oldstartpoint);
		}

		@Override
		public String getErrorMessage() {
			
			return errorMessage;
		}
		
		
		public double Distance(final double[] cordone, final double[] cordtwo) {

			double distance = 0;

			for (int d = 0; d < cordone.length; ++d) {

				distance += Math.pow((cordone[d] - cordtwo[d]), 2);

			}
			return Math.sqrt(distance);
		}
	}



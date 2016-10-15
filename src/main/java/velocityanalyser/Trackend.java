package velocityanalyser;
import java.util.ArrayList;
import java.util.Iterator;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.SimpleWeightedGraph;

import graphconstructs.Logger;
import graphconstructs.Staticproperties;
import labeledObjects.Subgraphs;

public class Trackend implements Linetracker {

	

		private final ArrayList<ArrayList<Staticproperties>> Allstartandend;
		private final long maxframe;
		private SimpleWeightedGraph< double[], DefaultWeightedEdge > graph;
		private ArrayList<Subgraphs> Framedgraph;
		protected Logger logger = Logger.DEFAULT_LOGGER;
		protected String errorMessage;

		public Trackend(
				final ArrayList<ArrayList<Staticproperties>> Allstartandend,  
				final long maxframe){
			this.Allstartandend = Allstartandend;
			this.maxframe = maxframe;
			
			
		}
		
		public ArrayList<Subgraphs> getFramedgraph() {

			return Framedgraph;
		}
		
		

		@Override
		public boolean process() {

			reset();
			
			/*
			 * Outputs
			 */

			graph = new SimpleWeightedGraph<double[], DefaultWeightedEdge>(DefaultWeightedEdge.class);
			Framedgraph = new ArrayList<Subgraphs>();
			for (int frame = 1; frame < maxframe   ; ++frame){
			
			
				ArrayList<Staticproperties> Baseframestartend = Allstartandend.get(frame - 1);
				
				
				
				Iterator<Staticproperties> baseobjectiterator = Baseframestartend.iterator();
				
				SimpleWeightedGraph<double[], DefaultWeightedEdge> subgraph = new SimpleWeightedGraph<double[], DefaultWeightedEdge>(
						DefaultWeightedEdge.class);
		      
				
				while(baseobjectiterator.hasNext()){
					
					final Staticproperties source = baseobjectiterator.next();
					
					
					double sqdist = Distance(source.oldendpoint, source.newendpoint);
					
					synchronized (graph) {
						
						graph.addVertex(source.oldendpoint);
						graph.addVertex(source.newendpoint);
						final DefaultWeightedEdge edge = graph.addEdge(source.oldendpoint, source.newendpoint);
						graph.setEdgeWeight(edge, sqdist);
						
						
					}
					subgraph.addVertex(source.oldendpoint);
					subgraph.addVertex(source.newendpoint);
					final DefaultWeightedEdge subedge = subgraph.addEdge(source.oldendpoint, source.newendpoint);
					subgraph.setEdgeWeight(subedge, sqdist);

					Subgraphs currentframegraph = new Subgraphs(frame - 1, frame, subgraph);
					Framedgraph.add(currentframegraph);
			       
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
			final StringBuilder errrorHolder = new StringBuilder();;
			final boolean ok = checkInput();
			if (!ok) {
				errorMessage = errrorHolder.toString();
			}
			return ok;
		}
		
		public void reset() {
			graph = new SimpleWeightedGraph<double[], DefaultWeightedEdge>(DefaultWeightedEdge.class);

			graph.addVertex(Allstartandend.get(0).get(0).oldendpoint);
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



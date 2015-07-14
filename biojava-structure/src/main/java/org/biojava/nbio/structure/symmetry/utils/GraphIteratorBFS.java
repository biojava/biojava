package org.biojava.nbio.structure.symmetry.utils;

import java.util.Iterator;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.Set;
import java.util.TreeSet;

/**
 * Graph traverse by breath first search (BFS).
 * <p>
 * Graph iterators can traverse any {@link Graph} 
 * implementation starting from any vertex in the graph 
 * as the source.
 * 
 * @author Aleix Lafita
 * @since 4.1.1
 *
 */
public class GraphIteratorBFS<V> implements Iterator<V> {

	private Graph<V> graph;
	private Queue<Integer> queue; //stores the vertices to visit next
	private Set<Integer> visited; //stores the vertices already visited

	/**
	 * Constructor. Builds an iterator of the graph starting
	 * from the vertex source specified.
	 * 
	 * @param graph
	 * @param source vertex index to start the search
	 */
	public GraphIteratorBFS(Graph<V> graph, int source){

		this.graph = graph;
		queue = new PriorityQueue<Integer>();
		visited = new TreeSet<Integer>();
		queue.add(source);
		visited.add(source);
	}

	@Override
	public boolean hasNext() {
		return !queue.isEmpty();
	}

	@Override
	public V next() {

		Integer vertex = queue.poll();
		List<Integer> children = graph.getChildren(vertex);

		for (Integer child : children){
			
			if (!visited.contains(child)){
				queue.add(child);
				visited.add(child);
			}
		}
		return graph.getVertex(vertex);
	}

	@Override
	public void remove() {}
	
}

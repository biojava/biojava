package org.biojava.nbio.structure.symmetry.utils;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.Stack;
import java.util.TreeSet;

/**
 * Graph traverse by depth first search (DFS).
 * <p>
 * Graph iterators that can traverse any {@link Graph} 
 * implementation starting from any vertex in the graph 
 * as the source.
 * 
 * @author Aleix Lafita
 * @since 4.1.1
 *
 */
public class GraphIteratorDFS<V> implements Iterator<V> {

	private Graph<V> graph;
	private Stack<Integer> path; //stores the current vertex path
	private Stack<List<Integer>> stack; //stores the vertices to visit next
	private Set<Integer> visited; //stores the vertices already visited

	/**
	 * Constructor. Builds an iterator of the graph starting
	 * from the vertex source specified.
	 * 
	 * @param graph
	 * @param source vertex index to start the search
	 */
	public GraphIteratorDFS(Graph<V> graph, int source){

		this.graph = graph;
		path = new Stack<Integer>();
		stack = new Stack<List<Integer>>();
		visited = new TreeSet<Integer>();

		//First index is the vertex index, second is the depth level
		List<Integer> vertex = new ArrayList<Integer>();
		vertex.add(source);
		vertex.add(0);
		stack.push(vertex);
	}

	@Override
	public boolean hasNext() {
		return !stack.isEmpty();
	}

	@Override
	public V next() {

		List<Integer> vertex = stack.pop();
		Integer index = vertex.get(0);
		int depth = vertex.get(1);

		while (depth < path.size()){
			path.pop();
		}

		if (!path.contains(index)){
			path.push(index);
			List<Integer> children = graph.getChildren(index);

			for (Integer c : children){
				//process all children of the vertex
				List<Integer> child = new ArrayList<Integer>();
				child.add(c);
				child.add(path.size());
				//conditions: not in path, not visited
				if (!path.contains(c) && !visited.contains(c)){
					stack.push(child);
					visited.add(c);
				}
			}
		}
		return graph.getVertex(index);
	}

	@Override
	public void remove() {}
	
	/**
	 * Returns the current exploring path of the DFS.
	 * The order of the elements returned is sorted by
	 * their depth (so that the source is the first 
	 * element of the List).
	 * <p>
	 * Useful for checking for back edges (cycles) in a
	 * graph.
	 *  
	 * @return List of vertices in the path.
	 */
	public List<V> getPath(){
		
		List<V> elements = new ArrayList<V>();
		for (Integer index : path){
			elements.add(graph.getVertex(index));
		}
		return elements;
	}
}

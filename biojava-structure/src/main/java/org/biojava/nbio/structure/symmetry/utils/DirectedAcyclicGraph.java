package org.biojava.nbio.structure.symmetry.utils;

import java.util.List;

/**
 * A directed acyclic graph (DAG) implementation based on a 
 * List of vertices and a List of Edges for each vertex 
 * (adjacency List).
 * <p>
 * Every time a new Edge is added to the Graph, a check is
 * made to ensure that no cycles are present. If a cycle
 * exists, an Exception is thrown.
 * 
 * @author Aleix Lafita
 * @since 4.1.1
 * 
 */
public class DirectedAcyclicGraph<V> extends DirectedGraph<V> {

    @Override
    public boolean addEdge(V vertex1, V vertex2) {
    	
    	//Check if the vertices exist
        int index1 = vertices.indexOf(vertex1);
        int index2 = vertices.indexOf(vertex2);
        if (index1 == -1 || index2 == -1) return false;
        
        //Check if the edge already existed
        boolean changed = vertexMap.get(index1).add(index2);
        if (!changed) return false;
        
        //Check for cycles in the DAG: use a DFS to check for back edges
        for (int source=0; source<size(); source++){
        	GraphIteratorDFS<V> iter = new GraphIteratorDFS<V>(this, source);
        	while (iter.hasNext()){
        		V vertex = iter.next();
        		List<V> path = iter.getPath();
        		for (Integer c : getChildren(indexOf(vertex))){
        			V child = getVertex(c);
        			if (path.contains(child)){
        				//We found a back edge
        				throw new IllegalStateException(
        						"There is a cycle in the DAG");
        			}
        		}
        	}
        }
        return true;
    }
    

    @Override
    public DirectedAcyclicGraph<V> clone() {
    	return (DirectedAcyclicGraph<V>) super.clone();
    }
    
    @Override
    public DirectedAcyclicGraph<V> extractSubGraph(List<Integer> indices) {
    	return (DirectedAcyclicGraph<V>) super.extractSubGraph(indices);
    }
    
}

package org.biojava.nbio.structure.symmetry.utils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

/**
 * A directed graph implementation based on a 
 * List of vertices and a List of Edges for each 
 * vertex (adjacency List).
 * 
 * @author Aleix Lafita
 * 
 */
public class DirectedGraph<V> implements Graph<V> {

    private List<V> vertices = new ArrayList<V>();
    private List<Set<Integer>> vertexMap = new ArrayList<Set<Integer>>();
    
    @Override
    public boolean addEdge(V vertex1, V vertex2) {
        int index1 = vertices.indexOf(vertex1);
        int index2 = vertices.indexOf(vertex2);
        if (index1 == -1 || index2 == -1) return false;
        vertexMap.get(index1).add(index2);
        return true;
    }
    
    @Override
    public boolean addEdge(Edge<V> edge) {
        return addEdge(edge.getVertex1(), edge.getVertex2());
    }
    
    @Override
    public int size() {
        return vertices.size();
    }
    
    @Override
    public int getEdgeCount() {
        int edgeCount = 0;
        for (int i = 0; i < size(); i++){
            edgeCount += getValence(i);
        }
        return edgeCount;
    }
    
    @Override
    public int getValence(V vertex) {
        return getValence(indexOf(vertex));
    }
    
    @Override
    public int getValence(int index) {
        if (index != -1) {
            return getNeighborIndices(index).size();
        }
        return 0;
    }
    
    @Override
    public boolean addVertex(V vertex) {
        if (vertices == null) vertices = new ArrayList<V>();
        if (containsVertex(vertex)) return false;
        vertices.add(vertex);
        if (vertexMap == null) vertexMap = new ArrayList<Set<Integer>>();
        vertexMap.add(new TreeSet<Integer>());
        return true;
    }
    
    @Override
    public boolean containsVertex(V vertex) {
        return vertices.contains(vertex);
    }
    
    @Override
    public void setVertices(List<V> list) {
        vertices = list;
        vertexMap = new ArrayList<Set<Integer>>(vertices.size());
        for (int i = 0; i < vertices.size(); i++) {
            vertexMap.add(new TreeSet<Integer>());
        }
    }
    
    @Override
    public List<V> getVertices() {
        return vertices;
    }
    
    @Override
    public List<Edge<V>> getEdges() {
        List<Edge<V>> edges = new ArrayList<Edge<V>>();
        for(int index1 = 0; index1 < size(); index1++) {
            for(int index2 : getChildren(index1)) {
                Edge<V> e = new Edge<V>(vertices.get(index1),
                		vertices.get(index2));
                edges.add(e);
            }
        }
        return edges;
    }
    
    @Override
    public V getVertex(int index) {
        return vertices.get(index);
    }
    
    @Override
	public String toString() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < vertices.size(); i++) {
            sb.append(vertices.get(i));
            sb.append(" [");
            List<Integer> children = getChildren(i);
            if (children != null) {
                for (int j = 0; j < children.size(); j++) {
                    sb.append(vertices.get(children.get(j)));
                    if (j < children.size()-1) {
                        sb.append(",");
                    }
                }
            }
            sb.append("]; ");
        }
        return sb.toString();
    }
    
    @Override
    public int indexOf(V vertex) {
        return vertices.indexOf(vertex);
    }
    
    @Override
    public List<Integer> getNeighborIndices(int index) {
    	
    	List<Integer> neighbors = new ArrayList<Integer>();
    	neighbors.addAll(getChildren(index));
    	neighbors.addAll(getParents(index));
    	
    	Collections.sort(neighbors);
        return neighbors;
    }
    
    @Override
    public boolean removeEdge(int index1, int index2) {
    	
        if (index1 < 0 || index2 < 0) return false;
        
        List<Integer> children = getChildren(index1);
        int index = children.indexOf(index2);
        
        if (index < 0) return false;
        else children.remove(children.indexOf(index2));
        return true;
    }
    
    public boolean removeEdge(V vertex1, V vertex2) {
        return removeEdge(indexOf(vertex1), indexOf(vertex2));
    }
    
    @Override
    public DirectedGraph<V> clone() {
        DirectedGraph<V> graph = new DirectedGraph<V>();
        //Add all vertices
        for (int i = 0; i < vertices.size(); i++) {
        	V vertex = vertices.get(i);
            graph.addVertex(vertex);
        }
        //Add all Edges
        for (int i = 0; i < vertices.size(); i++) {
        	V vertex = vertices.get(i);
            List<Integer> children = getChildren(i);
            for (int n : children) {
                graph.addEdge(vertex, vertices.get(n));
            }
        }
        return graph;
    }

    @Override
    public boolean containsEdge(int index1, int index2) {
    	
        if (index1 < 0 || index2 < 0) return false;

        List<Integer> children = getChildren(index1);
        int index = children.indexOf(index2);
        
        if (index < 0) return false;
        else return true;
    }
    
    @Override
    public DirectedGraph<V> extractSubGraph(List<Integer> indices) {
    	
    	DirectedGraph<V> graph = new DirectedGraph<V>();
    	//Add all subset vertices
    	for (int index : indices) {
    		V vertex = vertices.get(index);
    		graph.addVertex(vertex);
    	}
    	//Add edges between subset vertices
    	for (int index : indices) {
    		V vertex = vertices.get(index);
    		List<Integer> children = getChildren(index);
    		for (int n: children) {
    			if (indices.contains(n)) {
    				graph.addEdge(vertex, vertices.get(n));
    			}
    		}
    	}
    	return graph;
    }

	@Override
	public List<Integer> getParents(int index) {
		
		List<Integer> parents = new ArrayList<Integer>();
		//Loop through all vertices and find those connected to index
		for (int v=0; v<size(); v++){
			if (v==index) continue;
			for (int child : getChildren(v)){
				if (child == index) parents.add(v);
			}
		}
		return null;
	}

	@Override
	public List<Integer> getChildren(int index) {
		
		List<Integer> children = new ArrayList<Integer>();
    	children.addAll(vertexMap.get(index));
        return children;
	}
	
}

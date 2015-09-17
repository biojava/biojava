/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */
package org.biojava.nbio.structure.symmetry.utils;

import java.util.ArrayList;
import java.util.List;


/**
 * A simple nondirected graph implementation based on a list of edges for each vertex
 * @author Peter
 */
public class SimpleGraph<V> implements Graph<V>, Cloneable {
    private List<V> vertices = new ArrayList<V>();
    private List<ArrayList<Integer>> vertexMap = new ArrayList<ArrayList<Integer>>();
    
    @Override
	public boolean addEdge(V vertex1, V vertex2) {
        int index1 = vertices.indexOf(vertex1);
        int index2 = vertices.indexOf(vertex2);
        if (index1 == -1 || index2 == -1) return false;
        vertexMap.get(index1).add(index2);
        vertexMap.get(index2).add(index1);
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
        for (int i = 0; i < size(); i++)
            edgeCount += getValence(i);
        
        return edgeCount/2;
    }
    
    @Override
	public int getValence(V vertex) {
        return getValence(indexOf(vertex));
    }
    
    @Override
	public int getValence(int index) {
        if (index != -1) {
            return vertexMap.get(index).size();
        }
        return 0;
    }
    
    @Override
	public boolean addVertex(V vertex) {
        if (vertices == null) vertices = new ArrayList<V>();
        if (containsVertex(vertex)) return false;
        vertices.add(vertex);
        if (vertexMap == null) vertexMap = new ArrayList<ArrayList<Integer>>();
        vertexMap.add(new ArrayList<Integer>());
        return true;
    }
    
    @Override
	public boolean containsVertex(V vertex) {
        return vertices.contains(vertex);
    }
    
    @Override
	public void setVertices(List<V> list) {
        vertices = list;
        vertexMap = new ArrayList<ArrayList<Integer>>(vertices.size());
        for (int i = 0; i < vertices.size(); i++) {
            vertexMap.add(new ArrayList<Integer>());
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
            for(int index2 : getNeighborIndices(index1)) {
                if(index1 < index2) {
                    Edge<V> e = new Edge<V>(vertices.get(index1),vertices.get(index2));
                    edges.add(e);
                }
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
            List<Integer> neighbors = vertexMap.get(i);
            if (neighbors != null) {
                for (int j = 0; j < neighbors.size(); j++) {
                    sb.append(vertices.get(neighbors.get(j)));
                    if (j < neighbors.size()-1) {
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
        return vertexMap.get(index);
    }
    
    @Override
	public boolean removeEdge(int index1, int index2) {
        if (index1 < 0 || index2 < 0) return false;
        
        List<Integer> neighbors = getNeighborIndices(index1);
        
        int index = neighbors.indexOf(index2);
        if (index < 0) return false;
        neighbors.remove(neighbors.indexOf(index2));
        
        neighbors = getNeighborIndices(index2);
        index = neighbors.indexOf(index1);
        if (index < 0) return false;
        neighbors.remove(neighbors.indexOf(index1));
        return true;
    }
    
    @Override
	public boolean removeEdge(V vertex1, V vertex2) {
        return removeEdge(indexOf(vertex1), indexOf(vertex2));
    }
    
    @Override
    public Object clone() {
        // clone should not use constructor ??
        SimpleGraph<V> graph = new SimpleGraph<V>();
        
        for (int i = 0; i < vertices.size(); i++) {
            V vertex = vertices.get(i);
            graph.addVertex(vertex);
            List<Integer> neighbors = vertexMap.get(i);
            if (neighbors.size() > 0) {
                for (int n: neighbors) {
                    graph.addEdge(vertex, vertices.get(n));
                }
            }
        }
        return graph;
    }

    @Override
	public boolean containsEdge(int index1, int index2) {
        if (index1 < 0 || index2 < 0) return false;

        List<Integer> neighbors = getNeighborIndices(index1);

        int index = neighbors.indexOf(index2);
        if (index < 0) {
            return false;
        } else {
            return true;
        }
    }
    
    @Override
	public SimpleGraph<V> extractSubGraph(List<Integer> indices) {
    	SimpleGraph<V> graph = new SimpleGraph<V>();
    	for (int index: indices) {
    		V vertex = vertices.get(index);
    		graph.addVertex(vertex);
    		List<Integer> neighbors = getNeighborIndices(index);
    		for (int n: neighbors) {
    			if (indices.contains(n)) {
    				graph.addEdge(vertex, vertices.get(n));
    			}
    		}
    	}
    	return graph;
    }
}

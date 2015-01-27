package org.biojava.bio.structure.symmetry.utils;

import java.util.List;


/**
 *
 * @author Peter
 */
public interface Graph<V> {
    
    public boolean addVertex(V vertex);
    
    public boolean addEdge(V vertex1, V vertex2);
    
    public boolean addEdge(Edge<V> edge);
    
    public void setVertices(List<V> list);
    
    public boolean containsVertex(V vertex);
    
    public boolean containsEdge(int index1, int index2);
    
    public int size();
    
    public int getEdgeCount();
    
    public List<Edge<V>> getEdges();
    
    public int getValence(V vertex);
    
    public int getValence(int index);
    
    public List<V> getVertices();
    
    public V getVertex(int index);
    
    public int indexOf(V vertex);
    
    public List<Integer> getNeighborIndices(int index);
    
    public boolean removeEdge(int index1, int index2);
    
    public boolean removeEdge(V vertex1, V vertex2);
    
    public Graph<V> extractSubGraph(List<Integer> indices);
    
    public Object clone();
}

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

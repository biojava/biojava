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
 * Interface for the representation of a Graph. 
 * Generalization for any implementation, either
 * directed or undirected.
 * 
 * @author Peter
 * @author Aleix Lafita
 * 
 */
public interface Graph<V> {
    
	/**
	 * Add an unconnected vertex to the Graph.
	 * 
	 * @param vertex
	 * @return false if the vertex is already
	 * 			in the Graph, true otherwise.
	 */
    public boolean addVertex(V vertex);
    
    /**
     * Add an Edge between the two vertices to the Graph.
     * Both vertices need to exist already in the Graph, 
     * no new vertices will be added.
     * For undirected Graphs, the reverse Edge will also 
     * be added.
     * 
     * @param vertex1
     * @param vertex2
     * @return true if the Edge is successfully added, 
     * 			false if the Edge already existed or
     * 			any of the vertices are not in the Graph.
     */
    public boolean addEdge(V vertex1, V vertex2);
    
    /**
     * Add an Edge between to the Graph. 
     * Both vertices need to exist already in the Graph, 
     * no new vertices will be added.
     * For undirected Graphs, the reverse Edge will also 
     * be added.
     * 
     * @param edge
     * @return true if the Edge is successfully added, 
     * 			false if the Edge already existed or
     * 			any of the vertices are not in the Graph.
     */
    public boolean addEdge(Edge<V> edge);
    
    /**
     * Resets all vertices of the Graph to the ones in the 
     * input List. It also clears all the Edges of the Graph, 
     * so that the Graph is fully unconnected.
     * 
     * @param List of vertices
     */
    public void setVertices(List<V> list);
    
    /**
     * Checks if a vertex is contained in the Graph.
     * 
     * @param vertex
     * @return true if contained, false otherwise.
     */
    public boolean containsVertex(V vertex);
    
    /**
     * Checks if it exist an Edge between the two
     * vetices in the Graph.
     * 
     * @param index1
     * @param index2
     * @return true if the Edge exists, false otherwise
     * 			or when the indices are out of bounds.
     */
    public boolean containsEdge(int index1, int index2);
    
    /**
     * The size of a Graph is its number of verices.
     * 
     * @return number of vertices in the Graph.
     */
    public int size();
    
    /**
     * Return the number of Edges in the Graph. For
     * undirected Graphs, it returns the number of 
     * unique Edges (they are not counted twice).
     * 
     * @return number of Edges in the Graph.
     */
    public int getEdgeCount();
    
    /**
     * Return a List of all Edges in the Graph. For
     * undirected Graphs, it returns the unique Edges
     * only.
     * 
     * @return List of Edges in the Graph.
     */
    public List<Edge<V>> getEdges();
    
    /**
     * Return the valence (degree) of a vertex, that is,
     * the number of vertices connected to this vertex.
     * 
     * @param vertex
     * @return number of connected vertices to the vertex,
     * 			or 0 if the vertex is not in the Graph.
     */
    public int getValence(V vertex);
    
    /**
     * Return the valence (degree) of a vertex, that is,
     * the number of vertices connected to this vertex.
     * 
     * @param index the vertex index in the Graph.
     * @return number of connected vertices to the vertex,
     * 			or 0 if the vertex is not in the Graph.
     */
    public int getValence(int index);
    
    /**
     * Return a List of all Vertices in the Graph.
     * Their can be sorted in any undetermined way,
     * depending on the implementation used.
     * 
     * @return List of Vertices
     */
    public List<V> getVertices();
    
    /**
     * Return the vertex at the specified index of the
     * Graph.
     * 
     * @param index of the vertex
     * @return vertex at the index
     */
    public V getVertex(int index);
    
    /**
     * Return the index of a vertex in the Graph.
     * Returns -1 if the vertex is not in the Graph.
     * 
     * @param vertex
     * @return index of the vertex, or -1 if not in
     * 			the Graph.
     */
    public int indexOf(V vertex);
    
    /**
     * Return a List of all vertex indices connected to 
     * the vertex.<p>
     * For directed Graphs all child and parent nodes will
     * be returned.
     * 
     * @param index
     * @return List of indicies of connected vertices
     */
    public List<Integer> getNeighborIndices(int index);
    
    /**
     * Return a List of all parent vertex indices connected 
     * to the vertex. The parent vertices are the ones that
     * contain an edge pointing to this vertex.<p>
     * For undirected Graphs there is no distinction 
     * between child and parent nodes, so all neighboring
     * edges will be returned.
     * 
     * @param index
     * @return
     */
    public List<Integer> getParents(int index);
    
    /**
     * Return a List of all child vertex indices connected 
     * to the vertex. The children vertices are the ones that
     * contain an edge pointing out of this vertex.<p>
     * For undirected Graphs there is no distinction 
     * between child and parent nodes, so all neighboring
     * edges will be returned.
     * 
     * @param index
     * @return
     */
    public List<Integer> getChildren(int index);
    
    /**
     * Remove the Edge from the Graph. If the graph is
     * undirected, it also removes the reverse Edge.
     * 
     * @param index1
     * @param index2
     * @return true if the Edge is successfully removed,
     * 			false if the Edge or any or its vertices
     * 			are not in the Graph.
     */
    public boolean removeEdge(int index1, int index2);
    
    /**
     * Remove the Edge from the Graph. If the graph is
     * undirected, it also removes the reverse Edge.
     * 
     * @param vertex1
     * @param vertex2
     * @return true if the Edge is successfully removed,
     * 			false if the Edge or any or its vertices
     * 			are not in the Graph.
     */
    public boolean removeEdge(V vertex1, V vertex2);
    
    /**
     * Return a Graph object that only contains the 
     * specified List of vertices. Only the Edge between
     * these vertices are included in the new Graph.
     * 
     * @param indices of vertices to be included.
     * @return a new Graph with the subset of vertices.
     */
    public Graph<V> extractSubGraph(List<Integer> indices);
    
    /**
     * Return an identical copy of this Graph.
     * It does not clone all the vertex Objects
     * of the Graph, only references are copied.
     * 
     * @return identical copy of this Graph.
     */
    public Graph<V> clone();
}

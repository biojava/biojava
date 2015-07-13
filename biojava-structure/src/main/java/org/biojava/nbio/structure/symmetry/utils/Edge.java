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

/**
 * Stores an instance of an edge for graphs.
 * Represents an edge between two vertex objects in the graph.
 * 
 * @author Peter
 * 
 */
public class Edge<V> {
    V vertex1;
    V vertex2;
    
    /** 
     * Creates a new instance of Edge 
     * 
     * @param vertex1
     * @param vertex2
     */
    public Edge(V vertex1, V vertex2) {
        this.vertex1 = vertex1;
        this.vertex2 = vertex2;
    }
    
    /**
     * Return the first vertex of the Edge.
     */
    public V getVertex1() {
        return vertex1;
    }
    
    /**
     * Return the second vertex of the Edge.
     */
    public V getVertex2() {
        return vertex2;
    }
}

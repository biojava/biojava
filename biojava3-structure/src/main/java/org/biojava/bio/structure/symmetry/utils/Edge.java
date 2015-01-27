package org.biojava.bio.structure.symmetry.utils;

/**
 *
 * @author Peter
 */
public class Edge<V> {
    V vertex1;
    V vertex2;
    
    /** Creates a new instance of Edge */
    public Edge(V vertex1, V vertex2) {
        this.vertex1 = vertex1;
        this.vertex2 = vertex2;
    }
    
    public V getVertex1() {
        return vertex1;
    }
    
    public V getVertex2() {
        return vertex2;
    }
}

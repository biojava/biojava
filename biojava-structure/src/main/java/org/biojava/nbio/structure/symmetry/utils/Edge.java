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

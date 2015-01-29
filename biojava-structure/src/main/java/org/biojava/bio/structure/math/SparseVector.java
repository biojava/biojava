/**
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
 * Created on 5 Mar 2013
 * Created by Andreas Prlic
 *
 * @since 3.0.6
 */
package org.biojava.bio.structure.math;

import java.io.Serializable;


/**
 *  
 *  A sparse vector, implemented using a symbol table.
 *  
 *  Derived from http://introcs.cs.princeton.edu/java/44st/SparseVector.java.html
 *  
 *  For additional documentation, see <a href="http://introcs.cs.princeton.edu/44st">Section 4.4</a> of
 *  <i>Introduction to Programming in Java: An Interdisciplinary Approach</i> by Robert Sedgewick and Kevin Wayne. 
 */

public class SparseVector implements Serializable{
    /**
	 * 
	 */
	private static final long serialVersionUID = 1174668523213431927L;
	
	private final int N;             // length
	
    private SymbolTable<Integer, Double> symbolTable;  // the vector, represented by index-value pairs


    /** Constructor. initialize the all 0s vector of length N
     *  
     * @param N
     */
    public SparseVector(int N) {
        this.N  = N;
        this.symbolTable = new SymbolTable<Integer, Double>();
    }

   /** Setter method (should it be renamed to set?)
    * 
    * @param i set symbolTable[i]
    * @param value
    */
    public void put(int i, double value) {
        if (i < 0 || i >= N) throw new IllegalArgumentException("Illegal index " + i + " should be > 0 and < " + N);
        if (value == 0.0) symbolTable.delete(i);
        else              symbolTable.put(i, value);
    }

    /** get a value
     * 
     * @param i
     * @return  return symbolTable[i]
     */
    public double get(int i) {
        if (i < 0 || i >= N) throw new IllegalArgumentException("Illegal index " + i + " should be > 0 and < " + N);
        if (symbolTable.contains(i)) return symbolTable.get(i);
        else                return 0.0;
    }

    // return the number of nonzero entries
    public int nnz() {
        return symbolTable.size();
    }

    // return the size of the vector
    public int size() {
        return N;
    }

    /** Calculates the dot product of this vector a with b
     * 
     * @param b
     * @return
     */
    public double dot(SparseVector b) {
        SparseVector a = this;
        if (a.N != b.N) throw new IllegalArgumentException("Vector lengths disagree. " + a.N + " != " + b.N);
        double sum = 0.0;

        // iterate over the vector with the fewest nonzeros
        if (a.symbolTable.size() <= b.symbolTable.size()) {
            for (int i : a.symbolTable)
                if (b.symbolTable.contains(i)) sum += a.get(i) * b.get(i);
        }
        else  {
            for (int i : b.symbolTable)
                if (a.symbolTable.contains(i)) sum += a.get(i) * b.get(i);
        }
        return sum;
    }

    /** Calculates the 2-norm
     * 
     * @return
     */
    public double norm() {
        SparseVector a = this;
        return Math.sqrt(a.dot(a));
    }

    /** Calculates  alpha * a
     * 
     * @param alpha
     * @return
     */
    public SparseVector scale(double alpha) {
        SparseVector a = this;
        SparseVector c = new SparseVector(N);
        for (int i : a.symbolTable) c.put(i, alpha * a.get(i));
        return c;
    }

    /** Calcualtes return a + b
     * 
     * @param b
     * @return
     */
    public SparseVector plus(SparseVector b) {
        SparseVector a = this;
        if (a.N != b.N) throw new IllegalArgumentException("Vector lengths disagree : " + a.N + " != " + b.N);
        SparseVector c = new SparseVector(N);
        for (int i : a.symbolTable) c.put(i, a.get(i));                // c = a
        for (int i : b.symbolTable) c.put(i, b.get(i) + c.get(i));     // c = c + b
        return c;
    }

    @Override
    public String toString() {
        StringBuilder s = new StringBuilder();
        for (int i : symbolTable) {
            s.append("(");
            s.append(i);
            s.append(", ");
            s.append(symbolTable.get(i));
            s.append(") ");
        }
        return s.toString();
    }


}


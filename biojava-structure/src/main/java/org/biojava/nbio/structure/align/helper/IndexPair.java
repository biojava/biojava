/*
 *                  BioJava development code
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
 * Created on Mar 1, 2006
 *
 */
package org.biojava.nbio.structure.align.helper;

import java.io.Serializable;
import java.util.Objects;

public class IndexPair implements Serializable {

    /**
     *
     */
    private static final long serialVersionUID = 1832393751152650420L;
    private short row;
    private short col;

	public IndexPair(int row, int col) {
    	this.row = (short)row;
    	this.col = (short)col;
    	if (this.row!=row || this.col!=col)
    		throw new ArithmeticException();
	}

    public IndexPair(short row, short col) {
        this.row = row;
        this.col = col;
    }

    public short col() {
        return col;
    }

    public void setCol(short col) {
        this.col = col;
    }

    public short row() {
        return row;
    }

    public void setRow(short row) {
        this.row = row;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof IndexPair)) return false;
        IndexPair indexPair = (IndexPair) o;
        return row == indexPair.row && col == indexPair.col;
    }

    @Override
    public int hashCode() {
        return Objects.hash(row, col);
    }

    @Override
    public String toString() {
        return "[" + row + " " + col + "]";
    }
}


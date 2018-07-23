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
package org.biojava.nbio.structure.contact;

import java.io.Serializable;

/**
 * A simple class to store contacts in the form of pairs of indices and a distance associated to them.
 * @author Jose Duarte
 * @since 5.0
 */
public class Contact implements Serializable {

	private static final long serialVersionUID = 5059569852079048728L;

	private int i;
	private int j;
	private double distance;
	
	public Contact(int i, int j, double distance) {
		this.i = i;
		this.j = j;
		this.distance = distance;
	}
	
	public Pair<Integer> getIndexPair() {
		return new Pair<Integer>(i,j);
	}
	
	public int getI() {
		return i;
	}
	
	public int getJ() {
		return j;
	}
	
	public double getDistance() {
		return distance;
	}
}

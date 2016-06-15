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

import org.biojava.nbio.structure.Atom;

/**
 * A pair of atoms that are in contact
 *
 * @author duarte_j
 *
 */
public class AtomContact implements Serializable {

	private static final long serialVersionUID = 1L;

	private Pair<Atom> pair;

	private double distance;

	public AtomContact(Pair<Atom> pair, double distance) {
		this.pair = pair;
		this.distance = distance;
	}

	public Pair<Atom> getPair() {
		return pair;
	}

	public void setPair(Pair<Atom> pair) {
		this.pair = pair;
	}

	public double getDistance() {
		return distance;
	}

	public void setDistance(double distance) {
		this.distance = distance;
	}
}

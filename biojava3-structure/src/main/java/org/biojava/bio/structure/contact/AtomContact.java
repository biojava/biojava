package org.biojava.bio.structure.contact;

import org.biojava.bio.structure.Atom;

/**
 * A pair of atoms that are in contact
 * 
 * @author duarte_j
 *
 */
public class AtomContact {
	
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

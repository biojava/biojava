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
package org.biojava.nbio.structure.domain.pdp;

import java.io.Serializable;



public class Segment implements Serializable, Comparable<Segment> {
	/**
	 *
	 */
	private static final long serialVersionUID = 1393487067559539657L;
	private Integer from;
	private Integer to;
	double score;

	public Segment(){

	}



	@Override
	public String toString() {
		return "Segment [from=" + from + ", to=" + to + ", score=" + score
				+ "]";
	}



	public Integer getFrom() {
		return from;
	}

	public void setFrom(Integer from) {
		this.from = from;
	}

	public Integer getTo() {
		return to;
	}

	public void setTo(Integer to) {
		this.to = to;
	}

	public double getScore() {
		return score;
	}

	public void setScore(double score) {
		this.score = score;
	}



	@Override
	public int compareTo(Segment o) {

		Integer s1 = getFrom();
		Integer s2 = o.getFrom();

		int comp = s1.compareTo(s2);
		if ( comp != 0)
			return comp;

		Integer e1 = getTo();
		Integer e2 = o.getTo();

		return e1.compareTo(e2);

	}


}

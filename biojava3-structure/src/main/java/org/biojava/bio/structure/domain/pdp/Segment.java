package org.biojava.bio.structure.domain.pdp;

import java.io.Serializable;




public class Segment implements Serializable, Comparable<Segment> {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1393487067559539657L;
	private int from;
	private int to;
	double score;

	public Segment(){

	}
	
	

	@Override
	public String toString() {
		return "Segment [from=" + from + ", to=" + to + ", score=" + score
				+ "]";
	}


	
	public int getFrom() {
		return from;
	}

	public void setFrom(int from) {
		this.from = from;
	}

	public int getTo() {
		return to;
	}

	public void setTo(int to) {
		this.to = to;
	}

	public double getScore() {
		return score;
	}

	public void setScore(double score) {
		this.score = score;
	}



	public int compareTo(Segment o) {
		int s1 = getFrom();
		int s2 = o.getFrom();
		
		if (s1 > s2)
			return 1;
		if ( s1 == s2)
			return 0;
		return -1;
			
	}


}

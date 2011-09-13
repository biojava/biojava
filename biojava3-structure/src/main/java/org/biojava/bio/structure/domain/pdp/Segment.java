package org.biojava.bio.structure.domain.pdp;

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

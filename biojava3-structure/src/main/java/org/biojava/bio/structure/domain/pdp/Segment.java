package org.biojava.bio.structure.domain.pdp;

public class Segment {
	public int from;
	public int to;
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


}

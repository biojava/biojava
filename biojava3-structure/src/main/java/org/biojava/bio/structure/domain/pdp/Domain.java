package org.biojava.bio.structure.domain.pdp;

import java.io.Serializable;
import java.util.ArrayList;

import java.util.List;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlRootElement;

@XmlRootElement(name = "Domain", namespace ="http://www.biojava.org")
@XmlAccessorType(XmlAccessType.PUBLIC_MEMBER)

public class Domain implements Comparable<Domain>, Serializable{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = -1293994033102271366L;
	
	String id;
	int size;
	int nseg;
	double score;

	List<Segment>segments = new ArrayList<Segment>();

	public Domain(){
		
	}
	
	public String getId() {
		return id;
	}



	public void setId(String id) {
		this.id = id;
	}



	@Override
	public String toString() {
		return "Domain [size=" + size + ", nseg=" + nseg + ", score=" + score
				
				+ "]";
	}
	
	public List<Segment> getSegments() {
		
		return segments;
	}
	
	public Segment getSegmentAtPos(int pos){
		int size = segments.size();
		while ( pos >= size){
			segments.add(new Segment());
			size++;
			
		}
		return segments.get(pos);
	}


	public int compareTo(Domain other) {
		
		
		if ( this.getId() == null)
			return 1;
		if ( other.getId() == null)
			return -1;
		return this.getId().compareTo(other.getId());
	}

	public int getSize() {
		return size;
	}

	public void setSize(int size) {
		this.size = size;
	}

	public int getNseg() {
		return nseg;
	}

	public void setNseg(int nseg) {
		this.nseg = nseg;
	}

	public double getScore() {
		return score;
	}

	public void setScore(double score) {
		this.score = score;
	}

	public void setSegments(List<Segment> segments) {
		this.segments = segments;
	}
	
	
	
}



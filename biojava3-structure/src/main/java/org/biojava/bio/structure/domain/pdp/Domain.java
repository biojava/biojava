package org.biojava.bio.structure.domain.pdp;

import java.io.Serializable;
import java.util.ArrayList;

import java.util.List;
public class Domain implements Comparable<Domain>, Serializable{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = -1293994033102271366L;
	
	String id;
	public int size;
	public int nseg;
	public double score;

	public List<Segment>segments = new ArrayList<Segment>();

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
	
	public Segment getSegment(int pos){
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
	
	
}



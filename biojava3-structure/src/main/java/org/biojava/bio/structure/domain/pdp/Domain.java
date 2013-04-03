package org.biojava.bio.structure.domain.pdp;

import java.util.ArrayList;

import java.util.List;

public class Domain {
	
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
	
	
}



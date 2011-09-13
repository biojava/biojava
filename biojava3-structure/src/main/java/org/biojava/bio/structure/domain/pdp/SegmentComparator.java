package org.biojava.bio.structure.domain.pdp;

import java.util.Comparator;

public class SegmentComparator implements Comparator<Segment> {

	public int compare(Segment v1, Segment v2) {
		
		return v1.compareTo(v2);
	}	
}



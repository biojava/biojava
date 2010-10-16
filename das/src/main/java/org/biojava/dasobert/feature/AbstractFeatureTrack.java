/*
 *                  BioJava development code
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
 * Created on Feb 9, 2005
 *
 */
package org.biojava.dasobert.feature;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/** An Abstract class representing a Feature as being diplayed in the SeqFeaturePanel 
 *  A feature corresponds to everything that is visible in a "line" and can contain one or multiple Segments.
 * 
 * 
 * @author Andreas Prlic
 *
 */
public abstract class AbstractFeatureTrack implements FeatureTrack,Cloneable {

	String name;
	String method;
	String type;
	List <Segment>  segments ;
	String note;
	String link;
	String source;
	String score;
	String orientation;
	String typeID;
	String typeCategory;

	public AbstractFeatureTrack() {
		source = "Unknown";
		method = "Unknown";
		type   = "Unknown";
		note   = "";
		link   = "";
		score  = "";
		orientation = null;
		segments = new ArrayList<Segment>();

	}

	public abstract Object clone();


	public String toString() {
		String str = "Feature: method: " + method +" type: " + type ;
		if ( name != null)
			str += " name: " + name;
		
		if (( note != null) && (! note.equals("null")))
		{
			if (note.length() > 40)
				str += "note: " +note.substring(0,39) + "...";
			else
				str += " note: "+note;
		}
		str += " # segments: " + segments.size() ;
		return str ;
	}


	/** returns true if the specified sequence position is within the range of this Feature
	 * 
	 * @param seqPosition the position to check
	 * @return true if the position is within the ranges of the segments of this feature
	 */
	public boolean overlaps(int seqPosition){
		List segments = getSegments();
		Iterator iter =segments.iterator();

		while (iter.hasNext()){

			Segment seg = (Segment) iter.next();
			if ( seg.overlaps(seqPosition) )
				return true;                                     
		}

		return false;
	}

	public void setSource(String s) { source = s;}
	public String getSource() { return source; };


	public void setName(String nam) { name = nam; }
	public String getName() { return name; }

	public void setMethod(String methd) { method = methd ; }
	public String getMethod() { return method ; }

	public void setType(String typ) { type = typ ; }
	public String getType() { return type ; }

	public void setNote(String nte) { 
		if (nte != null)
			note = nte; 
		}
	public String getNote() { return note ; }

	public void setLink(String lnk) { link = lnk;}
	public String getLink() { return link;}

	public void setScore(String s){ score = s;}
	public String getScore() { return score;}

	/** add a segment to this feature */
	public void addSegment(int start, int end, String name) {
		Segment s = new SegmentImpl() ;
		s.setStart(start);
		s.setEnd(end) ;
		s.setName(name);
		s.setParent(this);
		segments.add(s);
	}

	public void addSegment( Segment s ){
		s.setParent(this);
		segments.add(s);
	}

	public List<Segment> getSegments() {
		return segments ;
	}


	public String getOrientation() {
		return orientation;
	}

	public void setOrientation(String orientation) {
		this.orientation = orientation;
	}

	/** test if two features are equivalent 
	 * important: only comares type,method and source.
	 * The individual segments are not compared!
	 * 
	 * */
	public  boolean equals(FeatureTrack feat) {
//		if ( note == null) {
		//  if (( feat.getNote() == null ) || 
		// ( feat.getNote().equals(""))) {
		//} else if ( this.note.equals(feat.getNote())){
		//  return true;
		//}
		if ( this.type.equals(feat.getType())){
			if ( this.method.equals(feat.getMethod())){
				if ( this.source.equals(feat.getSource())){										
					if (this.note.equals(feat.getNote())){
						return true;
					}
				}
			}
		}
		return false;

	}
	public String getTypeCategory() {
		// TODO Auto-generated method stub
		return typeCategory;
	}

	public String getTypeID() {
		// TODO Auto-generated method stub
		return typeID;
	}

	public void setTypeCategory(String typeCategory) {
		this.typeCategory = typeCategory;
		
	}

	public void setTypeID(String typeID) {
		this.typeID = typeID;
		
	}
	


}

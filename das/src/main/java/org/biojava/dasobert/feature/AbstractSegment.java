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
 * Created on May 22, 2007
 * 
 */

package org.biojava.dasobert.feature;

import java.awt.Color;

public abstract class AbstractSegment implements Segment,Cloneable {
	int start   ;
	int end     ;
	String name ;
	Color color ;
	FeatureTrack parent ;
	String txtColor ;
	String note;

	public abstract Object clone();
	
	/* (non-Javadoc)
	 * @see org.biojava.spice.feature.SegmentIF#toString()
	 */
	public String toString() {
		String str = "Segment: " +name + " " +start + " " + end ;
		if ( ( note != null ) && ( ! note.equals("null")))
			if ( note.length() >40)
				str += note.substring(0,39)+"...";
			else
				str += note;
		return str ;
	}

	/* (non-Javadoc)
	 * @see org.biojava.spice.feature.SegmentIF#getNote()
	 */
	public String getNote() {
		return note;
	}

	/* (non-Javadoc)
	 * @see org.biojava.spice.feature.SegmentIF#setNote(java.lang.String)
	 */
	public void setNote(String note) {
		this.note = note;
	}

	/* (non-Javadoc)
	 * @see org.biojava.spice.feature.SegmentIF#setStart(int)
	 */
	public void setStart(int strt) {start = strt ; }
	/* (non-Javadoc)
	 * @see org.biojava.spice.feature.SegmentIF#getStart()
	 */
	public int  getStart() {return start ;}

	/* (non-Javadoc)
	 * @see org.biojava.spice.feature.SegmentIF#setEnd(int)
	 */
	public void setEnd(int ed) { end = ed;}
	/* (non-Javadoc)
	 * @see org.biojava.spice.feature.SegmentIF#getEnd()
	 */
	public int getEnd() { return end;}

	/* (non-Javadoc)
	 * @see org.biojava.spice.feature.SegmentIF#setName(java.lang.String)
	 */
	public void setName(String nam) { name = nam;}
	/* (non-Javadoc)
	 * @see org.biojava.spice.feature.SegmentIF#getName()
	 */
	public String getName() { return name ; }

	/* (non-Javadoc)
	 * @see org.biojava.spice.feature.SegmentIF#setColor(java.awt.Color)
	 */
	public void setColor(Color col) { color = col; }
	/* (non-Javadoc)
	 * @see org.biojava.spice.feature.SegmentIF#getColor()
	 */
	public Color getColor() { return color ; }

	/* (non-Javadoc)
	 * @see org.biojava.spice.feature.SegmentIF#setParent(org.biojava.spice.feature.Feature)
	 */
	public void setParent(FeatureTrack f) { parent = f;}
	/* (non-Javadoc)
	 * @see org.biojava.spice.feature.SegmentIF#getParent()
	 */
	public FeatureTrack getParent(){ return parent;}

	/* (non-Javadoc)
	 * @see org.biojava.spice.feature.SegmentIF#setTxtColor(java.lang.String)
	 */
	public void setTxtColor(String str) { txtColor = str; }
	/* (non-Javadoc)
	 * @see org.biojava.spice.feature.SegmentIF#getTxtColor()
	 */
	public String getTxtColor() { return txtColor;}


	/* (non-Javadoc)
	 * @see org.biojava.spice.feature.SegmentIF#overlaps(int)
	 */
	public boolean overlaps(int seqPosition){
		if ( ( getStart() <= seqPosition) && ( getEnd() >= seqPosition)){
			return true;             
		}   
		return false;
	}



	/* (non-Javadoc)
	 * @see org.biojava.spice.feature.SegmentIF#overlaps(org.biojava.spice.feature.Segment)
	 */
	public boolean overlaps(Segment segment){
		if (! (this.start <= this.end )) 
			throw new IndexOutOfBoundsException("start > end for segment" + this);

		if ( ! (segment.getStart() <= segment.getEnd() ))
			throw new IndexOutOfBoundsException("start > end for segment" + segment);

		// start must be in region of other
		if ( this.start >= segment.getStart()){
			if ( this.start <= segment.getEnd()){
				return true;
			}
		}
		// or end must be in region of other..
		if ( this.end >= segment.getStart() ) {
			if ( this.end <= segment.getEnd()){
				return true;
			}
		}

		if ( this.start <= segment.getStart() ) {
			if ( this.end >= segment.getEnd() ) {
				return true;
			}
		}
		return false;
	}
}

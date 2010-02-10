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

/** a Histogram segment is an extension of the standard Segment with a score
 * 
 * @author Andreas Prlic
 *
 */ 
public class HistogramSegment extends AbstractSegment {

	double score;
	
	public HistogramSegment() {
		super();
		
	}

	public double getScore() {
		return score;
	}

	public void setScore(double score) {
		this.score = score;
	}
	
	
	public Object clone(){
        
        Segment s = new HistogramSegment();
        s.setStart(start);
        s.setEnd(end);
        s.setName(name);
        s.setColor(color);
        s.setTxtColor(txtColor);
        s.setNote(note);
        return s;
        
    }
	
	
	
	
}

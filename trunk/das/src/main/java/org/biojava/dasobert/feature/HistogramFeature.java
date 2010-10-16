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


import java.util.Iterator;

/** a class that represents Histogram Style features
 * in addition to normal features they know about Max and Minimum scores for the whole line
 * Histogram feautes have only one (Histogram) Segment, which contains the scores for each position
 *  
 * @author Andreas Prlic
 *
 */
public class HistogramFeature 
	extends AbstractFeatureTrack {
	
	double max;
	double min;
	
	public HistogramFeature() {
		super();
		// TODO Auto-generated constructor stub
	}
	
	
	
	public double getMax() {
		return max;
	}



	public void setMax(double max) {
		this.max = max;
	}



	public double getMin() {
		return min;
	}



	public void setMin(double min) {
		this.min = min;
	}



	public Object clone() {
		
		HistogramFeature f = new HistogramFeature();
		
		f.setName(name);
    	f.setMethod(method);
    	f.setType(type);
    	f.setNote(note);
    	f.setLink(link);
    	f.setSource(source);
    	f.setScore(score);
    	
    	Iterator iter = segments.iterator();
    	
    	while (iter.hasNext()){
    		Segment s = (Segment) iter.next();
    		f.addSegment((Segment)s.clone());
    	}
    	
    	return f;
		
		
	}
	
	
	
	
	

}

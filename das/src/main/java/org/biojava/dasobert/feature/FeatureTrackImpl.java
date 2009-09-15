/*
 *                    BioJava development code
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
 * Created on 22.09.2004
 * @author Andreas Prlic
 *
 */

package org.biojava.dasobert.feature;

import java.util.Iterator;

/** a class to store FeatureData and to visualize them
 * coordinate system of features is always UniProt !
 * PDBresnum features served by DAS need to be converted into UniProt coord sys first.
 *
 * a feature consists of one or several segments.
 * segmetns cotnains <start> and <end> information.
 *
 * @author Andreas Prlic
 */
public class FeatureTrackImpl 
extends AbstractFeatureTrack
implements FeatureTrack

{
   
	
	 public Object clone(){
	    	
	    	FeatureTrack f = new FeatureTrackImpl();
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





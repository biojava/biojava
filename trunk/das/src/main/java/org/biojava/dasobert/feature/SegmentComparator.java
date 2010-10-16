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

import java.util.Comparator;

public class SegmentComparator implements Comparator {

	public int compare(Object arg0, Object arg1) {
		
		Segment s1 = (Segment) arg0;
		Segment s2 = (Segment) arg1;
		
		if (s1.getStart() < s2.getStart())
			return -1;
		if (s1.getStart() > s2.getStart())
			return 1;		
		
		return 0;
	}

}

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
 * Created on 23.09.2004
 * @author Andreas Prlic
 *
 */


package org.biojava.dasobert.feature ;

import java.util.Comparator ;
import java.util.Iterator;
import java.util.List;


/** a comparator to sort Features byt type
 * @author Andreas Prlic
 */

public class FeatureComparator 
implements Comparator
{

	public FeatureComparator() {
	}

	public int compare(Object a, Object b) {
		FeatureTrack x = (FeatureTrack) a;
		FeatureTrack y = (FeatureTrack) b;

		String typea = x.getType();
		String typeb = y.getType();

		if ( ! typea.equals(typeb))
			return typea.compareTo(typeb);

		List s1 = x.getSegments();
		List s2 = y.getSegments();

		Iterator iter1 = s1.iterator();
		Iterator iter2 = s2.iterator();

		while (iter1.hasNext()){
			Segment seg1 = (Segment)iter1.next();
			int start1 = seg1.getStart();

			while (iter2.hasNext()){
				Segment seg2 = (Segment)iter2.next();
				int start2 = seg2.getStart();

				if ( start1 < start2){
					return -1;
				} if ( start1 > start2){
					return 1;
				}

			}
		}

		return 0;

	}

}

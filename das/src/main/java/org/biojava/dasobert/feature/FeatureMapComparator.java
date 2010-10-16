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
import java.util.Map ;

/** a comparator to sort Features if they are still in a Map ( sorts by type )
 * @author Andreas Prlic
 */

public class FeatureMapComparator 
    implements Comparator
{

    public FeatureMapComparator() {
    }

    public int compare(Object a, Object b) {
	Map x = (Map) a;
	Map y = (Map) b;

	String typea = (String)x.get("TYPE");
	String typeb = (String)y.get("TYPE");
	
	

	if ( isSecstruc(typea) && isSecstruc(typeb)) {
	    return 0 ;
	}
	return typea.compareTo(typeb);
    }

    public boolean isSecstruc(String type) {
	if ( type.equals("HELIX") 
	     ||
	     type.equals("STRAND") 
	     ||
	     type.equals("TURN") 
	     ) {
	    return true ;
	}
	return false ;
    }

}

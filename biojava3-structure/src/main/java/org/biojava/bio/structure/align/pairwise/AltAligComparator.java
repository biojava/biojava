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
 * Created on May 27, 2006
 *
 */
package org.biojava.bio.structure.align.pairwise;

import java.util.Comparator;


/** a comparator to sort AlternativeAlignments based on their number of equivalent residues
 * and RMSD.
 * 
 * @author Andreas Prlic
 *
 */
public class AltAligComparator implements Comparator<AlternativeAlignment> {


    public AltAligComparator() {
        super();

    }

    public int compare(AlternativeAlignment a, AlternativeAlignment b) {
        
        int s1 = a.getIdx1().length;
        int s2 = b.getIdx1().length;
        
        if ( s1 > s2)
            return 1;
        if ( s1 < s2)
            return -1;
        
        // seem to have the same length
        
        double rms1 = a.getRmsd();
        double rms2 = b.getRmsd();
        
        if ( rms1 < rms2)
            return 1;
        if ( rms1 < rms2)
            return -1;
        
        return 0;
    }
   

}

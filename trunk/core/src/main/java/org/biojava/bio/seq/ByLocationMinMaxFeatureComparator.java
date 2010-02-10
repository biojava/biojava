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
 */
package org.biojava.bio.seq;



/**
 * Comparator that compares the min and max positions of Features
 * 
 * Required by org.biojava.bio.gui.sequence.AbstractPeptideDigestRenderer instances.
 *
 * @author Mark Southern
 * @since 1.5
 */
public class ByLocationMinMaxFeatureComparator implements java.util.Comparator {
    private ByLocationMinMaxComparator lc = new ByLocationMinMaxComparator();

    public ByLocationMinMaxFeatureComparator() {
    }

    public int compare(Object o1, Object o2) {
        Feature f1 = ( Feature ) o1;
        Feature f2 = ( Feature ) o2;

        return lc.compare(f1.getLocation(), f2.getLocation());
    }
}

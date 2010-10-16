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

import java.util.Comparator;

import org.biojava.bio.symbol.Location;


/**
 * A Comparator similar to Feature.ByLocationComparator except that the min and max positions of
 * the location are both compared
 * 
 *
 * @author Mark Southern
 * @since 1.5
 */
public class ByLocationMinMaxComparator implements Comparator {
    private static Comparator comparator = new ByLocationMinMaxComparator();

    public ByLocationMinMaxComparator() {
    }

    public static Comparator comparator() {
        return comparator;
    }

    public int compare(Object o1, Object o2) {
        Location loc1 = ( Location ) o1;
        Location loc2 = ( Location ) o2;

        Integer i1 = new Integer(loc1.getMin());
        Integer i2 = new Integer(loc2.getMin());
        int iComp = i1.compareTo(i2);

        if (iComp != 0) {
            return iComp;
        }

        i1 = new Integer(loc1.getMax());
        i2 = new Integer(loc2.getMax());

        return i1.compareTo(i2);
    }
}

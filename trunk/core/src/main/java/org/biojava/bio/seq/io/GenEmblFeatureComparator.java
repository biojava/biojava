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

package org.biojava.bio.seq.io;

import java.util.Comparator;

import org.biojava.bio.seq.Feature;

/**
 * <code>GenEmblFeatureComparator</code> compares
 * <code>Feature</code>s by the minimum base position of their
 * <code>Location</code>, but places <code>Feature</code>s with
 * type "source" first.
 * @deprecated Use org.biojavax.bio.seq.io framework instead
 * @author Keith James
 * @since 1.2
 */
public final class GenEmblFeatureComparator implements Comparator
{
    public static final Comparator INSTANCE = new GenEmblFeatureComparator();

    private GenEmblFeatureComparator() { }

    public int compare(Object o1, Object o2)
    {
        Feature f1 = (Feature) o1;
        Feature f2 = (Feature) o2;

        boolean source1 = f1.getType().equals("source") ? true : false;
        boolean source2 = f2.getType().equals("source") ? true : false;

        // We don't subtract one coordinate from another as one or
        // both may be set to Integer.MAX_VALUE/Integer.MIN_VALUE
        // and the result could wrap around. Convert to Long if
        // necessary.
        if (! source1 && source2)
        {
            return 1;
        }
        else if (source1 && ! source2)
        {
            return -1;
        }
        else if (f1.getLocation().getMin() > f2.getLocation().getMin())
            return 1;
        else if (f1.getLocation().getMin() < f2.getLocation().getMin())
            return -1;
        else
            return 0;
    }
}

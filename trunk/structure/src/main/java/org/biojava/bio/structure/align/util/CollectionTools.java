package org.biojava.bio.structure.align.util;

import java.util.*;

/**
 * Utilities for working with collections.
 *
 * @author Thomas Down
 */
public class CollectionTools {
    public static int[] toIntArray(Collection<Number> l) {
        int[] a = new int[l.size()];
        int i = 0;
        for (Iterator<Number> j = l.iterator(); j.hasNext(); ) {
            a[i++] = ( j.next()).intValue();
        }
        return a;
    }
    
    public static double[] toDoubleArray(Collection<Number> l) {
        double[] a = new double[l.size()];
        int i = 0;
        for (Iterator<Number> j = l.iterator(); j.hasNext(); ) {
            a[i++] = ((Number) j.next()).doubleValue();
        }
        return a;
    }
    
    public static Object randomPick(Collection<Number> col) {
        Object[] objs = col.toArray(new Object[col.size()]);
        return objs[(int) Math.floor(Math.random() * objs.length)];
    }
}

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
package org.biojava.nbio.structure.align.util;

import java.util.Collection;
import java.util.Iterator;

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
			a[i++] = j.next().doubleValue();
		}
		return a;
	}

	public static Object randomPick(Collection<Number> col) {
		Object[] objs = col.toArray(new Object[col.size()]);
		return objs[(int) Math.floor(Math.random() * objs.length)];
	}
}

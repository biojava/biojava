/*
 *                    PDB web development code
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
 *
 * Created on May 23, 2009
 * Created by Andreas Prlic
 *
 */

package org.biojava.nbio.structure.align.helper;

import java.io.Serializable;
import java.util.Comparator;

public class IdxComparator implements Comparator<int[]>, Serializable
{
    private static final long serialVersionUID = 1;

	@Override
	public int compare(int[] o1, int[] o2)
	{
		if (((o1[0]) == (o2[0])) &&
				((o2[1]) == (o2[1])))
			return 0;
		if ( o1[0] < o2[0])
			return -1;

		return 1;
	}

}

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
package org.biojava.nbio.genome;

import static org.junit.Assert.*;

import org.biojava.nbio.genome.parsers.gff.Location;
import org.junit.Test;

public class TestIssue355 {

	@Test
	public void testIssue1() {
		Location l1 = Location.fromBio(51227320, 51227381, '+');
		Location l2 = Location.fromBio(51227323, 51227382, '+');

		Location union = l1.union(l2);
		assertEquals(51227320,union.bioStart());
		assertEquals(51227382,union.bioEnd());
	}
	
	@Test
	public void testIssue2() {
		Location l1 = Location.fromBio(100, 200, '+');
		Location l2 = Location.fromBio(1, 99, '+');
		Location intersection = l1.intersection(l2);
		assertNull(intersection);
	}

}


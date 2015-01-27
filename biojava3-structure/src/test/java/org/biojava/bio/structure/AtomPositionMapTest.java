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
 * Created on 2012-12-01
 *
 */

package org.biojava.bio.structure;

import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.LocalPDBDirectory;
import org.junit.Before;
import org.junit.Test;

import java.io.IOException;
import java.util.NavigableMap;

import static org.junit.Assert.assertEquals;

/**
 * A unit test for {@link org.biojava.bio.structure.AtomPositionMap}.
 * @author dmyerstu
 */
public class AtomPositionMapTest {

	@Before
	public void setUp() throws Exception {
		cache = new AtomCache(); // TODO Should mock instead of depending on real data from AtomCache
		cache.setObsoleteBehavior(LocalPDBDirectory.ObsoleteBehavior.FETCH_OBSOLETE);
	}
	
	private AtomCache cache;
	
	/**
	 * Tests with no insertion codes.
	 * @throws StructureException 
	 * @throws IOException 
	 */
	@Test
	public void testEasy() throws IOException, StructureException { // no insertion codes
		String pdbId = "1w0p";
		int length = 92;
		ResidueNumber start = new ResidueNumber("A", 25, null);
		ResidueNumber end = new ResidueNumber("A", 117, null);
		AtomPositionMap map = new AtomPositionMap(cache.getAtoms(pdbId));
		NavigableMap<ResidueNumber,Integer> navMap = map.getNavMap();
		for (ResidueNumber n : navMap.keySet()) {
			assertEquals("An element is missing", map.getPosition(n).intValue(), navMap.get(n).intValue());
		}
		int realLength = map.calcLength(start, end);
		assertEquals("Real atom length is wrong", length, realLength);
	}

	/**
	 * Tests with insertion codes.
	 * @throws StructureException 
	 * @throws IOException 
	 */
	@Test
	public void testHard() throws IOException, StructureException {

		String pdbId = "1qdm";
		// has 2 insertion code regions:
		// P at the beginning starting at 6P and ending at 27P, where 2 starts
		// S between 247 (before:1S) and 248 (after:104S)

		AtomPositionMap map = new AtomPositionMap(cache.getAtoms(pdbId));
		NavigableMap<ResidueNumber,Integer> navMap = map.getNavMap();

		for (ResidueNumber n : navMap.keySet()) {
			assertEquals("An element is missing", map.getPosition(n).intValue(), navMap.get(n).intValue());
		}
		
		int length1 = 59;
		int length2 = 131;
		ResidueNumber start = new ResidueNumber("A", 246, null);
		ResidueNumber mid = new ResidueNumber("A", 85, 'S');
		ResidueNumber end = new ResidueNumber("A", 300, null);
		int realLength1 = map.calcLength(start, mid);
		assertEquals("Real atom length is wrong", length1, realLength1);
		int realLength2 = map.calcLength(start, end);
		assertEquals("Real atom length is wrong", length2, realLength2);
	}

}

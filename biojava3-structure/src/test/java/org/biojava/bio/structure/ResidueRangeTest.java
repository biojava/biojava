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
 * Created on 2012-11-20
 *
 */

package org.biojava.bio.structure;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.NavigableMap;

import org.biojava.bio.structure.align.util.AtomCache;
import org.junit.Before;
import org.junit.Test;

/**
 * A unit test for {@link ResidueRange}.
 * @author dmyerstu
 *
 */
public class ResidueRangeTest {

	private AtomCache cache;
	
	@Before
	public void setUp() throws Exception {
		cache = new AtomCache();
	}
	
	/**
	 * Tests creating ResidueRanges and calculating their lengths.
	 */
	@Test
	public void testBasic() {
		String[] ids = new String[] {"1w0p", "3qq3", "3chc", "2ei7"}; // more: , "2qbr"
		char[] chains = new char[] {'A', 'B', 'A', 'L'};
		ResidueNumber[] starts = new ResidueNumber[] {new ResidueNumber("A", 5, ' '), new ResidueNumber("B", 10, 's'), new ResidueNumber("A", 15, 'm'), new ResidueNumber("L", 44, ' ')};
		ResidueNumber[] ends = new ResidueNumber[] {new ResidueNumber("A", 117, ' '), new ResidueNumber("B", 200, 's'), new ResidueNumber("A", 464, 'q'), new ResidueNumber("L", 254, 't')};
		Integer[] lengths = new Integer[] {117-5, 200-10, 111, null};
		int totalLength = 0;
		List<ResidueRange> ranges = new ArrayList<ResidueRange>(ids.length);
		for (int i = 0; i < ids.length; i++) {
			ResidueRange rr = new ResidueRange(chains[i], starts[i], ends[i], lengths[i]);
			assertEquals("The chain is incorrect", chains[i], rr.getChain());
			assertEquals("The start is incorrect", starts[i], rr.getStart());
			assertEquals("The end is incorrect", ends[i], rr.getEnd());
			assertEquals("The length is incorrect", lengths[i], rr.getLength());
			ranges.add(rr);
			if (lengths[i] != null) {
				totalLength += lengths[i];
				assertEquals("Total length is wrong", totalLength, ResidueRange.calcLength(ranges));
			} else {
				try {
					ResidueRange.calcLength(ranges); // should fail
					fail("Lengths should be undefined");
				} catch (IllegalArgumentException e) {}
			}
		}
	}
	
	/**
	 * Tests {@link ResidueRange#parseMultiple(String)}.
	 */
	@Test
	public void testParseAndEqual() {

		String pdbId1 = "2eke";
		String string1 = "C_1023-1063,C_1064-1084";
		List<ResidueRange> list1 = ResidueRange.parseMultiple(string1);
		assertEquals(new ResidueRange('C', new ResidueNumber("C", 1023, null), new ResidueNumber("C", 1063, null), null), list1.get(0));
		assertEquals(new ResidueRange('C', new ResidueNumber("C", 1064, null), new ResidueNumber("C", 1084, null), null), list1.get(1));
		assertEquals(null, list1.get(0).getLength());
		assertEquals(null, list1.get(1).getLength());

		String pdbId = "1qdm";
		String string2 = "A_3S-37S,A_65S-99S";
		List<ResidueRange> list2 = ResidueRange.parseMultiple(string2);
		assertEquals(new ResidueRange('A', new ResidueNumber("A", 3, 'S'), new ResidueNumber("A", 37, 'S'), null), list2.get(0));
		assertEquals(new ResidueRange('A', new ResidueNumber("A", 65, 'S'), new ResidueNumber("A", 99, 'S'), null), list2.get(1));
	}

	/**
	 * Tests {@link ResidueRange#parseMultiple(String, NavigableMap)}.
	 * @throws StructureException 
	 * @throws IOException 
	 */
	@Test
	public void testParseAndEqualWithLengths() throws IOException, StructureException {

		AtomPositionMap map;
		
		String pdbId1 = "2eke";
		map = new AtomPositionMap(cache.getAtoms(pdbId1)); // TODO this could probably be mocked
		String string1 = "C_1023-1063,C_1064-1084";
		List<ResidueRange> list1 = ResidueRange.parseMultiple(string1, map);
		assertEquals(new ResidueRange('C', new ResidueNumber("C", 1023, null), new ResidueNumber("C", 1063, null), null), list1.get(0));
		assertEquals(new ResidueRange('C', new ResidueNumber("C", 1064, null), new ResidueNumber("C", 1084, null), null), list1.get(1));
		assertEquals(1063-1023, list1.get(0).getLength().intValue()); // no insertion codes
		assertEquals(1084-1064, list1.get(1).getLength().intValue());

		list1 = ResidueRange.parseMultiple(string1, map);
		assertEquals(new ResidueRange('C', new ResidueNumber("C", 1023, null), new ResidueNumber("C", 1063, null), null), list1.get(0));
		assertEquals(new ResidueRange('C', new ResidueNumber("C", 1064, null), new ResidueNumber("C", 1084, null), null), list1.get(1));
		
	}
	
}

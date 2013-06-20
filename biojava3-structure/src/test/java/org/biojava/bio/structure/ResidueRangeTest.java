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
import java.util.Iterator;
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
	
	@Test
	public void testIterator() throws IOException, StructureException {
		String pdbId = "2eke";
		String[] expected = new String[] {"C_1023", "C_1024", "C_1025", "C_1026", "C_1027", "C_1028", "C_1029", "C_1030", "C_1031", "C_1032", "C_1033", "C_1034", "C_1035", "C_1036", "C_1037", "C_1038", "C_1039", "C_1040", "C_1041", "C_1042", "C_1043", "C_1044", "C_1045", "C_1046", "C_1047", "C_1048", "C_1049", "C_1050", "C_1051", "C_1052", "C_1053", "C_1054", "C_1055", "C_1056", "C_1057", "C_1058", "C_1059", "C_1060", "C_1061", "C_1062", "C_1063"};
		ResidueRange rr = ResidueRange.parse("C_1023-1063");
		AtomPositionMap map = new AtomPositionMap(cache.getAtoms(pdbId));
		Iterator<ResidueNumber> iter = rr.iterator(map);
		int i = 0;
		while (iter.hasNext()) {
			ResidueNumber rn = iter.next();
			assertEquals(expected[i], rn.printFull());
			i++;
		}
	}

	@Test
	public void testMultiIterator() throws IOException, StructureException {
		String pdbId = "1qdm";
		String[] expected = new String[] {"A_3S", "A_4S", "A_5S", "A_6S", "A_7S", "A_8S", "A_9S", "A_10S", "A_11S", "A_12S", "A_13S", "A_14S", "A_15S", "A_16S", "A_17S", "A_18S", "A_19S", "A_20S", "A_21S", "A_22S", "A_23S", "A_24S", "A_25S", "A_26S", "A_27S", "A_28S", "A_29S", "A_30S", "A_31S", "A_32S", "A_33S", "A_34S", "A_35S", "A_36S", "A_37S", "A_65S", "A_66S", "A_67S", "A_68S", "A_69S", "A_70S", "A_71S", "A_72S", "A_73S", "A_74S", "A_75S", "A_76S", "A_77S", "A_78S", "A_79S", "A_80S", "A_81S", "A_82S", "A_83S", "A_84S", "A_85S", "A_86S", "A_87S", "A_88S", "A_89S", "A_90S", "A_91S", "A_92S", "A_93S", "A_94S", "A_95S", "A_96S", "A_97S", "A_98S", "A_99S"};
		List<ResidueRange> rrs = ResidueRange.parseMultiple("A_3S-37S,A_65S-99S");
		AtomPositionMap map = new AtomPositionMap(cache.getAtoms(pdbId));
		Iterator<ResidueNumber> iter = ResidueRange.multiIterator(map, rrs);
		int i = 0;
		while (iter.hasNext()) {
			ResidueNumber rn = iter.next();
			assertEquals(expected[i], rn.printFull());
			i++;
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

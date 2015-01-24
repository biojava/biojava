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

import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.LocalPDBDirectory.ObsoleteBehavior;
import org.junit.Before;
import org.junit.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

/**
 * A unit test for {@link ResidueRange}.
 * @author dmyerstu
 *
 */
public class ResidueRangeTest {

	private AtomCache cache;
	
	@Before
	public void setUp() throws Exception {
		cache = new AtomCache(); // TODO Should mock instead of depending on real data from AtomCache
		cache.setObsoleteBehavior(ObsoleteBehavior.FETCH_OBSOLETE);
	}

	@Test
	public void testWholeChainBasic() {
		String range = "B:";
		ResidueRange rr = ResidueRange.parse(range);
		assertEquals("Wrong chain Id", "B", rr.getChainId());
		assertNull("Start residue should be null", rr.getStart());
		assertNull("End residue should be null", rr.getEnd());
	}

	@Test
	public void testPartialChainWithMap() throws IOException, StructureException {
		String pdbId = "1cph";
		AtomPositionMap map = new AtomPositionMap(cache.getAtoms(pdbId));
		String range = "B:1-30";
		ResidueRange rr = ResidueRangeAndLength.parse(range, map);
		ResidueNumber start = new ResidueNumber("B", 1, null);
		ResidueNumber end = new ResidueNumber("B", 30, null);
		assertEquals("Wrong chain Id", "B", rr.getChainId());
		assertEquals("Wrong start", start, rr.getStart());
		assertEquals("Wrong end", end, rr.getEnd());
	}

	@Test
	public void testWholeChainWithMap() throws IOException, StructureException {
		String pdbId = "1cph";
		AtomPositionMap map = new AtomPositionMap(cache.getAtoms(pdbId));
		String range = "B:";
		ResidueRange rr = ResidueRangeAndLength.parse(range, map);
		assertEquals("Wrong chain Id", "B", rr.getChainId());
		assertNull("Wrong start", rr.getStart());
		assertNull("Wrong end", rr.getEnd());
	}

	/**
	 * Tests creating ResidueRanges and calculating their lengths.
	 */
	@Test
	public void testWithLengths() throws IOException, StructureException {
		String[] ids = new String[] {"1w0p", "3qq3", "3chc", "2ei7"}; // more: , "2qbr"
		String[] chains = new String[] {"A", "B", "A", "L"};
		ResidueNumber[] starts = new ResidueNumber[] {new ResidueNumber("A", 5, ' '), new ResidueNumber("B", 10, 's'), new ResidueNumber("A", 15, 'm'), new ResidueNumber("L", 44, ' ')};
		ResidueNumber[] ends = new ResidueNumber[] {new ResidueNumber("A", 117, ' '), new ResidueNumber("B", 200, 's'), new ResidueNumber("A", 464, 'q'), new ResidueNumber("L", 254, 't')};
		int[] lengths = new int[] {117-5, 200-10, 111, 55};
		int totalLength = 0;
		List<ResidueRangeAndLength> ranges = new ArrayList<ResidueRangeAndLength>(ids.length);
		for (int i = 0; i < ids.length; i++) {
			ResidueRangeAndLength rr = new ResidueRangeAndLength(chains[i], starts[i], ends[i], lengths[i]);
			assertEquals("The chain is incorrect", chains[i], rr.getChainId());
			assertEquals("The start is incorrect", starts[i], rr.getStart());
			assertEquals("The end is incorrect", ends[i], rr.getEnd());
			assertEquals("The length is incorrect", lengths[i], rr.getLength());
			ranges.add(rr);
			totalLength += lengths[i];
			assertEquals("Total length is wrong", totalLength, ResidueRangeAndLength.calcLength(ranges));
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

		//String pdbId1 = "2eke";
		String string1 = "C_1023-1063,C_1064-1084";
		List<ResidueRange> list1 = ResidueRange.parseMultiple(string1);
		assertEquals(new ResidueRange("C", new ResidueNumber("C", 1023, null), new ResidueNumber("C", 1063, null)), list1.get(0));
		assertEquals(new ResidueRange("C", new ResidueNumber("C", 1064, null), new ResidueNumber("C", 1084, null)), list1.get(1));

		//String pdbId = "1qdm";
		String string2 = "A_3S-37S,A_65S-99S";
		List<ResidueRange> list2 = ResidueRange.parseMultiple(string2);
		assertEquals(new ResidueRange("A", new ResidueNumber("A", 3, 'S'), new ResidueNumber("A", 37, 'S')), list2.get(0));
		assertEquals(new ResidueRange("A", new ResidueNumber("A", 65, 'S'), new ResidueNumber("A", 99, 'S')), list2.get(1));
	}

	/**
	 * Tests {@link org.biojava.bio.structure.ResidueRangeAndLength#parseMultiple(String, org.biojava.bio.structure.AtomPositionMap)}.
	 * @throws StructureException 
	 * @throws IOException 
	 */
	@Test
	public void testParseAndEqualWithLengths() throws IOException, StructureException {

		AtomPositionMap map;
		
		String pdbId1 = "2eke";
		map = new AtomPositionMap(cache.getAtoms(pdbId1));
		String string1 = "C_1023-1063,C_1064-1084";
		List<ResidueRangeAndLength> list1 = ResidueRangeAndLength.parseMultiple(string1, map);
		assertEquals(new ResidueRangeAndLength("C", new ResidueNumber("C", 1023, null), new ResidueNumber("C", 1063, null), 1063-1023), list1.get(0));
		assertEquals(new ResidueRangeAndLength("C", new ResidueNumber("C", 1064, null), new ResidueNumber("C", 1084, null), 1084-1064), list1.get(1));

	}
	
	@Test
	public void testLooksLikeRange() {
		String[] yes = new String[] {"A_", "A:", "ABC:", "abc:", "A_5-100", "A_5-100S", "A_5S-100", "A_5S-100S", "A_-5-100", "A_-5--100", "A_-5S--100S", "ABC:-5--200S"};
		for (String s : yes) {
			assertTrue(s + " was not considered a valid range format", ResidueRange.looksLikeRange(s));
		}
		String[] no = new String[] {"A", "A1", "A_1", "A_1-", "A_1S-", "A_1-100-", "A_-10-1000_", "", "3A:1-100"};
		for (String s : no) {
			assertFalse(s + " was considered a valid range format", ResidueRange.looksLikeRange(s));
		}
	}
	
}

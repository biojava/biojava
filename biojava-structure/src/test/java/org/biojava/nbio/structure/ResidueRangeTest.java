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

package org.biojava.nbio.structure;

import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.LocalPDBDirectory.ObsoleteBehavior;
import org.junit.Before;
import org.junit.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import static org.junit.Assert.*;

/**
 * A unit test for {@link ResidueRange}.
 *
 * @author dmyerstu
 *
 */
public class ResidueRangeTest {

	private AtomCache cache;

	@Before
	public void setUp() throws Exception {
		cache = new AtomCache(); // TODO Should mock instead of depending on
		// real data from AtomCache
		cache.setObsoleteBehavior(ObsoleteBehavior.FETCH_OBSOLETE);
	}

	@Test
	public void testWholeChainBasic() {
		String range = "B:";
		ResidueRange rr = ResidueRange.parse(range);
		assertEquals("Wrong chain Id", "B", rr.getChainName());
		assertNull("Start residue should be null", rr.getStart());
		assertNull("End residue should be null", rr.getEnd());
	}

	@Test
	public void testPartialChainWithMap() throws IOException,
	StructureException {
		String pdbId = "1cph";
		AtomPositionMap map = new AtomPositionMap(cache.getAtoms(pdbId));
		String range = "B:1-30";
		ResidueRange rr = ResidueRangeAndLength.parse(range, map);
		ResidueNumber start = new ResidueNumber("B", 1, null);
		ResidueNumber end = new ResidueNumber("B", 30, null);
		assertEquals("Wrong chain Id", "B", rr.getChainName());
		assertEquals("Wrong start", start, rr.getStart());
		assertEquals("Wrong end", end, rr.getEnd());
	}

	@Test
	public void testWholeChainWithMap() throws IOException, StructureException {
		String pdbId = "1cph";
		AtomPositionMap map = new AtomPositionMap(cache.getAtoms(pdbId));
		String range = "B:";
		ResidueRange rr = ResidueRangeAndLength.parse(range, map);
		assertEquals("Wrong chain Id", "B", rr.getChainName());
		assertEquals("Wrong start", new ResidueNumber("B",1,null),rr.getStart());
		assertEquals("Wrong end", new ResidueNumber("B",30,null),rr.getEnd());
	}

	/**
	 * Tests creating ResidueRanges and calculating their lengths.
	 */
	@Test
	public void testWithLengths() throws IOException, StructureException {
		String[] ids = new String[] { "1w0p", "3qq3", "3chc", "2ei7" }; // more:
		// ,
		// "2qbr"
		String[] chains = new String[] { "A", "B", "A", "L" };
		ResidueNumber[] starts = new ResidueNumber[] {
				new ResidueNumber("A", 5, ' '),
				new ResidueNumber("B", 10, 's'),
				new ResidueNumber("A", 15, 'm'),
				new ResidueNumber("L", 44, ' ') };
		ResidueNumber[] ends = new ResidueNumber[] {
				new ResidueNumber("A", 117, ' '),
				new ResidueNumber("B", 200, 's'),
				new ResidueNumber("A", 464, 'q'),
				new ResidueNumber("L", 254, 't') };
		int[] lengths = new int[] { 117 - 5, 200 - 10, 111, 55 };
		int totalLength = 0;
		List<ResidueRangeAndLength> ranges = new ArrayList<ResidueRangeAndLength>(
				ids.length);
		for (int i = 0; i < ids.length; i++) {
			ResidueRangeAndLength rr = new ResidueRangeAndLength(chains[i],
					starts[i], ends[i], lengths[i]);
			assertEquals("The chain is incorrect", chains[i], rr.getChainName());
			assertEquals("The start is incorrect", starts[i], rr.getStart());
			assertEquals("The end is incorrect", ends[i], rr.getEnd());
			assertEquals("The length is incorrect", lengths[i], rr.getLength());
			ranges.add(rr);
			totalLength += lengths[i];
			assertEquals("Total length is wrong", totalLength,
					ResidueRangeAndLength.calcLength(ranges));
		}
	}

	@Test
	public void testLengths() throws IOException, StructureException {
		String pdbId = "1w0p";
		String rangeStr = "A:25-26";
		Atom[] atoms = cache.getAtoms(pdbId);
		AtomPositionMap map = new AtomPositionMap(atoms);
		ResidueRangeAndLength range = ResidueRangeAndLength
				.parse(rangeStr, map);
		assertEquals(2, range.getLength());
	}

	@Test
	public void testIterator() throws IOException, StructureException {
		String pdbId = "2eke";
		String[] expected = new String[] { "C_1023", "C_1024", "C_1025",
				"C_1026", "C_1027", "C_1028", "C_1029", "C_1030", "C_1031",
				"C_1032", "C_1033", "C_1034", "C_1035", "C_1036", "C_1037",
				"C_1038", "C_1039", "C_1040", "C_1041", "C_1042", "C_1043",
				"C_1044", "C_1045", "C_1046", "C_1047", "C_1048", "C_1049",
				"C_1050", "C_1051", "C_1052", "C_1053", "C_1054", "C_1055",
				"C_1056", "C_1057", "C_1058", "C_1059", "C_1060", "C_1061",
				"C_1062", "C_1063" };
		ResidueRange rr = ResidueRange.parse("C_1023-1063");
		AtomPositionMap map = new AtomPositionMap(cache.getAtoms(pdbId));
		Iterator<ResidueNumber> iter = rr.iterator(map);
		int i = 0;
		while (iter.hasNext()) {
			ResidueNumber rn = iter.next();
			assertTrue(i<expected.length);
			assertEquals(expected[i], rn.printFull());
			i++;
		}
		assertEquals(expected.length,i);
	}

	@Test
	public void testMultiIterator() throws IOException, StructureException {
		String pdbId = "1qdm";
		String[] expected = new String[] { "A_3S", "A_4S", "A_5S", "A_6S",
				"A_7S", "A_8S", "A_9S", "A_10S", "A_11S", "A_12S", "A_13S",
				"A_14S", "A_15S", "A_16S", "A_17S", "A_18S", "A_19S", "A_20S",
				"A_21S", "A_22S", "A_23S", "A_24S", "A_25S", "A_26S", "A_27S",
				"A_28S", "A_29S", "A_30S", "A_31S", "A_32S", "A_33S", "A_34S",
				"A_35S", "A_36S", "A_37S", "A_65S", "A_66S", "A_67S", "A_68S",
				"A_69S", "A_70S", "A_71S", "A_72S", "A_73S", "A_74S", "A_75S",
				"A_76S", "A_77S", "A_78S", "A_79S", "A_80S", "A_81S", "A_82S",
				"A_83S", "A_84S", "A_85S", "A_86S", "A_87S", "A_88S", "A_89S",
				"A_90S", "A_91S", "A_92S", "A_93S", "A_94S", "A_95S", "A_96S",
				"A_97S", "A_98S", "A_99S" };
		List<ResidueRange> rrs = ResidueRange
				.parseMultiple("A_3S-37S,A_65S-99S");
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
		String rangeStr;
		List<ResidueRange> ranges;
		ResidueRange range;

		// String pdbId1 = "2eke";
		rangeStr = "C_1023-1063,C_1064-1084";
		ranges = ResidueRange.parseMultiple(rangeStr);
		assertEquals(new ResidueRange("C", new ResidueNumber("C", 1023, null),
				new ResidueNumber("C", 1063, null)), ranges.get(0));
		assertEquals(new ResidueRange("C", new ResidueNumber("C", 1064, null),
				new ResidueNumber("C", 1084, null)), ranges.get(1));

		// String pdbId = "1qdm";
		rangeStr = "A_3S-37S,A_65S-99S";
		ranges = ResidueRange.parseMultiple(rangeStr);
		assertEquals(new ResidueRange("A", new ResidueNumber("A", 3, 'S'),
				new ResidueNumber("A", 37, 'S')), ranges.get(0));
		assertEquals(new ResidueRange("A", new ResidueNumber("A", 65, 'S'),
				new ResidueNumber("A", 99, 'S')), ranges.get(1));

		// Multi-character chains
		rangeStr = "AB,A1,ABCD_1-55,NotAG00dID:-5-1R";
		ranges = ResidueRange.parseMultiple(rangeStr);
		range = ranges.get(0);
		assertEquals("Error parsing " + rangeStr, "AB", range.getChainName());
		assertNull("Error parsing " + rangeStr, range.getStart());
		assertNull("Error parsing " + rangeStr, range.getEnd());
		range = ranges.get(1);
		assertEquals("Error parsing " + rangeStr, "A1", range.getChainName());
		assertNull("Error parsing " + rangeStr, range.getStart());
		assertNull("Error parsing " + rangeStr, range.getEnd());
		range = ranges.get(2);
		assertEquals("Error parsing " + rangeStr, "ABCD", range.getChainName());
		assertEquals("Error parsing " + rangeStr, new ResidueNumber("ABCD", 1,
				null), range.getStart());
		assertEquals("Error parsing " + rangeStr, new ResidueNumber("ABCD", 55,
				null), range.getEnd());
		range = ranges.get(3);
		assertEquals("Error parsing " + rangeStr, "NotAG00dID",
				range.getChainName());
		assertEquals("Error parsing " + rangeStr, new ResidueNumber(
				"NotAG00dID", -5, null), range.getStart());
		assertEquals("Error parsing " + rangeStr, new ResidueNumber(
				"NotAG00dID", 1, 'R'), range.getEnd());

		// Wildcard chains
		rangeStr = "_,__,_:1-5,_:+1-+5";
		ranges = ResidueRange.parseMultiple(rangeStr);
		range = ranges.get(0);
		assertEquals("Error parsing " + rangeStr, "_", range.getChainName());
		assertNull("Error parsing " + rangeStr, range.getStart());
		assertNull("Error parsing " + rangeStr, range.getEnd());
		range = ranges.get(1);
		assertEquals("Error parsing " + rangeStr, "_", range.getChainName());
		assertNull("Error parsing " + rangeStr, range.getStart());
		assertNull("Error parsing " + rangeStr, range.getEnd());
		range = ranges.get(2);
		assertEquals("Error parsing " + rangeStr, "_", range.getChainName());
		assertEquals("Error parsing " + rangeStr, new ResidueNumber("_", 1,
				null), range.getStart());
		assertEquals("Error parsing " + rangeStr, new ResidueNumber("_", 5,
				null), range.getEnd());
		range = ranges.get(3);
		assertEquals("Error parsing " + rangeStr, "_", range.getChainName());
		assertEquals("Error parsing " + rangeStr, new ResidueNumber("_", 1,
				null), range.getStart());
		assertEquals("Error parsing " + rangeStr, new ResidueNumber("_", 5,
				null), range.getEnd());

	}

	@Test(expected=IllegalArgumentException.class)
	public void testBadSyntax() throws IOException, StructureException {
		ResidueRange.parse("-");
	}

	@Test
	public void testPartialRange() throws IOException, StructureException {
		String rangeStr = "C_1023-";
		ResidueRange range = ResidueRange.parse(rangeStr);
		assertEquals(rangeStr,1023,(int)range.getStart().getSeqNum());
		assertNull(rangeStr,range.getEnd());

		rangeStr = "C_-";
		range = ResidueRange.parse(rangeStr);
		assertNull(rangeStr,range.getStart());
		assertNull(rangeStr,range.getEnd());
		
		rangeStr = "A_-+55";
		range = ResidueRange.parse(rangeStr);
		assertNull(rangeStr,range.getStart());
		assertEquals(rangeStr,55,(int)range.getEnd().getSeqNum());

	}
	@Test
	public void testPartialRangeLength() throws IOException, StructureException {
		AtomPositionMap map = new AtomPositionMap(cache.getAtoms("2eke"));
		String rangeStr = "C_1023-";
		ResidueRangeAndLength range = ResidueRangeAndLength.parse(rangeStr, map);
		
		assertEquals(rangeStr,1023,(int)range.getStart().getSeqNum());
		assertEquals(rangeStr,1095,(int)range.getEnd().getSeqNum());
		assertEquals(rangeStr, 73, range.getLength());
		
		
	}

	/**
	 * Tests
	 * {@link org.biojava.nbio.structure.ResidueRangeAndLength#parseMultiple(String, org.biojava.nbio.structure.AtomPositionMap)}
	 * .
	 *
	 * @throws StructureException
	 * @throws IOException
	 */
	@Test
	public void testParseAndEqualWithLengths() throws IOException,
	StructureException {
		String rangeStr;
		List<ResidueRangeAndLength> ranges;
		ResidueRangeAndLength range;

		AtomPositionMap map;

		String pdbId = "2eke";
		map = new AtomPositionMap(cache.getAtoms(pdbId));
		rangeStr = "C_1023-1063,C_1064-1084,C";//C is 105-112,1013-1095
		ranges = ResidueRangeAndLength.parseMultiple(rangeStr, map);
		assertEquals(new ResidueRangeAndLength("C", new ResidueNumber("C",
				1023, null), new ResidueNumber("C", 1063, null), 1063 - 1023+1),
				ranges.get(0));
		assertEquals(new ResidueRangeAndLength("C", new ResidueNumber("C",
				1064, null), new ResidueNumber("C", 1084, null), 1084 - 1064+1),
				ranges.get(1));
		assertEquals(new ResidueRangeAndLength("C", new ResidueNumber("C", 105,
				null), new ResidueNumber("C", 1095, null), 91), ranges.get(2));

		// Wildcard chains

		pdbId = "4r61"; // A:8-52,58-109,119-161
		map = new AtomPositionMap(cache.getAtoms(pdbId));
		rangeStr = "_,__,_:52-58";
		ranges = ResidueRangeAndLength.parseMultiple(rangeStr, map);
		range = ranges.get(0);
		assertEquals("Error parsing " + rangeStr, "A", range.getChainName());
		assertEquals("Error parsing " + rangeStr, new ResidueNumber("A",8,null),range.getStart());
		assertEquals("Error parsing " + rangeStr, new ResidueNumber("A",161,null),range.getEnd());
		range = ranges.get(1);
		assertEquals("Error parsing " + rangeStr, "A", range.getChainName());
		assertEquals("Error parsing " + rangeStr, new ResidueNumber("A",8,null),range.getStart());
		assertEquals("Error parsing " + rangeStr, new ResidueNumber("A",161,null),range.getEnd());
		range = ranges.get(2);
		assertEquals("Error parsing " + rangeStr, "A", range.getChainName());
		assertEquals("Error parsing " + rangeStr, new ResidueNumber("A", 52,null), range.getStart());
		assertEquals("Error parsing " + rangeStr, new ResidueNumber("A", 58,null), range.getEnd());

		// wildcards not converted without the map
		ResidueRange range2 = ResidueRange.parse("_");
		assertEquals("Error parsing " + rangeStr, "_", range2.getChainName());
		assertNull("Error parsing " + rangeStr,range2.getStart());
		assertNull("Error parsing " + rangeStr,range2.getEnd());

	}

	@Test
	public void testRangeRegex() {
		// Valid ranges
		String[] yes = new String[] { "A_", "A:", "ABC:", "abc:", "A_5-100",
				"A_5-100S", "A_5S-100", "A_5S-100S", "A_-5-100", "A_-5--100",
				"A_-5S--100S", "ABC:-5--200S", "A", "ABCD", "A_1",
				"A1", // valid multi-char chain name
				"3A:1-100", // Weird chain name
				"_", "_:1-10", "__-2--1", "__", // catch-all chain
				"A:-3-+1","A:-3-+1","A:+1-6", // Positive numbers ok, although weird
				"A:1-","A:1S-","A:--5","A:-+5", // Partial ranges
		};
		for (String s : yes) {
			assertTrue(s + " was not considered a valid range format",
					ResidueRange.RANGE_REGEX.matcher(s).matches());
		}
		// invalid ranges
		String[] no = new String[] {  "A_1-100-",
				 "", "-", "___", "__:","A_-10-1000_",
				
		};
		for (String s : no) {
			assertFalse(s + " was considered a valid range format",
					ResidueRange.RANGE_REGEX.matcher(s).matches());
		}
	}
	
	@Test
	public void testTerminalSymbols() {
		String rangeStr;
		ResidueRange range;
		
		rangeStr = "A:1-$";
		range = ResidueRange.parse(rangeStr);
		assertEquals(rangeStr,1,(int)range.getStart().getSeqNum());
		assertNull(rangeStr,range.getEnd());
		
		rangeStr = "A:^-1";
		range = ResidueRange.parse(rangeStr);
		assertNull(rangeStr,range.getStart());
		assertEquals(rangeStr,1,(int)range.getEnd().getSeqNum());
		
		rangeStr = "A:^-$";
		range = ResidueRange.parse(rangeStr);
		assertNull(rangeStr,range.getStart());
		assertNull(rangeStr,range.getEnd());
	}

}

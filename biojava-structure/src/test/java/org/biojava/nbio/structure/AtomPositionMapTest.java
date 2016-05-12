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

package org.biojava.nbio.structure;

import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.LocalPDBDirectory;
import org.junit.Before;
import org.junit.Test;

import java.io.IOException;
import java.util.NavigableMap;

import static org.junit.Assert.*;
import static org.junit.Assume.*;

/**
 * A unit test for {@link org.biojava.nbio.structure.AtomPositionMap}.
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
		// Straightforward case. Density for residues 25-777 (743 residues)
		String pdbId = "1w0p";
		int length = 93;
		ResidueNumber start = new ResidueNumber("A", 25, null);
		ResidueNumber end = new ResidueNumber("A", 117, null);
		AtomPositionMap map = new AtomPositionMap(cache.getAtoms(pdbId));
		NavigableMap<ResidueNumber,Integer> navMap = map.getNavMap();
		for (ResidueNumber n : navMap.keySet()) {
			assertEquals("An element is missing", map.getPosition(n).intValue(), navMap.get(n).intValue());
		}
		int realLength = map.getLength(start, end);
		assertEquals("Real atom length is wrong", length, realLength);
	}

	@Test
	public void testLengths() throws IOException, StructureException {
		// Two identical chains, residues 1-68, no insertion codes or missing residues
		String pdbId = "3w0e";

		Atom[] atoms = cache.getAtoms(pdbId);
		AtomPositionMap map = new AtomPositionMap(atoms);

		// Double check that the chain length is correct
		int chainAlen = cache.getStructure(pdbId).getPolyChainByPDB("A").getAtomGroups(GroupType.AMINOACID).size();
		assumeTrue(68==chainAlen);


		int len;
		int start,end;// 0-based

		// Single residue
		start = 0;
		end = 0;
		len = map.getLength(new ResidueNumber("A",start+1,null), new ResidueNumber("A",end+1,null));
		assertEquals("Bad length for ("+start+","+end+")",1, len);
		len = map.getLength(start, end,"A");
		assertEquals("Bad length for ("+start+","+end+")",1, len);
		len = map.getLengthDirectional(start, end, "A");
		assertEquals("Bad length for ("+start+","+end+")",1, len);
		len = map.getLengthDirectional(end, start, "A");
		assertEquals("Bad length for ("+start+","+end+")",1, len);

		// Short range
		start = 2;
		end = 4;
		len = map.getLength(new ResidueNumber("A",start+1,null), new ResidueNumber("A",end+1,null));
		assertEquals("Bad length for ("+start+","+end+")",3, len);
		len = map.getLength(start, end,"A");
		assertEquals("Bad length for ("+start+","+end+")",3, len);
		len = map.getLengthDirectional(start, end, "A");
		assertEquals("Bad length for ("+start+","+end+")",3, len);
		len = map.getLengthDirectional(end, start, "A");
		assertEquals("Bad length for ("+start+","+end+")",-3, len);

		//Full chain
		start = 0;
		end = chainAlen-1;
		len = map.getLength(new ResidueNumber("A",start+1,null), new ResidueNumber("A",end+1,null));
		assertEquals("Bad length for ("+start+","+end+")",chainAlen, len);
		len = map.getLength(start, end,"A");
		assertEquals("Bad length for ("+start+","+end+")",chainAlen, len);
		len = map.getLengthDirectional(start, end, "A");
		assertEquals("Bad length for ("+start+","+end+")",chainAlen, len);
		len = map.getLengthDirectional(end, start, "A");
		assertEquals("Bad length for ("+start+","+end+")",-chainAlen, len);

		// Chain spanning
		start = chainAlen-1;
		end = chainAlen;
		try {
			len = map.getLength(new ResidueNumber("A",chainAlen,null), new ResidueNumber("B",1,null));
			fail("Not the same chain");
		} catch( IllegalArgumentException e) {}
		len = map.getLength(start, end,"A");
		assertEquals("Bad length for ("+start+","+end+")",1, len);
		len = map.getLengthDirectional(start, end, "A");
		assertEquals("Bad length for ("+start+","+end+")",1, len);
		len = map.getLengthDirectional(end, start, "A");
		assertEquals("Bad length for ("+start+","+end+")",-1, len);
		len = map.getLengthDirectional(start,end, "B");
		assertEquals("Bad length for ("+start+","+end+")",1, len);

		start = chainAlen-2; //2 residues of A
		end = chainAlen+2; // 3 residues of B
		len = map.getLengthDirectional(start, end, "A");
		assertEquals("Bad length for ("+start+","+end+")",2, len);
		len = map.getLengthDirectional(end, start, "A");
		assertEquals("Bad length for ("+start+","+end+")",-2, len);
		len = map.getLengthDirectional(start, end, "B");
		assertEquals("Bad length for ("+start+","+end+")",3, len);
		len = map.getLengthDirectional(end, start, "B");
		assertEquals("Bad length for ("+start+","+end+")",-3, len);

		start = 0;
		end = chainAlen;
		try {
			len = map.getLength(new ResidueNumber("A",start+1,null), new ResidueNumber("A",end+1,null));
			fail("Residue found from the wrong chain");
		} catch( IllegalArgumentException e) {
			// end residue should be B1, not A142
		}
		// Chain Spanning
		len = map.getLength(start, end,"B");
		assertEquals("Bad length for ("+start+","+end+")",1, len);
		len = map.getLengthDirectional(start, end, "B");
		assertEquals("Bad length for ("+start+","+end+")",1, len);
		len = map.getLengthDirectional(end, start, "B");
		assertEquals("Bad length for ("+start+","+end+")",-1, len);
		len = map.getLength(start, end,"A");
		assertEquals("Bad length for ("+start+","+end+")",chainAlen, len);
		len = map.getLengthDirectional(start, end, "A");
		assertEquals("Bad length for ("+start+","+end+")",chainAlen, len);
		len = map.getLengthDirectional(end, start, "A");
		assertEquals("Bad length for ("+start+","+end+")",-chainAlen, len);
	}

	/**
	 * Tests with insertion codes.
	 * @throws StructureException
	 * @throws IOException
	 */
	@Test
	public void testInsertionCodes() throws IOException, StructureException {
		String pdbId = "1qdm";
		// has 2 insertion code regions, lettered P and S, as well as disordered regions:
		// 6P-26P,2-163,169-247,1S-37S,65S-104S,248-338
		// Len:21,  162,     79,    37,      40,     91 = 430

		AtomPositionMap map = new AtomPositionMap(cache.getAtoms(pdbId));
		NavigableMap<ResidueNumber,Integer> navMap = map.getNavMap();

		for (ResidueNumber n : navMap.keySet()) {
			assertEquals("An element is missing", map.getPosition(n).intValue(), navMap.get(n).intValue());
		}

		int length1 = 60; // 2+37+21
		int length2 = 132;// 2+37+40+53
		ResidueNumber start = new ResidueNumber("A", 246, null);
		ResidueNumber mid = new ResidueNumber("A", 85, 'S');
		ResidueNumber end = new ResidueNumber("A", 300, null);
		int realLength1 = map.getLength(start, mid);
		assertEquals("Real atom length is wrong", length1, realLength1);
		int realLength2 = map.getLength(start, end);
		assertEquals("Real atom length is wrong", length2, realLength2);

		int realLength = map.getLength(new ResidueNumber("A",6,'P'),new ResidueNumber("A",338,null));
		assertEquals("Full length wrong",430,realLength);
	}

	@Test
	public void testTrim() throws IOException, StructureException {
		// Two identical chains, residues 1-68, no insertion codes or missing residues
		String pdbId = "3w0e";
		AtomPositionMap map = new AtomPositionMap(cache.getAtoms(pdbId));

		ResidueRangeAndLength trimmed;
		ResidueRange untrimmed;

		untrimmed = new ResidueRange("A", "1", "68");
		trimmed = map.trimToValidResidues(untrimmed);
		assertEquals("Wrong start after trimming "+untrimmed,new ResidueNumber("A", 1, null), trimmed.getStart());
		assertEquals("Wrong end after trimming "+untrimmed,new ResidueNumber("A", 68, null), trimmed.getEnd());
		assertEquals("Wrong length after trimming "+untrimmed,68, trimmed.getLength());

		untrimmed = new ResidueRange("A", "1", "1");
		trimmed = map.trimToValidResidues(untrimmed);
		assertEquals("Wrong start after trimming "+untrimmed,new ResidueNumber("A", 1, null), trimmed.getStart());
		assertEquals("Wrong end after trimming "+untrimmed,new ResidueNumber("A", 1, null), trimmed.getEnd());
		assertEquals("Wrong length after trimming "+untrimmed,1, trimmed.getLength());

		untrimmed = new ResidueRange("A", "-1", "70");
		trimmed = map.trimToValidResidues(untrimmed);
		assertEquals("Wrong start after trimming "+untrimmed,new ResidueNumber("A", 1, null), trimmed.getStart());
		assertEquals("Wrong end after trimming "+untrimmed,new ResidueNumber("A", 68, null), trimmed.getEnd());
		assertEquals("Wrong length after trimming "+untrimmed,68, trimmed.getLength());

		// Fully out of range
		untrimmed = new ResidueRange("A", "-4", "-1");
		trimmed = map.trimToValidResidues(untrimmed);
		assertNull("Should be empty range "+untrimmed,trimmed);

		// Start before end. Arguably should be invalid, but currently works
		untrimmed = new ResidueRange("A", "4", "1");
		trimmed = map.trimToValidResidues(untrimmed);
		assertEquals("Wrong start after trimming "+untrimmed,new ResidueNumber("A", 4, null), trimmed.getStart());
		assertEquals("Wrong end after trimming "+untrimmed,new ResidueNumber("A", 1, null), trimmed.getEnd());
		assertEquals("Wrong length after trimming "+untrimmed,4, trimmed.getLength());

		// However, doesn't work if the endpoints are invalid, since searches wrong direction
		untrimmed = new ResidueRange("A", "4", "-1");
		trimmed = map.trimToValidResidues(untrimmed);
		assertNull("Should be empty range "+untrimmed,trimmed);

		untrimmed = new ResidueRange("A", "70", "50");
		trimmed = map.trimToValidResidues(untrimmed);
		assertNull("Should be empty range "+untrimmed,trimmed);

	}

	@Test
	public void testChains() throws IOException, StructureException {
		String pdbId = "1qdm";
		AtomPositionMap map = new AtomPositionMap(cache.getAtoms(pdbId));

		ResidueNumber start,end;
		try {
			start = new ResidueNumber("A",6,'P');
			end = new ResidueNumber("B",338,null);
			map.getLength(start, end);
			fail("Chain missmatch");
		} catch(IllegalArgumentException e) {
			// Expected
		}
		try {
			start = new ResidueNumber("A",6,'P');
			end = new ResidueNumber("B",338,null);
			map.getLengthDirectional(start, end);
			fail("Chain missmatch");
		} catch(IllegalArgumentException e) {
			// Expected
		}

		// With integers, only count matching chain atoms
		start = new ResidueNumber("A",338,null);
		end = new ResidueNumber("B",6,'P');
		int len = map.getLength(map.getPosition(start),map.getPosition(end),"A");
		assertEquals(1, len);
	}
}

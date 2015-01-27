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

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;import static org.junit.Assume.*;

import java.io.IOException;
import java.util.NavigableMap;

import static org.junit.Assert.assertEquals;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.LocalPDBDirectory;
import org.biojava.bio.structure.io.LocalPDBDirectory.FetchBehavior;
import org.junit.Before;
import org.junit.Test;

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
	
	@Test
	public void testLengths() throws IOException, StructureException {
		String pdbId = "4hhb";
		cache.setUseMmCif(false);
		cache.setFetchBehavior(FetchBehavior.LOCAL_ONLY);
		Atom[] atoms = cache.getAtoms(pdbId);
		AtomPositionMap map = new AtomPositionMap(cache.getAtoms(pdbId));

		int len;
		int start,end;
		
		// Single residue
		start = 1;
		end = 1;
		len = map.calcLength(new ResidueNumber("A",start,null), new ResidueNumber("A",end,null));
		assertEquals("Bad length for ("+start+","+end+")",end-start, len);
		len = map.calcLength(start, end,"A");
		assertEquals("Bad length for ("+start+","+end+")",end-start, len);
		len = map.calcLengthDirectional(start, end, "A");
		assertEquals("Bad length for ("+start+","+end+")",end-start, len);
		len = map.calcLengthDirectional(end, start, "A");
		assertEquals("Bad length for ("+start+","+end+")",0, len);

		// Short range
		start = 1;
		end = 3;
		len = map.calcLength(new ResidueNumber("A",start,null), new ResidueNumber("A",end,null));
		assertEquals("Bad length for ("+start+","+end+")",end-start, len);
		len = map.calcLength(start, end,"A");
		assertEquals("Bad length for ("+start+","+end+")",end-start, len);
		len = map.calcLengthDirectional(start, end, "A");
		assertEquals("Bad length for ("+start+","+end+")",end-start, len);
		len = map.calcLengthDirectional(end, start, "A");
		assertEquals("Bad length for ("+start+","+end+")",0, len);
		
		// Double check that the chain length is correct
		int chainAlen = cache.getStructure(pdbId).getChainByPDB("A").getAtomGroups(GroupType.AMINOACID).size();
		assumeTrue(141==chainAlen);

		//Full chain
		start = 1;
		end = chainAlen;
		len = map.calcLength(new ResidueNumber("A",start,null), new ResidueNumber("A",end,null));
		assertEquals("Bad length for ("+start+","+end+")",end-start, len);
		len = map.calcLength(start, end,"A");
		assertEquals("Bad length for ("+start+","+end+")",end-start, len);
		len = map.calcLengthDirectional(start, end, "A");
		assertEquals("Bad length for ("+start+","+end+")",end-start, len);
		len = map.calcLengthDirectional(end, start, "A");
		assertEquals("Bad length for ("+start+","+end+")",0, len);
		

		//Fuller chain
		try {
			start = 1;
			end = chainAlen+1;
			len = map.calcLength(new ResidueNumber("A",start,null), new ResidueNumber("A",end,null));
			assertEquals("Bad length for ("+start+","+end+")",end-start, len);
			len = map.calcLength(start, end,"A");
			assertEquals("Bad length for ("+start+","+end+")",end-start, len);
			len = map.calcLengthDirectional(start, end, "A");
			assertEquals("Bad length for ("+start+","+end+")",end-start, len);
			len = map.calcLengthDirectional(end, start, "A");
			assertEquals("Bad length for ("+start+","+end+")",0, len);
		} catch( NullPointerException e) {
			// WTF
		}

		
		// Chain spanning
		start = chainAlen;
		end = chainAlen+1;
		len = map.calcLength(new ResidueNumber("A",start,null), new ResidueNumber("B",end-chainAlen,null));
		assertEquals("Bad length for ("+start+","+end+")",-1/*end-start*/, len);
		len = map.calcLength(start, end,"A");
		assertEquals("Bad length for ("+start+","+end+")",-1/*end-start*/, len);
		len = map.calcLengthDirectional(start, end, "A");
		assertEquals("Bad length for ("+start+","+end+")",-1/*end-start*/, len);
		len = map.calcLengthDirectional(end, start, "A");
		assertEquals("Bad length for ("+start+","+end+")",chainAlen-1, len);

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

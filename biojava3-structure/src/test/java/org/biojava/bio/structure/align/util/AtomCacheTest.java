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
package org.biojava.bio.structure.align.util;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import org.biojava.bio.structure.AtomPositionMap;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.ResidueRange;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.scop.ScopFactory;
import org.junit.Before;
import org.junit.Test;


/**
 * A test for {@link AtomCache}.
 * @author dmyerstu
 * @since 3.0.6
 */
public class AtomCacheTest {

	private AtomCache cache;
	
	@Before
	public void setUp() {
		cache = new AtomCache();
		cache.setFetchFileEvenIfObsolete(true);
		cache.setStrictSCOP(false);
		// Use a fixed SCOP version for stability
		ScopFactory.setScopDatabase(ScopFactory.VERSION_1_75B);
	}
	
	/**
	 * Tests {@link AtomCache#getStructureForDomain(String)} on a multi-chain domain with no ligands but an explicit range (not whole-chain).
	 */
	@Test
	public void testGetStructureForDomain1() throws IOException, StructureException {
		String ranges = "A:328-396,B:518-527";
		Structure whole = cache.getStructure("1h6w");
		AtomPositionMap map = new AtomPositionMap(StructureTools.getAllAtomArray(whole), AtomPositionMap.ANYTHING_MATCHER);
		List<ResidueRange> rrs = ResidueRange.parseMultiple(ranges, map);
		int expectedLengthA = rrs.get(0).getLength();
		int expectedLengthB = rrs.get(1).getLength();
		Structure structure = cache.getStructureForDomain("d1h6w.2");
		assertEquals(2, structure.getChains().size());
		Chain a = structure.getChainByPDB("A");
		Chain b = structure.getChainByPDB("B");
		// we're one off because getAtomGroups().size() includes the last Group
		assertEquals(expectedLengthA, a.getAtomGroups().size()-1);
		assertEquals(expectedLengthB, b.getAtomGroups().size()-1);
	}

	/**
	 * Tests {@link AtomCache#getStructureForDomain(String)} on a multi-chain domain with two zinc ligands that occurs after the TER. The ligands are in chains E and F, so they should not be included in the domain.
	 */
	@Test
	public void testGetStructureForDomain2() throws IOException, StructureException {
		String ranges = "A:,B:";
		Structure whole = cache.getStructure("1I3O");
		AtomPositionMap map = new AtomPositionMap(StructureTools.getAllAtomArray(whole), AtomPositionMap.ANYTHING_MATCHER);
		List<ResidueRange> rrs = ResidueRange.parseMultiple(ranges, map);
		int expectedLengthA = rrs.get(0).getLength();
		int expectedLengthB = rrs.get(1).getLength();
		Structure structure = cache.getStructureForDomain("d1i3o.1");
		assertEquals(2, structure.getChains().size());
		Chain a = structure.getChainByPDB("A");
		Chain b = structure.getChainByPDB("B");
		// we're one off because getAtomGroups().size() includes the last Group
		assertEquals(expectedLengthA, a.getAtomGroups().size()-1);
		assertEquals(expectedLengthB, b.getAtomGroups().size()-1);
		List<Group> ligandsA = StructureTools.filterLigands(b.getAtomGroups());
		assertEquals(0, ligandsA.size());
		List<Group> ligandsB = StructureTools.filterLigands(b.getAtomGroups());
		assertEquals(0, ligandsB.size());
	}

	/**
	 * Tests {@link AtomCache#getStructureForDomain(String)} on a single-chain domain with two zinc ligands that occurs after the TER. 
	 */
	@Test
	public void testGetStructureForDomain3() throws IOException, StructureException {
		String ranges = "E:";
		Structure whole = cache.getStructure("1I3O");
		AtomPositionMap map = new AtomPositionMap(StructureTools.getAllAtomArray(whole), AtomPositionMap.ANYTHING_MATCHER);
		List<ResidueRange> rrs = ResidueRange.parseMultiple(ranges, map);
		int expectedLengthE = rrs.get(0).getLength();
		Structure structure = cache.getStructureForDomain("d1i3oe_");
		assertEquals(1, structure.getChains().size());
		Chain e = structure.getChainByPDB("E");
		// we're one off because getAtomGroups().size() includes the last Group
		assertEquals(expectedLengthE, e.getAtomGroups().size()-1);
		List<Group> ligandsE = StructureTools.filterLigands(e.getAtomGroups());
		assertEquals(1, ligandsE.size());
	}
	
}

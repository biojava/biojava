/**
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
 * Created on Nov 28, 2012
 * Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package org.biojava.nbio.structure;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertNotNull;

import java.io.IOException;
import java.util.Iterator;
import java.util.List;

import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.junit.Test;

public class TestCloning {

	@Test
	public void test1a4wCloning() throws StructureException, IOException {

		Structure s;

		final AtomCache cache = new AtomCache();
		final FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(true);
		cache.setFileParsingParams(params);

		StructureIO.setAtomCache(cache);

		s = StructureIO.getStructure("1a4w");

		final Structure c = s.clone();

		compareCloned(s, c);

	}

	@Test
	public void testAsymUnitCloning() throws StructureException, IOException {

		Structure s;

		final AtomCache cache = new AtomCache();
		final FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(false);
		cache.setFileParsingParams(params);

		StructureIO.setAtomCache(cache);

		s = StructureIO.getStructure("1stp");

		final Structure c = s.clone();

		compareCloned(s, c);
	}

	@Test
	public void testBioUnitCloning() throws StructureException, IOException {

		Structure s;
		s = StructureIO.getBiologicalAssembly("1stp", 1);

		final Structure c = s.clone();

		compareCloned(s, c);

	}

	/**
	 * A Structure with alt locs, we make sure they are being cloned too
	 *
	 * @throws StructureException
	 * @throws IOException
	 */
	@Test
	public void test3piuCloning() throws StructureException, IOException {

		final AtomCache cache = new AtomCache();
		final FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(true);
		cache.setFileParsingParams(params);

		StructureIO.setAtomCache(cache);

		final Structure s = StructureIO.getStructure("3piu");

		final Structure c = s.clone();

		compareCloned(s, c);
	}

	private void compareCloned(final Structure s, final Structure c) throws StructureException {

		assertEquals(s.getChains().size(), c.getChains().size());

		for (final Chain chain : s.getChains()) {

			final Chain test = c.getChain(chain.getId());

			assertEquals("Could not correctly clone seqres for chain " + chain.getId(), chain.getSeqResLength(),
					test.getSeqResLength());

			assertEquals("Could not correctly clone atom records for chain " + chain.getId(), chain.getAtomLength(),
					test.getAtomLength());

			Iterator<Group> it = test.getAtomGroups().iterator();
			for (final Group g : chain.getAtomGroups()) {
				final Group testGroup = it.next();
				// if (g.hasAltLoc()) {
				// System.out.println(g.toString());
				// }
				assertEquals(g.getAltLocs().size(), testGroup.getAltLocs().size());
			}

			it = test.getSeqResGroups().iterator();
			for (final Group g : chain.getSeqResGroups()) {
				final Group testGroup = it.next();
				assertEquals(g.getAltLocs().size(), testGroup.getAltLocs().size());
			}
		}

		final Atom[] allAtoms = StructureTools.getAllAtomArray(s);

		final Atom[] allAtomsCloned = StructureTools.getAllAtomArray(c);

		assertEquals(allAtoms.length, allAtomsCloned.length);

	}

	@Test
	public void testBondCloning() throws IOException, StructureException {

		final AtomCache cache = new AtomCache();
		cache.setUseMmCif(true);

		final FileParsingParameters params = cache.getFileParsingParams();
		params.setCreateAtomBonds(true);
		cache.setFileParsingParams(params);

		final Structure s = cache.getStructure("2I13");
		List<Bond> bonds = s.getNonPolyChain("G").getAtomGroup(0).getAtom(0).getBonds();
		assertNotNull(bonds);

		Structure s2 = s.clone();
		List<Bond> bonds2 = s2.getNonPolyChain("G").getAtomGroup(0).getAtom(0).getBonds();
		assertNotNull(bonds2);

		assertEquals(bonds.toString(), bonds2.toString());
		// But the objects should be different as the atoms are clones
		assertNotEquals(bonds.toArray(), bonds2.toArray());

		// Also test for polymeric chains
		bonds = s.getPolyChain("E").getAtomGroup(0).getAtom(0).getBonds();
		assertNotNull(bonds);

		s2 = s.clone();
		bonds2 = s2.getPolyChain("E").getAtomGroup(0).getAtom(0).getBonds();
		assertNotNull(bonds2);
	}

}

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
package org.biojava.nbio.structure;

import org.biojava.nbio.structure.align.util.AtomCache;
import org.junit.Test;

import java.io.IOException;
import java.util.Iterator;

import static org.junit.Assert.*;

public class TestExperimentalTechniques {

	@Test
	public void test6F2Q() throws IOException, StructureException {

		// a multiple experimental techniques PDB entry (X-RAY + NEUTRON DIFFRACTION)

		AtomCache cache = new AtomCache();

		StructureIO.setAtomCache(cache);

		cache.setUseMmCif(false);
		Structure sPdb = StructureIO.getStructure("6F2Q");
		cache.setUseMmCif(true);
		Structure sCif = StructureIO.getStructure("6F2Q");

		comparePdbToCif(sPdb, sCif);

		assertEquals(2, sPdb.getPDBHeader().getExperimentalTechniques().size());
		assertEquals(2, sCif.getPDBHeader().getExperimentalTechniques().size());



		assertTrue(sPdb.isCrystallographic());
		assertTrue(sCif.isCrystallographic());

		assertFalse(sPdb.isNmr());
		assertFalse(sCif.isNmr());


	}

	@Test
	public void test3ZPK() throws IOException, StructureException {

		// a multiple experimental techniques PDB entry (EM + SOLUTION NMR)

		AtomCache cache = new AtomCache();

		StructureIO.setAtomCache(cache);

		cache.setUseMmCif(false);
		Structure sPdb = StructureIO.getStructure("3ZPK");
		cache.setUseMmCif(true);
		Structure sCif = StructureIO.getStructure("3ZPK");

		comparePdbToCif(sPdb, sCif);

		assertEquals(2, sPdb.getPDBHeader().getExperimentalTechniques().size());
		assertEquals(2, sCif.getPDBHeader().getExperimentalTechniques().size());



		assertFalse(sPdb.isCrystallographic());
		assertFalse(sCif.isCrystallographic());

		assertTrue(sPdb.isNmr());
		assertTrue(sCif.isNmr());


	}

	@Test
	public void test2B6O() throws IOException, StructureException {

		// a single experimental technique ELECTRON CRYSTALLOGRAPHY entry

		AtomCache cache = new AtomCache();

		StructureIO.setAtomCache(cache);

		cache.setUseMmCif(false);
		Structure sPdb = StructureIO.getStructure("2B6O");
		cache.setUseMmCif(true);
		Structure sCif = StructureIO.getStructure("2B6O");

		comparePdbToCif(sPdb, sCif);

		assertEquals(1, sPdb.getPDBHeader().getExperimentalTechniques().size());
		assertEquals(1, sCif.getPDBHeader().getExperimentalTechniques().size());



		assertTrue(sPdb.isCrystallographic());
		assertTrue(sCif.isCrystallographic());

		assertFalse(sPdb.isNmr());
		assertFalse(sCif.isNmr());


	}

	@Test
	public void test4CSO() throws IOException, StructureException {

		// a single experimental technique (X-RAY) entry

		AtomCache cache = new AtomCache();

		StructureIO.setAtomCache(cache);

		cache.setUseMmCif(false);
		Structure sPdb = StructureIO.getStructure("4CSO");
		cache.setUseMmCif(true);
		Structure sCif = StructureIO.getStructure("4CSO");

		comparePdbToCif(sPdb, sCif);

		assertEquals(1, sPdb.getPDBHeader().getExperimentalTechniques().size());
		assertEquals(1, sCif.getPDBHeader().getExperimentalTechniques().size());



		assertTrue(sPdb.isCrystallographic());
		assertTrue(sCif.isCrystallographic());

		assertFalse(sPdb.isNmr());
		assertFalse(sCif.isNmr());


	}

	private void comparePdbToCif(Structure sPdb, Structure sCif) {
		assertNotNull(sPdb.getPDBHeader().getExperimentalTechniques());
		assertNotNull(sCif.getPDBHeader().getExperimentalTechniques());

		Iterator<ExperimentalTechnique> itCif = sCif.getPDBHeader().getExperimentalTechniques().iterator();
		for (ExperimentalTechnique et:sPdb.getPDBHeader().getExperimentalTechniques()) {
			assertEquals(et, itCif.next());
		}
	}

}

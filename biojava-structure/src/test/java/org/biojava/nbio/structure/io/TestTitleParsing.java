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
package org.biojava.nbio.structure.io;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.junit.Test;

import java.io.IOException;

import static org.junit.Assert.*;

/**
 * Testing for title parsing in PDB and mmCIF files
 * especially multi-line titles in PDB
 *
 * @author duarte_j
 *
 */
public class TestTitleParsing {

	@Test
	public void test2W6E() throws IOException, StructureException {

		// an entry with a title in multiple lines in PDB file

		AtomCache cache = new AtomCache();

		StructureIO.setAtomCache(cache);

		cache.setUseMmCif(false);
		Structure sPdb = StructureIO.getStructure("2W6E");
		cache.setUseMmCif(true);
		Structure sCif = StructureIO.getStructure("2W6E");


		// we can only compare titles by first forcing lower case, since the cases don't coincide cif vs pdb
		assertEquals(sPdb.getPDBHeader().getTitle().toLowerCase().replace(" ", ""), sCif.getPDBHeader().getTitle().toLowerCase().replace(" ", ""));

		assertEquals(1, sPdb.getPDBHeader().getExperimentalTechniques().size());
		assertEquals(1, sCif.getPDBHeader().getExperimentalTechniques().size());

		assertEquals(sPdb.getPDBHeader().getResolution(),sCif.getPDBHeader().getResolution(),0.001);

		// an x-ray entry
		assertTrue(sPdb.isCrystallographic());
		assertTrue(sCif.isCrystallographic());

		assertFalse(sPdb.isNmr());
		assertFalse(sCif.isNmr());


	}



}

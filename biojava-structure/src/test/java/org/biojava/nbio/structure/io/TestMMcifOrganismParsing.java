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
 * created at Aug 12, 2013
 * Author: ap3
 */

package org.biojava.nbio.structure.io;


import org.biojava.nbio.structure.Compound;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.IOException;

import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;


public class TestMMcifOrganismParsing {



	@BeforeClass
	public static void setUp() throws Exception {

		AtomCache cache = new AtomCache();
		cache.setUseMmCif(true);
		StructureIO.setAtomCache(cache);
	}

	@Test
	public void test1STP() throws IOException, StructureException{
		String pdbId = "1stp";

		checkPDB(pdbId);

	}

	// removed this test, since entity 3 of 1a4w has no organism tax_id
	public void test1a4w() throws IOException, StructureException{
		String pdbId = "1a4w";

		checkPDB(pdbId);

	}

	@Test
	public void test4hhb() throws IOException, StructureException{
		String pdbId = "4hhb";

		checkPDB(pdbId);

	}

	@Test
	public void test3ZD6() throws IOException, StructureException {
		// a PDB ID that contains a synthetic entity
		String pdbId = "3ZD6";

		checkPDB(pdbId);

	}




	private void checkPDB(String pdbId) throws IOException, StructureException {
		Structure s = StructureIO.getStructure(pdbId);

		assertNotNull(s.getCompounds());
		assertTrue(s.getCompounds().size() > 0);

		for ( Compound c : s.getCompounds()) {
			assertNotNull(c.getOrganismTaxId());
		}

	}

}

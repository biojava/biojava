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


import org.biojava.nbio.structure.EntityInfo;
import org.biojava.nbio.structure.EntityType;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.IOException;

import static org.junit.Assert.assertEquals;
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

		checkPDB(pdbId, "1895");

	}

	// removed this test, since entity 3 of 1a4w has no organism tax_id
	public void test1a4w() throws IOException, StructureException{
		String pdbId = "1a4w";

		checkPDB(pdbId, "9606");

	}

	@Test
	public void test4hhb() throws IOException, StructureException{
		String pdbId = "4hhb";

		checkPDB(pdbId, "9606");

	}

	@Test
	public void test3ZD6() throws IOException, StructureException {
		// a PDB ID that contains a synthetic entity
		String pdbId = "3zd6";

		checkPDB(pdbId, "9606");

	}


	private void checkPDB(String pdbId, String organismTaxId) throws IOException, StructureException {
		Structure s = StructureIO.getStructure(pdbId);

		assertNotNull(s.getEntityInfos());
		assertTrue(s.getEntityInfos().size() > 0);

		for ( EntityInfo c : s.getEntityInfos()) {
			if(EntityType.POLYMER.equals(c.getType())) { 
				assertNotNull(c.getOrganismTaxId());
				if(pdbId.equals("3zd6")){
					if(c.getMolId()==2) {
						assertEquals(c.getOrganismTaxId(), "32630");
						continue;
					}
				}
				assertEquals(c.getOrganismTaxId(), organismTaxId);
			
			}
		}

	}

}

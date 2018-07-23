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
package org.biojava.nbio.structure.secstruc;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.zip.GZIPInputStream;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 * Test the correctness of the DSSP implementation in BioJava
 * for the calculation of secondary structure in a Structure object.
 *
 * EXAMPLES:
 * 			Big structures: 4v7r, 4V60 (use mmCif parser)
 * 			Helical: 4hhb, 4lup
 * 			Mixed small: 5pti
 * 			First sheet: 1ze3, 3k19
 * 			Insertion code: 1how
 *          More than 2 Beta-Bridges: 2k4t
 *
 * @author Aleix Lafita
 *
 */
public class TestSecStrucCalc {

	@Test
	public void testSecStrucPred() throws StructureException, IOException {

		//List of names to test the DSSP prediction
		List<String> names = Arrays.asList(
				"5pti", "1tim", "4hhb", "1how", "4i4q", "2k4t", "1deu");
		SecStrucCalc sec = new SecStrucCalc();
		//Predict with BioJava the SS -> Anthony has moved this out of the loop.
		//SecStrucCalc does not need to be reinitialised every time
		for (String name : names) {

			AtomCache cache = new AtomCache();
			Structure structure = cache.getStructure(name);


			List<SecStrucState> biojava = sec.calculate(structure, true);

			//Download the original DSSP implementation output
			List<SecStrucState> dssp = DSSPParser.parseInputStream(new GZIPInputStream(
					this.getClass().getResourceAsStream("/org/biojava/nbio/structure/secstruc/"+name+".dssp.gz")), structure, false);

			assertEquals("SS assignment lengths do not match",
					biojava.size(), dssp.size()*structure.nrModels());

			for (int i=0; i<dssp.size(); i++){
				assertEquals("SS assignment position "+(i+1)+" does not match",
						biojava.get(i), dssp.get(i));
			}
		}
	}

	
	/**
	 * Test that calculating the secondary structure for multi-model systems works.
	 * Combine two PDBs into one multi-model system
	 * Calculate the secondary structure
	 * Combine with the combined list fetching from the server
	 * @throws StructureException
	 * @throws IOException
	 */
	@Test
	public void testMultiModelPred() throws StructureException, IOException {

		String pdbId = "5pti";
		String pdbIdTwo = "4hhb";
		SecStrucCalc sec = new SecStrucCalc();
		// Combine these into one structure with two models
		AtomCache cache = new AtomCache();
		Structure structure = cache.getStructure(pdbId);
		Structure structureTwo = cache.getStructure(pdbIdTwo);
		// Join them together
		structure.addModel(structureTwo.getChains());
		
		List<SecStrucState> biojava = sec.calculate(structure, true);

		// Download the original DSSP implementation output
		List<SecStrucState> dssp = DSSPParser.parseInputStream(new GZIPInputStream(
				this.getClass().getResourceAsStream("/org/biojava/nbio/structure/secstruc/"+pdbId+".dssp.gz")),cache.getStructure(pdbId), false);
		dssp.addAll(DSSPParser.parseInputStream(new GZIPInputStream(
				this.getClass().getResourceAsStream("/org/biojava/nbio/structure/secstruc/"+pdbIdTwo+".dssp.gz")), cache.getStructure(pdbIdTwo), false));
		
		assertEquals("SS assignment lengths do not match",
				biojava.size(), dssp.size());

		for (int i=0; i<dssp.size(); i++){
			assertEquals("SS assignment position "+(i+1)+" does not match",
					biojava.get(i), dssp.get(i));
		}
	}
}

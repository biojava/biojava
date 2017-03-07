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
package org.biojava.nbio.structure.cluster;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.List;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.junit.Test;

/**
 * Test the {@link SubunitExtractor} correctness on different real structures
 * with different types of difficulties.
 * 
 * @author Aleix Lafita
 *
 */
public class TestSubunitExtractor {

	/**
	 * Some collagen structures have very short Chains, so the minimum sequence
	 * length is adjusted: 1A3I.
	 * 
	 * @see SubunitClustererParameters#getMinimumSequenceLengthFraction()
	 */
	@Test
	public void testCollagen() throws StructureException, IOException {

		Structure s = StructureIO.getStructure("1A3I");

		List<Subunit> subunits = SubunitExtractor.extractSubunits(s, 5, 0.75, 20);

		// We expect all 3 subunits to be returned
		assertEquals(subunits.size(), 3);

		subunits = SubunitExtractor.extractSubunits(s, 8, 0.75, 9);

		// Now we expect only the long Subunit to be returned
		assertEquals(subunits.size(), 1);
		assertEquals(subunits.get(0).size(), 9);
	}

	/**
	 * Make sure that only aminoacid chains are extracted: 5B2I.
	 */
	@Test
	public void testHistone() throws StructureException, IOException {

		Structure s = StructureIO.getStructure("5B2I");

		List<Subunit> subunits = SubunitExtractor.extractSubunits(s, 5, 0.75, 20);

		// We expect all 8 histone subunits to be returned
		assertEquals(subunits.size(), 8);
		assertEquals(subunits.get(0).size(), 99);
		assertEquals(subunits.get(1).size(), 82);
		assertEquals(subunits.get(2).size(), 106);
	}

	/**
	 * Test that all chains from biological assemblies are extracted.
	 */
	@Test
	public void testBioAssembly() throws StructureException, IOException {

		Structure s = StructureIO.getStructure("BIO:4E3E:1");

		List<Subunit> subunits = SubunitExtractor.extractSubunits(s, 5, 0.75, 20);

		// We expect all 3 equal double hot-dog subunits to be returned
		assertEquals(subunits.size(), 3);
		assertEquals(subunits.get(0).size(), subunits.get(1).size());
		assertEquals(subunits.get(0).size(), subunits.get(2).size());
	}
}

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
package org.biojava.nbio.structure.test.io;

import static org.junit.Assert.*;

import java.io.IOException;

import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.junit.Test;

/**
 * Test StructureIO methods.
 * 
 */
public class StructureIOTest {
	
	@Test
	public void testStructureIO() throws IOException, StructureException {

		String pdbId = "1gav";
		int nrAssembls = StructureIO.getBiologicalAssemblies(pdbId).size();
		assertEquals(1,nrAssembls);

		pdbId = "1hv4";
		nrAssembls = StructureIO.getBiologicalAssemblies(pdbId).size();
		assertEquals(2,nrAssembls);

	}
}

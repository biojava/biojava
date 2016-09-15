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

import static org.junit.Assert.*;

import java.io.IOException;

import org.biojava.nbio.structure.align.util.AtomCache;
import org.junit.Test;

public class TestBioAssemblyIdentifier {


	@Test
	public void test() throws IOException, StructureException {

		AtomCache cache = new AtomCache();
		BioAssemblyIdentifier id;
		Structure s;

		// first assembly
		id = new BioAssemblyIdentifier("BIO:2ehz:1");
		s = cache.getStructure(id);
		assertEquals("Number of models",1, s.nrModels());
		assertEquals("Number of chains",88, s.getChains().size());
		// equivalent
		id = new BioAssemblyIdentifier("BIO:2ehz");
		s = cache.getStructure(id);
		assertEquals("Number of models",1, s.nrModels());
		assertEquals("Number of chains",8,s.getPolyChains().size());
		// No second
		id = new BioAssemblyIdentifier("BIO:2ehz:2");
		try {
			s = cache.getStructure(id);
			fail("Expected exception for invalid assembly number");
		} catch( StructureException e) {}
		// AU
		id = new BioAssemblyIdentifier("BIO:2ehz:0");
		s = cache.getStructure(id);
		assertEquals("Number of models",1, s.nrModels());
		assertEquals("Number of chains per model",1,s.getPolyChains(0).size());

	}

}

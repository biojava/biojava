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

import org.junit.Test;
import static org.junit.Assert.*;


import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;

public class TestStructWithMultiparentChemComp {
	
	@Test
	public void test4Q7U() throws Exception {
		
		AtomCache cache = new AtomCache();
		cache.setUseMmCif(true);
		StructureIO.setAtomCache(cache);
		
 		Structure s = StructureIO.getStructure("4q7u");
 		assertNotNull(s);
 		Chain c = s.getPolyChain("A");
 		
 		String seq = c.getSeqResSequence();
 		assertNotNull(seq);
 		assertEquals(245, seq.length());
 		assertFalse(seq.contains("?"));
 		assertTrue(seq.contains("X"));
		
	}

}

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
 * Created on Apr 20, 2012
 * Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package org.biojava.bio.structure;

import java.io.IOException;
import java.util.List;

import junit.framework.TestCase;

import org.biojava.bio.structure.align.util.AtomCache;

/** This class tests the correct loading of Nucleotides
 * 
 * @author Andreas Prlic
 * @since 3.0.3
 */
public class TestNucleotides extends TestCase{

	AtomCache cache = new AtomCache();

	public void test3T5N(){

		Structure s = null;
		try {
			s = getStructure("3T5N");

		} catch (Exception e) {
			e.printStackTrace();
			fail(e.getMessage());

		}

		assertEquals(2,s.getChains().size());

		Chain c = s.getChains().get(1);
		assertEquals("C", c.getChainID());
		List<Group> ngr = c.getAtomGroups("nucleotide");
		assertEquals(6,ngr.size());
	}


	public void test1OFX(){
		Structure s = null;
		try {
			s = getStructure("1OFX");

		} catch (Exception e) {
			e.printStackTrace();
			fail(e.getMessage());

		}

		assertEquals(2,s.getChains().size());

		Chain a = s.getChains().get(0);
		assertEquals("A", a.getChainID());
		List<Group> ngrA = a.getAtomGroups("nucleotide");
		assertEquals(10,ngrA.size());
		
		Chain b = s.getChains().get(1);
		assertEquals("B", b.getChainID());
		List<Group> ngrB = b.getAtomGroups("nucleotide");
		assertEquals(10,ngrB.size());
	}

private Structure getStructure(String pdbId) throws IOException, StructureException {
	return cache.getStructure(pdbId);
}
}

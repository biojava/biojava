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

import java.io.IOException;
import java.util.List;
import static org.junit.Assert.*;

import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.junit.Test;

public class TestParseOnAsymId {



	@Test
	public void test4cup() throws IOException, StructureException {

		AtomCache cache = new AtomCache();
		cache.setUseMmCif(true);
		FileParsingParameters params = cache.getFileParsingParams();

		cache.setFileParsingParams(params);
		StructureIO.setAtomCache(cache);
		// Set the test lists
		String[] asymChainList = {"A","B","C","D","E","F"};
		String[] authChainList = {"A","A","A","A","A","A"};
		String[] asymChainListTest = new String[6];
		String[] authChainListTest = new String[6];
		// Get the structure
		Structure bioJavaStruct = StructureIO.getStructure("4cup");
		List<Chain> chainList = bioJavaStruct.getChains();
		assertEquals(6,chainList.size());
		for(int i=0; i<chainList.size();i++){
			Chain c = chainList.get(i);
			authChainListTest[i] = c.getName();
			asymChainListTest[i] = c.getId();
		}
		// Now check both lists are the same
		assertArrayEquals(authChainListTest, authChainList);
		assertArrayEquals(asymChainListTest, asymChainList);
	}

}

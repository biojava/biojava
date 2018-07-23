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
package org.biojava.nbio.structure.io.mmcif;


import java.io.IOException;

import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 * Created by andreas on 5/3/16.
 */
public class TestParseInternalChainId {

	@Test
	public void test2I13() throws IOException, StructureException {


		AtomCache cache = new AtomCache();
		cache.setUseMmCif(true);


		Structure s = cache.getStructure("2I13");

		assertEquals(6, s.getPolyChains().size());
		assertEquals(15, s.getNonPolyChains().size());
		assertEquals(6, s.getWaterChains().size());

		/** the four nucleic chains
		 *
		 */
		Chain asymA = s.getPolyChain("A");
		Chain asymB = s.getPolyChain("B");
		Chain asymC = s.getPolyChain("C");
		Chain asymD = s.getPolyChain("D");

		/** the protein chains
		 *
		 */
		Chain asymE = s.getPolyChain("E");
		Chain asymF = s.getPolyChain("F");

		Chain[] nucleicChains = new Chain[]{asymA,asymB,asymC,asymD};
		Chain[] proteinChains = new Chain[]{asymE,asymF};

		for ( Chain c : proteinChains){
			assertNotNull("Chain is null!",c);
		}

		for ( Chain c : nucleicChains){
			assertNotNull("Chain is null!", c);
		}


		assertEquals("A", asymA.getId());
		assertEquals("C", asymA.getName());

		assertEquals("B", asymB.getId());
		assertEquals("D", asymB.getName());

		assertEquals("C", asymC.getId());
		assertEquals("E", asymC.getName());

		assertEquals("D", asymD.getId());
		assertEquals("F", asymD.getName());

		assertEquals("E", asymE.getId());
		assertEquals("A", asymE.getName());

		assertEquals("F", asymF.getId());
		assertEquals("B", asymF.getName());

		Chain chainG = s.getNonPolyChain("G");
		// before mmcif v5, this used to be "A", now it's "C" - JD 2017-07-13
		assertEquals("C", chainG.getName());

	}


}

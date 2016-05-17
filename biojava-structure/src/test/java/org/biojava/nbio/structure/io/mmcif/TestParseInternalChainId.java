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

		System.out.println(s);


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
			System.out.println(c);
			assertNotNull("Chain is null!",c);
		}

		for ( Chain c : nucleicChains){
			System.out.println(c);
			assertNotNull("Chain is null!", c);
		}


		assertTrue(asymA.getId().equals("A"));
		assertTrue(asymA.getName().equals("C"));

		assertTrue(asymB.getId().equals("B"));
		assertTrue(asymB.getName().equals("D"));

		assertTrue(asymC.getId().equals("C"));
		assertTrue(asymC.getName().equals("E"));

		assertTrue(asymD.getId().equals("D"));
		assertTrue(asymD.getName().equals("F"));

		assertTrue(asymE.getId().equals("E"));
		assertTrue(asymE.getName().equals("A"));

		assertTrue(asymF.getId().equals("F"));
		assertTrue(asymF.getName().equals("B"));

		Chain chainG = s.getNonPolyChain("G");
		assertTrue(chainG.getName().equals("A"));

	}


}

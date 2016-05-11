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
import java.util.List;

import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.junit.Test;

/**
 * Created by andreas on 9/16/15.
 */
public class Test2JA5 {

	@Test
	public void test2JA5() throws IOException, StructureException {

		FileParsingParameters fileParsingParameters = new FileParsingParameters();
		fileParsingParameters.setAlignSeqRes(true); // Need to align seqres to match chains.
		fileParsingParameters.setHeaderOnly(false); // Need header only off to have chains to match.

		AtomCache cache = new AtomCache();
		cache.setUseMmCif(false);
		cache.setFileParsingParams(fileParsingParameters);

		StructureIO.setAtomCache(cache);

		Structure s1 = StructureIO.getStructure("2ja5");

		// This is not applicable anymore, we need to parse atoms to have chains to match.
		// assertTrue(StructureTools.getNrAtoms(s1) == 0);

		// SeqRes contains 15 chains, but since we cannot align Chain N to AtomGroups => 14.
		assertEquals(14, s1.getPolyChains().size());

		Chain nChain = s1.getPolyChain("N");
		
		assertNotNull(nChain);

		List<Chain> chains = s1.getPolyChainsByPDB("N");
		assertTrue(chains.size() == 0);
	}

	@Test
	public void test2JA5noHeader() throws IOException, StructureException {

		FileParsingParameters fileParsingParameters = new FileParsingParameters();
		fileParsingParameters.setHeaderOnly(true);

		AtomCache cache = new AtomCache();
		cache.setUseMmCif(false);
		cache.setFileParsingParams(fileParsingParameters);

		StructureIO.setAtomCache(cache);


		Structure s1 = StructureIO.getStructure("2ja5");

		// This is not applicable anymore, we need to parse atoms to have chains to match.
		assertEquals(0, StructureTools.getNrAtoms(s1));

		// All 15 seqres chains will be store.
		assertEquals(15, s1.getPolyChains().size());

		Chain nChain = s1.getPolyChainsByPDB("N").get(0);
		
		assertNotNull(nChain);
	}
}

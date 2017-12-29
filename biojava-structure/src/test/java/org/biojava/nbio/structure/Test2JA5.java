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

		// SeqRes contains 14 chains, but since we cannot align Chain N to AtomGroups => 14.
		// 2ja5 has been remediated on March 2017, now it has 14 chains in seqres matching the 14 chains in atoms (chain N has been removed)
		assertEquals(14, s1.getPolyChains().size());

		Chain nChain = s1.getPolyChain("N");
		
		assertNotNull(nChain);

		Chain chain = s1.getPolyChainByPDB("N");
		assertNull(chain);
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

		// 2ja5 has been remediated on March 2017, now it has 14 chains in seqres matching the 14 chains in atoms (chain N has been removed)
		assertEquals(14, s1.getPolyChains().size());

		Chain nChain = s1.getPolyChainByPDB("N");
		
		assertNull(nChain);
	}
}

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

import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

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
	fileParsingParameters.setLoadChemCompInfo(true);
	fileParsingParameters.setHeaderOnly(false); // Need header only off to have chains to match.

	AtomCache cache = new AtomCache();
	cache.setUseMmCif(false);
	cache.setFileParsingParams(fileParsingParameters);

	StructureIO.setAtomCache(cache);


	Structure s1 = StructureIO.getStructure("2ja5");

	// This is not applicable anymore, we need to parse atoms to have chains to match.
	// assertTrue(StructureTools.getNrAtoms(s1) == 0);

	// SeqRes contains 15 chains, but since we cannot align Chain N to AtomGroups => 14.
	assertTrue(s1.getChains().size() == 14);

	Chain nChain = null;
	try {
		nChain = s1.getChainByPDB("N");
	} catch (StructureException e){
		// this is expected here, since there is no chain N
	}
	assertNull(nChain);
	}

	@Test
	public void test2JA5noHeader() throws IOException, StructureException {

	FileParsingParameters fileParsingParameters = new FileParsingParameters();
	fileParsingParameters.setLoadChemCompInfo(true);
	fileParsingParameters.setHeaderOnly(true);

	AtomCache cache = new AtomCache();
	cache.setUseMmCif(false);
	cache.setFileParsingParams(fileParsingParameters);

	StructureIO.setAtomCache(cache);


	Structure s1 = StructureIO.getStructure("2ja5");

	// This is not applicable anymore, we need to parse atoms to have chains to match.
	assertTrue(StructureTools.getNrAtoms(s1) == 0);

	// All 15 seqres chains will be store.
	assertTrue(s1.getChains().size() == 15);

	Chain nChain = null;
	try {
		nChain = s1.getChainByPDB("N");
	} catch (StructureException e){
		// this is expected here, since there is no chain N
	}
	assertNotNull(nChain);
	}
}

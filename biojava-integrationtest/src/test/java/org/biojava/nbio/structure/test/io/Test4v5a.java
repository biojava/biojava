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

//import static org.junit.Assert.*;

import java.io.IOException;


import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.junit.Ignore;
import org.junit.Test;

/**
 * Test for a difficult large mmCIF file (a ribosome) with many 2-letter chain ids.
 * Both tests are set to Ignore because the parsing of the file takes too long.
 *
 * @author Jose Duarte
 *
 */
public class Test4v5a {

	@Test @Ignore
	public void test4v5a() throws StructureException, IOException {
		AtomCache cache = new AtomCache();

		FileParsingParameters params = cache.getFileParsingParams();
		params.setCreateAtomBonds(true);
		StructureIO.setAtomCache(cache);

		StructureIO.getStructure("4v5a");


	}

	@Test @Ignore
	public void test4v5aInternalChainIds() throws StructureException, IOException {
		AtomCache cache = new AtomCache();

		FileParsingParameters params = cache.getFileParsingParams();
		params.setCreateAtomBonds(true);
		StructureIO.setAtomCache(cache);


		// this would throw an exception when making ss bonds from struct_conn because
		// chain id "AD" is both an author chain id and an asym chain id and the chains were
		// getting mixed
		StructureIO.getStructure("4v5a");

		//for (Chain c : s.getChains()) {
		//	System.out.println(c.getChainID()+" "+c.getInternalChainID()+" "+c.getAtomLength());
		//}


	}
}

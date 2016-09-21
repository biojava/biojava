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

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 * Created by ap3 on 31/07/2015.
 */
public class TestURLBasedFileParsing {

	@Test
	public void testMMcifURL() throws StructureException, IOException{

		String u = "http://ftp.wwpdb.org/pub/pdb/data/biounit/mmCIF/divided/nw/4nwr-assembly1.cif.gz";

		Structure s = StructureIO.getStructure(u);

		assertNotNull(s);


	}
}

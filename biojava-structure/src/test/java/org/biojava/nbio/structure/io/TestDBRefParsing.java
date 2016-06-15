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

import static org.junit.Assert.assertEquals;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.junit.Test;

public class TestDBRefParsing {

	@Test
	public void testShortLine() throws IOException, StructureException {
		Structure s,ref;
		PDBFileParser pdbPars = new PDBFileParser();

		String shortLine = "DBREF  2W6E A  -42   510  UNP    P19483   ATPA1_BOVIN      1    553";
		// Parse short
		try(InputStream is = new ByteArrayInputStream(shortLine.getBytes())) {
			s = pdbPars.parsePDBFile(is);
		}
		// Parse padded
		String longline = String.format("%1$-80s", shortLine);
		try(InputStream is = new ByteArrayInputStream(longline.getBytes()) ){
			ref = pdbPars.parsePDBFile(is);
		}
		assertEquals(ref.getDBRefs().get(0), s.getDBRefs().get(0));
	}

	@Test
	public void testToPdbLength() throws IOException {
		Structure s;
		String shortLine = "DBREF  2W6E A  -42   510  UNP    P19483   ATPA1_BOVIN      1    553";
		PDBFileParser pdbPars = new PDBFileParser();
		// Parse short
		try(InputStream is = new ByteArrayInputStream(shortLine.getBytes()) ) {
			s = pdbPars.parsePDBFile(is);
		}
		//Make sure that the record is true to the input
		assertEquals(shortLine, s.getDBRefs().get(0).toPDB().trim());
		//And is padded to the correct lenght
		assertEquals(80, s.getDBRefs().get(0).toPDB().length());
	}
}

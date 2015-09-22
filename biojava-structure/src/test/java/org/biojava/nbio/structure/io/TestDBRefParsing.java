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
		String shortLine = "DBREF  2W6E A  -42   510  UNP    P19483   ATPA1_BOVIN      1    553";
		InputStream is = new ByteArrayInputStream(shortLine.getBytes());
		PDBFileParser pdbPars = new PDBFileParser();
		Structure s;
		try {
			s = pdbPars.parsePDBFile(is);
		} catch (Exception e) {
			is.close();
			throw new AssertionError("Unable to process truncated DBRef line");
		}
		is = new ByteArrayInputStream(String.format("%1$-80s", shortLine)
				.getBytes());
		Structure ref = pdbPars.parsePDBFile(is);
		is.close();
		assertEquals(ref.getDBRefs().get(0), s.getDBRefs().get(0));
	}

	@Test
	public void testToPdbLength() throws IOException {
		String shortLine = "DBREF  2W6E A  -42   510  UNP    P19483   ATPA1_BOVIN      1    553";
		InputStream is = new ByteArrayInputStream(shortLine.getBytes());
		PDBFileParser pdbPars = new PDBFileParser();
		Structure s;
		try {
			s = pdbPars.parsePDBFile(is);
		} catch (Exception e) {
			is.close();
			throw new AssertionError("Unable to process truncated DBRef line");
		}
		//Make sure that the record is true to the input
		assertEquals(shortLine, s.getDBRefs().get(0).toPDB().trim());
		//And is padded to the correct lenght
		assertEquals(80, s.getDBRefs().get(0).toPDB().length());
	}
}

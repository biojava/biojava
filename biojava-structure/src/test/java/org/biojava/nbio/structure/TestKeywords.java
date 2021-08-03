package org.biojava.nbio.structure;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.io.InputStream;
import java.util.List;
import java.util.zip.GZIPInputStream;

import org.biojava.nbio.structure.io.PDBFileParser;
import org.junit.Test;

public class TestKeywords {

	@Test
	public void testKeywordsOnFiveLines () throws IOException {
		String fileName = "/3cdl.pdb";
		InputStream inStream = this.getClass().getResourceAsStream(fileName);

		PDBFileParser pdbpars = new PDBFileParser();
		Structure structure = pdbpars.parsePDBFile(inStream);
		List<String> keywords = structure.getKeywords();
		assertEquals(keywords.size(), 12);
		assertEquals(keywords.get(11), "TRANSCRIPTION REGULATOR");
	}

	@Test
	public void testDash() throws IOException {
		String fileName;
		fileName = "/pdb6elw-26lines.ent.gz";
		InputStream resourceAsStream = getClass().getResourceAsStream(fileName);
		GZIPInputStream inStream = new GZIPInputStream(resourceAsStream);
		
		Structure structure =  new PDBFileParser().parsePDBFile(inStream);
		
		List<String> keywords = structure.getKeywords();
		assertEquals(keywords.size(), 6);
		assertEquals(keywords.get(3), "THIOREDOXIN-FOLD");
		assertEquals(keywords.get(4), "ANTI-OXIDATVE DEFENSE SYSTEM");
	}
}

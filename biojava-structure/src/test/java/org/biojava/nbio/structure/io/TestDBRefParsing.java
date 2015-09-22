package org.biojava.nbio.structure.io;

import static org.junit.Assert.*;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.junit.Test;

public class TestDBRefParsing {
	@Test
	public void test2W6E() throws IOException, StructureException {

		// an entry with a title in multiple lines in PDB file

		AtomCache cache = new AtomCache();

		StructureIO.setAtomCache(cache);

		cache.setUseMmCif(false);
		Structure sPdb = StructureIO.getStructure("2W6E");
		// System.out.println(sPdb.getName());
		System.out.println(sPdb.getDBRefs().get(0).toPDB());
	}

	@Test
	public void testShortLine() throws IOException,StructureException{
		String shortLine = "DBREF  2W6E A  -42   510  UNP    P19483   ATPA1_BOVIN      1    553";
		InputStream is = new ByteArrayInputStream(shortLine.getBytes());
		PDBFileParser pdbPars = new PDBFileParser();
		Structure s;
		try{
			s=pdbPars.parsePDBFile(is);
		}catch (Exception e){
			is.close();
			throw new AssertionError("Unable to process truncated DBRef line");
		}
		
		
		is = new ByteArrayInputStream(String.format("%1$-80s", shortLine).getBytes());
		Structure ref = pdbPars.parsePDBFile(is);
		is.close();
		assertEquals(ref.getDBRefs().get(0), s.getDBRefs().get(0));
	}
}

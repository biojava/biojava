package org.biojava.nbio.structure.secstruc;

import java.io.IOException;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 * Test the correctness of the DSSP implementation for calculating the
 * secondary structure of a Structure object.
 * 
 * @author Aleix Lafita
 *
 */
public class TestSecStruc {

	@Test
	public void testDSSP() throws IOException, StructureException {
		
			AtomCache cache = new AtomCache();
			cache.getFileParsingParams().setParseSecStruc(true);
			cache.setUseMmCif(false);
			Structure s = cache.getStructure("5pti");

			SecStruc sec = new SecStruc();
			sec.assign(s);
			
			String actual = sec.toString();
			String expected = "";
			
			assertEquals(expected, actual);
	}
}

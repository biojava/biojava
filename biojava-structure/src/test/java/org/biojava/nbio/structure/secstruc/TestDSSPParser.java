package org.biojava.nbio.structure.secstruc;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 * Test all the options for writting and fetching DSSP file
 * formats. Also compare that files match the prediction (assumed
 * to be correct).
 * 
 * @author Aleix Lafita
 *
 */
public class TestDSSPParser {
	
	@Test
	public void testDSSPParser() throws IOException, StructureException {
		
		//List of names to test the DSSP prediction
		List<String> names = Arrays.asList("5pti");
		
		for (String name : names) {
			
			AtomCache cache = new AtomCache();		
			Structure s = cache.getStructure(name);
			
			//Test loading from file
			List<SecStrucState> file = DSSPParser.parseFile(
					"src/test/resources/"+name+".dssp", s, false);
			
			//Test fetching from PDB
			List<SecStrucState> pdb = DSSPParser.fetch(name, s, false);
			
			//Test predicting, writting and parsing back
			SecStrucPred sec = new SecStrucPred();
			List<SecStrucState> pred = sec.predict(s, false);
			
			List<SecStrucState> parseBack = DSSPParser.parseString(
					sec.toString(), s, false);
			
			assertTrue("SS assignment lengths do not match",
					file.size() == pdb.size() && 
					pred.size() == parseBack.size() &&
					pred.size() == file.size());
			
			for (int i=0; i<file.size(); i++){
				assertEquals("SS assignment position "+(i+1)+" does not match",
						file.get(i), pdb.get(i));
				assertEquals("SS assignment position "+(i+1)+" does not match",
						pred.get(i), parseBack.get(i));
				assertEquals("SS assignment position "+(i+1)+" does not match",
						file.get(i), pred.get(i));
			}
		}
	}
}

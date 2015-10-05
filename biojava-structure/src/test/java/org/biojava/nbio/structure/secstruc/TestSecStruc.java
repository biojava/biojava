package org.biojava.nbio.structure.secstruc;

import java.io.IOException;
import java.util.List;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.junit.Before;
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
	
	private Structure s;
	
	@Before
	public void setUp() throws IOException, StructureException {
		
		AtomCache cache = new AtomCache();
		cache.getFileParsingParams().setParseSecStruc(true);
		cache.setUseMmCif(false);
		
		s = cache.getStructure("5pti");
		
	}

	@Test
	public void testDSSPImplementation() throws StructureException {

			SecStrucPred sec = new SecStrucPred();
			sec.predict(s, false);
			
			System.out.println(sec);
	}
	
	@Test
	public void testDSSPParser() throws IOException, StructureException {
		
		List<SecStrucState> secstruc = 
				DSSPParser.parseDSSP("src/test/resources/5pti.dssp", s, false);
		
		for (int i=0; i<secstruc.size(); i++){
			//System.out.println(secstruc.get(i).printDSSPline(i));
		}
	}
	
}

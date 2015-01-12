package org.biojava.bio.structure.io;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.Map;

import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.ResidueNumber;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.quaternary.BioAssemblyInfo;
import org.biojava3.structure.StructureIO;
import org.junit.Test;

/**
 * Testing parsing of some difficult mmCIF files. 
 * For instance those containing multi-line quoting using ";\n" as delimiters
 * Feel free to add any other difficult case here
 * 
 * 
 * @author duarte_j
 *
 */
public class TestDifficultMmCIFFiles {

	@Test
	public void test2BI6() throws IOException, StructureException {
		
		// In this entry _struct_conf contains multiline quoting (quoting with "\n;" ) in a non-loop field
		
		// It seems that at the moment the field is not parsed by the mmCIF parser, anyway let's 
		// keep this here if in the future it is  

		
		AtomCache cache = new AtomCache();
		
		StructureIO.setAtomCache(cache); 

		cache.setUseMmCif(true);
		Structure sCif = StructureIO.getStructure("2BI6");
		
		assertNotNull(sCif);
		
		// an NMR entry
		assertFalse(sCif.isCrystallographic());

		assertTrue(sCif.isNmr());
	
	}
	
	@Test
	public void test1GQO() throws IOException, StructureException {
		
		// In this entry _pdbx_struct_assembly_gen contains multiline quoting (quoting with "\n;" ) in loop field
		
		AtomCache cache = new AtomCache();
		
		StructureIO.setAtomCache(cache); 
				
		FileParsingParameters params = cache.getFileParsingParams();
		params.setParseBioAssembly(true);
		StructureIO.setAtomCache(cache);
	
		cache.setUseMmCif(false);
		Structure sPdb = StructureIO.getStructure("1GQO");
		
		cache.setUseMmCif(true);
		Structure sCif = StructureIO.getStructure("1GQO");
		
		assertNotNull(sCif);

		assertNotNull(sPdb.getPDBHeader().getBioAssemblies());
		assertNotNull(sCif.getPDBHeader().getBioAssemblies());
		
		Map<String,BioAssemblyInfo> mapPdb = sPdb.getPDBHeader().getBioAssemblies();
		Map<String,BioAssemblyInfo> mapCif = sCif.getPDBHeader().getBioAssemblies();
		
		
		
		assertEquals(mapPdb.size(),mapCif.size());
		
		assertEquals(60, mapCif.get("1").getTransforms().size());
		assertEquals(60, mapCif.get("2").getTransforms().size());
		
		// an X-RAY entry
		assertTrue(sPdb.isCrystallographic());
		assertTrue(sCif.isCrystallographic());

		assertFalse(sPdb.isNmr());
		assertFalse(sCif.isNmr());

		
	}

	@Test
	public void testResidueNumbers() throws IOException, StructureException {
		AtomCache cache = new AtomCache();
		cache.setUseMmCif(true);

		Structure s = cache.getStructure("2PTC");
		Chain c = s.getChain(0);

		assertEquals("Wrong first chain",c.getChainID(),"E");

		Group res = c.getAtomGroup(0);
		ResidueNumber resNum = res.getResidueNumber();

		assertEquals("Groups have wrong chain in resnum",resNum.getChainId(),"E");
	}
}

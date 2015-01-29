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
package org.biojava.bio.structure.io;

import static org.junit.Assert.*;
import static org.junit.Assume.*;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.List;
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
		
		Map<Integer,BioAssemblyInfo> mapPdb = sPdb.getPDBHeader().getBioAssemblies();
		Map<Integer,BioAssemblyInfo> mapCif = sCif.getPDBHeader().getBioAssemblies();
		
		
		
		assertEquals(mapPdb.size(),mapCif.size());
		
		assertEquals(60, mapCif.get(1).getTransforms().size());
		assertEquals(60, mapCif.get(2).getTransforms().size());
		
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

	@Test
	public void test4letterChains() throws IOException, StructureException, URISyntaxException {
		String filename = "/1hh0_4char.cif.gz";
		URL url = getClass().getResource(filename);
		assumeNotNull("Can't find resource "+filename,url);

		File file = new File(url.toURI());
		assumeNotNull(file);
		assumeTrue(file.exists());

		MMCIFFileReader reader = new MMCIFFileReader();
		Structure s = reader.getStructure(file);

		assertNotNull("Failed to load structure from jar",s);

		List<Chain> chains = s.getChains();
		assertEquals("Wrong number of chains",chains.size(), 1);

		Chain chain = chains.get(0);
		assertEquals("Wrong chain ID",chain.getChainID(),"ABCD");

		Chain chain2 = s.getChainByPDB("ABCD");
		assertNotNull(chain2);
		assertEquals(chain2, chain);
	}
}

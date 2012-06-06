/**
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
 * Created on Apr 20, 2012
 * Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package org.biojava.bio.structure;

import java.io.IOException;
import java.util.List;

import junit.framework.TestCase;

import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava.bio.structure.io.PDBFileReader;
import org.biojava.bio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.bio.structure.io.mmcif.ChemCompProvider;
import org.biojava.bio.structure.io.mmcif.DownloadChemCompProvider;
import org.biojava.bio.structure.io.mmcif.chem.PolymerType;


/** This class tests the correct loading of Nucleotides
 * 
 * @author Andreas Prlic
 * @since 3.0.3
 */
public class TestNucleotides extends TestCase{

	AtomCache cache = new AtomCache();

	
	public void test3T5N(){
		
		String pdbId = "3T5N";
		Structure s = null;
		try {
			s = getStructure(pdbId);

		} catch (Exception e) {
			e.printStackTrace();
			fail(e.getMessage());

		}

		assertEquals(2,s.getChains().size());

		Chain c = s.getChains().get(1);
		assertEquals("C", c.getChainID());
		List<Group> ngr = c.getAtomGroups("nucleotide");
		assertEquals(6,ngr.size());
		
		
		// now test if we download all definitions correctly for this one...
		PDBFileReader reader = new PDBFileReader();
		reader.setAutoFetch(true);
		FileParsingParameters params = new FileParsingParameters();
		params.setParseSecStruc(true);
		params.setAlignSeqRes(true);
		params.setParseCAOnly(false);
		params.setLoadChemCompInfo(true);
		reader.setFileParsingParameters(params);
		
		ChemCompProvider chemProv = ChemCompGroupFactory.getChemCompProvider();
		try {
			 
			
			DownloadChemCompProvider download = new DownloadChemCompProvider();
			
			ChemCompGroupFactory.setChemCompProvider(download);
			
			Structure s1 = reader.getStructureById(pdbId);
						
			assertNotNull(s1);
			
			assertEquals(2,s1.getChains().size());

			Chain c1 = s1.getChains().get(1);
			
			assertEquals("C", c1.getChainID());
			
			Group g = c1.getAtomGroup(0);
			assertNotNull(g);
			assertNotNull(g.getChemComp());
			assertNotNull(g.getChemComp().getPolymerType());
			assertNotNull(g.getChemComp().getPolymerType().name());
			
			assertTrue("Found an unknown polymertype!", (! g.getChemComp().getPolymerType().equals(PolymerType.unknown)));
			//System.out.println(g.getChemComp().getPolymerType());
			
			List<Group> ngr1 = c1.getAtomGroups("nucleotide");
			assertEquals(6,ngr1.size());
		
		} catch (Exception e){
			e.printStackTrace();
			fail(e.getMessage());
		}
		ChemCompGroupFactory.setChemCompProvider(chemProv);
		
		
	}


	public void test1OFX(){
		Structure s = null;
		try {
			s = getStructure("1OFX");

		} catch (Exception e) {
			e.printStackTrace();
			fail(e.getMessage());

		}

		assertEquals(2,s.getChains().size());

		Chain a = s.getChains().get(0);
		assertEquals("A", a.getChainID());
		List<Group> ngrA = a.getAtomGroups("nucleotide");
		assertEquals(10,ngrA.size());

		Chain b = s.getChains().get(1);
		assertEquals("B", b.getChainID());
		List<Group> ngrB = b.getAtomGroups("nucleotide");
		assertEquals(10,ngrB.size());
	}

	private Structure getStructure(String pdbId) throws IOException, StructureException {
		//System.out.println("cache: " + ChemCompGroupFactory.getChemCompProvider().getClass().getName());
		
		//System.out.println("cache: download chem comp:" + cache.getFileParsingParams().isLoadChemCompInfo());
		return cache.getStructure(pdbId);
	}

	public void test1REP(){
		
		PDBFileReader reader = new PDBFileReader();
		reader.setAutoFetch(true);
		FileParsingParameters params = new FileParsingParameters();
		params.setParseSecStruc(true);
		params.setAlignSeqRes(true);
		params.setParseCAOnly(false);
		params.setLoadChemCompInfo(true);
		reader.setFileParsingParameters(params);
		
		
		try {
			Structure s = reader.getStructureById("1REP");
			//System.out.println(s);
			//System.out.println(s.toPDB());
			Chain b = s.getChainByPDB("B");
			
			assertEquals(22,b.getSeqResGroups().size());
			assertEquals(23,b.getAtomGroups().size());
			
			Group n1 = b.getSeqResGroup(0);
			Group n2 = b.getAtomGroup(0);
			//System.out.println(n1);
			//System.out.println(n2);
			//System.out.println(n1.getChemComp());
			

			assertNotNull("Could not acces Chem Comp file!" , n1.getChemComp());
			assertTrue("ChemComp is not DC",n1.getChemComp().getId().equals("DC"));
			assertNotNull("Could not determine polymer type " , n1.getChemComp().getPolymerType());
			//System.out.println(n1.getChemComp().getPolymerType());
			assertTrue(n1.getChemComp().getPolymerType().equals(PolymerType.dna));
			
			assertNotNull(n1.getPDBName());
			assertNotNull(n1.getResidueNumber());
			assertNotNull(n2.getResidueNumber());
			assertEquals("23", n2.getResidueNumber().toString());
			assertTrue(n1.getResidueNumber().equals(n2.getResidueNumber()));
			
		} catch (Exception e){
			e.printStackTrace();
		}
	}
}
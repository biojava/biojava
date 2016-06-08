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
 * Created on Oct 5, 2009
 * Author: Andreas Prlic
 *
 */

package org.biojava.nbio.structure.test;

import junit.framework.TestCase;
import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.PDBFileParser;
import org.biojava.nbio.structure.io.mmcif.MMcifParser;
import org.biojava.nbio.structure.io.mmcif.SimpleMMcifConsumer;
import org.biojava.nbio.structure.io.mmcif.SimpleMMcifParser;

import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;

public class Test1a4w extends TestCase{

	private static Structure structure = null;


	@Override
	protected void setUp() throws Exception {
		// TODO Auto-generated method stub
		super.setUp();

		if ( structure != null )
			return;

		InputStream inStream = this.getClass().getResourceAsStream("/1a4w.pdb");
		assertNotNull(inStream);

		PDBFileParser pdbpars = new PDBFileParser();
		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(true);

		pdbpars.setFileParsingParameters(params);

		structure = pdbpars.parsePDBFile(inStream) ;

	}



	public void test1a4wPDBFile() throws IOException
	{

		//		structure = null;


		assertNotNull(structure);

		assertEquals("structure does not contain 3 chains ", 3 ,structure.size());

		testStructure(structure);



		Structure structure2 = null;

		InputStream inStream = this.getClass().getResourceAsStream("/1a4w.cif");
		assertNotNull(inStream);

		MMcifParser pdbpars = new SimpleMMcifParser();
		SimpleMMcifConsumer consumer = new SimpleMMcifConsumer();
		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(true);
		consumer.setFileParsingParameters(params);
		pdbpars.addMMcifConsumer(consumer);

		pdbpars.parse(inStream) ;
		structure2 = consumer.getStructure();



		assertNotNull(structure2);

		assertEquals("structure does not contain four chains ", 3 ,structure2.size());

		testStructure(structure2);

		assertEquals(structure.getPDBHeader().toPDB().toLowerCase(),structure2.getPDBHeader().toPDB().toLowerCase());

		for ( int i = 0 ; i < 3 ; i++){
			Chain c1 = structure.getChainByIndex(i);
			Chain c2 = structure.getChainByIndex(i);
			testEqualChains(c1, c2);
		}
	}



	private void testStructure(Structure structure){
		List<Chain> chains = structure.getPolyChains();
		assertEquals("1a4w should have 3 chains. " , 3 , chains.size());

		Chain a = chains.get(0);
		assertEquals("1a4w first chain should be L. " , a.getName(), "L");

		Chain b = chains.get(1);
		assertEquals("1a4w second chain should be H. " , b.getName(), "H");

		Chain c = chains.get(2);
		assertEquals("1a4w third chain should be I. " , c.getName(), "I");

		//System.out.println(structure);
		assertTrue("chain " + a.getName() + " length should be 26. was: " + a.getAtomGroups(GroupType.AMINOACID).size(), ( a.getAtomGroups(GroupType.AMINOACID).size() == 26 ) );

		assertTrue("chain " + a.getName() + " seqres length should be 36. was: " + a.getSeqResLength(), a.getSeqResLength() == 36);

		assertTrue("chain " + b.getName() + " length should be 248. was: " + b.getAtomGroups(GroupType.AMINOACID).size(), ( b.getAtomGroups(GroupType.AMINOACID).size() == 248 ) );

		assertTrue("chain " + b.getName() + " seqres length should be 259. was: " + b.getSeqResLength(), b.getSeqResLength() == 259);

		assertTrue("chain " + c.getName() + " length should be 9. was: " + c.getAtomGroups(GroupType.AMINOACID).size(), ( c.getAtomGroups(GroupType.AMINOACID).size() == 9 ) );

		assertTrue("chain " + c.getName() + " seqres length should be 12. was: " + c.getSeqResLength(), c.getSeqResLength() == 12);

		assertEquals("chain " + c.getName() + " seqres sequences is not correct!", "NGDFEEIPEEYL", c.getSeqResSequence());
	}

	private void testEqualChains(Chain a,Chain b){

		assertEquals("length of seqres " + a.getName() + " and "+b.getName()+" should be same. " , a.getSeqResLength(), b.getSeqResLength() );
		assertEquals("length of atom "   + a.getName() + " and "+b.getName()+" should be same. " , a.getAtomGroups(GroupType.AMINOACID).size(), b.getAtomGroups(GroupType.AMINOACID).size());
		assertEquals("sequences should be identical. " , a.getAtomSequence(),   b.getAtomSequence());
		assertEquals("sequences should be identical. " , a.getSeqResSequence(), b.getSeqResSequence());
	}

	public void testChemComps() throws IOException, StructureException { 
		AtomCache cache = new AtomCache();
		FileParsingParameters params = cache.getFileParsingParams();
		params.setAlignSeqRes(true);
		Structure s = cache.getStructure("1a4w");

		assertEquals(3, s.getPolyChains().size());

		Chain c2 = s.getChainByIndex(1);
		assertEquals("H", c2.getName());

		
		List<Group> hChainLigandGroups = new ArrayList<>();
		
		for (Chain ch : s.getNonPolyChains()) {
			if (ch.getName().equals("H")) {
				hChainLigandGroups.addAll(ch.getAtomGroups());
			}
		}


		boolean noWater = true;
		boolean darPresent = false;


		for ( Group g : hChainLigandGroups){
			String pdbName = g.getPDBName();
			if ( pdbName.equals("QWE"))
				darPresent = true;

			else if ( pdbName.equals("H2O"))
				noWater = false;
		}

		assertTrue("Found water in ligands list!", noWater );

		assertTrue("Did not find QWE in ligands list!", darPresent);

		assertEquals("Did not find the correct nr of ligands in chain! " , 3, hChainLigandGroups.size());

	}

	public void testLigandLoading(){
		Chain c2 = structure.getChainByIndex(1);
		assertEquals("H", c2.getName());

		

		List<Group> hChainLigandGroups = new ArrayList<>();
		
		for (Chain ch : structure.getNonPolyChains()) {
			if (ch.getName().equals("H")) {
				hChainLigandGroups.addAll(ch.getAtomGroups());
			}
		}

		System.out.println("LIGANDS:" + hChainLigandGroups);
		assertEquals("Did not find the correct nr of ligands in chain! " , 6,hChainLigandGroups.size());

		List<Group> lignads2 = StructureTools.filterLigands(hChainLigandGroups);

		assertEquals("Did not get the same nr of ligands from different access methods! ",hChainLigandGroups.size(), lignads2.size());

	}

	public void testSiteGroups(){

		assertNotNull(structure);

		//			Structure s = TmpAtomCache.cache.getStructure("1a4w");

		//                    test1a4wPDBFile();
		Structure s = structure;
		//			for (Chain chain : s.getChains()) {
		//				System.out.println("Chain: " + chain.getChainID());
		//			}
		Chain c2 = s.getChainByIndex(1);
		assertEquals("H", c2.getName());

		//			if (s == null) {
		//				System.out.println("No structure set");
		//			}
		List<Site> sites = s.getSites();
		//System.out.println("sites " + sites);
		assertEquals(7, sites.size());

		boolean noWater = true;
		boolean darPresent = false;
		boolean glyPresent = false;
		Site testSite = null;
		for (Site site : sites) {
			if (site.getSiteID().equals("AC3")) {
				testSite = site;
				for ( Group g : site.getGroups()){
					// before biojava 5.0 we were asserting for same chain, but now the DAR residue is in its own nonpoly chain, can't do it anymore 
					//assertEquals(c2, g.getChain());
					String pdbName = g.getPDBName();
					if ( pdbName.equals("DAR")) {
						darPresent = true;
						//System.out.println("darPresent");
					}
					else if ( pdbName.equals("GLY"))
						glyPresent = true;
					else if ( pdbName.equals("H2O"))
						noWater = false;
				}
			}
		}


		assertTrue("Found water in site list!", noWater );

		assertTrue("Did not find DAR in site list!", darPresent);

		assertTrue("Did not find GLY in site list!", glyPresent);

		//System.out.println(ligands);
		assertEquals("Did not find the correct nr of ligands in chain! " , 8, testSite.getGroups().size());
		

	}

}

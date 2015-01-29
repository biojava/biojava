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

package org.biojava.bio.structure;

import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava.bio.structure.io.PDBFileParser;
import org.biojava.bio.structure.io.mmcif.MMcifParser;
import org.biojava.bio.structure.io.mmcif.SimpleMMcifConsumer;
import org.biojava.bio.structure.io.mmcif.SimpleMMcifParser;
import org.junit.Test;

import java.io.IOException;
import java.io.InputStream;
import java.util.List;
import java.util.zip.GZIPInputStream;

import static org.junit.Assert.*;

/** A collection of tests to make sure 4hhb is represented correctly.
 *  
 * @author Andreas Prlic
 *
 */
public class Test4hhb {

	@Test
	public void test4hhbPDBFile() throws IOException
	{

		Structure structure = null;

		InputStream inStream = new GZIPInputStream(this.getClass().getResourceAsStream("/4hhb.pdb.gz"));
		assertNotNull(inStream);

		PDBFileParser pdbpars = new PDBFileParser();
		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(true);
		pdbpars.setFileParsingParameters(params);
		//pdbpars.setLoadChemCompInfo(true);

		structure = pdbpars.parsePDBFile(inStream) ;


		assertNotNull(structure);

		assertEquals("structure does not contain four chains ", 4 ,structure.size());

		testStructure(structure);


		Structure structure2 = null;

		inStream = new GZIPInputStream(this.getClass().getResourceAsStream("/4hhb.cif.gz"));
		assertNotNull(inStream);

		MMcifParser mmcifpars = new SimpleMMcifParser();
		SimpleMMcifConsumer consumer = new SimpleMMcifConsumer();
		params = new FileParsingParameters();
		params.setAlignSeqRes(true);
		consumer.setFileParsingParameters(params);
		mmcifpars.addMMcifConsumer(consumer);

		mmcifpars.parse(inStream) ;
		structure2 = consumer.getStructure();
		

		assertNotNull(structure2);

		assertEquals("structure does not contain four chains ", 4 ,structure2.size());

		testStructure(structure2);
		
		assertEquals(structure.getPDBHeader().toPDB().toLowerCase(),structure2.getPDBHeader().toPDB().toLowerCase());

		for ( int i = 0 ; i < 4 ; i++){
			Chain c1 = structure.getChain(i);
			Chain c2 = structure.getChain(i);
			testEqualChains(c1, c2);
		}
		
		// test crystallographic info parses correctly
		testCrystallographicInfo(structure, structure2);
		
	}


	private void testStructure(Structure structure){

		List<Chain> chains = structure.getChains();
		assertEquals("4HHB should have 4 chains. " , 4 , chains.size());

		Chain a = chains.get(0);
		assertEquals("4HHB first chain should be A. " , a.getChainID(), "A");

		Chain b = chains.get(1);
		assertEquals("4HHB second chain should be B. " , b.getChainID(), "B");

		Chain c = chains.get(2);
		assertEquals("4HHB third chain should be C. " , c.getChainID(), "C");

		Chain d = chains.get(3);
		assertEquals("4HHB fourth chain should be D. " , d.getChainID(), "D");

		assertTrue("chain " + a.getChainID() + " length should be 141. was: " + a.getAtomGroups(GroupType.AMINOACID).size(), ( a.getAtomGroups(GroupType.AMINOACID).size() == 141 ));
		assertTrue("chain " + b.getChainID() + " length should be 146. was: " + b.getAtomGroups(GroupType.AMINOACID).size(), ( b.getAtomGroups(GroupType.AMINOACID).size() == 146 ));		
		assertTrue("chain " + c.getChainID() + " length should be 141. was: " + c.getAtomGroups(GroupType.AMINOACID).size(), ( c.getAtomGroups(GroupType.AMINOACID).size() == 141 ));
		assertTrue("chain " + d.getChainID() + " length should be 146. was: " + d.getAtomGroups(GroupType.AMINOACID).size(), ( d.getAtomGroups(GroupType.AMINOACID).size() == 146 ));

		assertTrue("chain " + a.getChainID() + " length should be 141, but is " + a.getSeqResLength(), ( a.getSeqResLength() == 141 ));
		assertTrue("chain " + b.getChainID() + " length should be 146.", ( b.getSeqResLength() == 146 ));

		assertTrue("chain " + c.getChainID() + " length should be 141.", ( a.getSeqResLength() == 141 ));
		assertTrue("chain " + d.getChainID() + " length should be 146.", ( b.getSeqResLength() == 146 ));

		testEqualChains(a,c);

		testEqualChains(b,d);

		Chain[] chs = new Chain[]{a,b,c,d};

		for (Chain x : chs){
			testContainsHem(x);
		}

	}

	private void testContainsHem(Chain x) {

		boolean containsHem = false;
		for ( Group g : x.getAtomGroups()){
			if ( g.getPDBName().equals("HEM")){
				containsHem = true;
				break;
			}
		}

		assertTrue("Chain " + x.getChainID() + " does not contain a HEM group. " , containsHem);

	}


	private void testEqualChains(Chain a,Chain b){

		assertEquals("length of seqres " + a.getChainID() + " and "+b.getChainID()+" should be same. " , a.getSeqResLength(), b.getSeqResLength() );
		assertEquals("length of atom "   + a.getChainID() + " and "+b.getChainID()+" should be same. " , a.getAtomGroups(GroupType.AMINOACID).size(), b.getAtomGroups(GroupType.AMINOACID).size());
		assertEquals("sequences should be identical. " , a.getAtomSequence(),   b.getAtomSequence());
		assertEquals("sequences should be identical. " , a.getSeqResSequence(), b.getSeqResSequence());
	}
	
	private void testCrystallographicInfo(Structure s1, Structure s2) {
		PDBCrystallographicInfo xtalInfo = s1.getPDBHeader().getCrystallographicInfo();
		PDBCrystallographicInfo xtalInfo2 = s2.getPDBHeader().getCrystallographicInfo();
		assertEquals(xtalInfo.getA(), xtalInfo2.getA(), 0.0001);
		assertEquals(xtalInfo.getB(), xtalInfo2.getB(), 0.0001);
		assertEquals(xtalInfo.getC(), xtalInfo2.getC(), 0.0001);
		assertEquals(xtalInfo.getAlpha(), xtalInfo2.getAlpha(), 0.0001);
		assertEquals(xtalInfo.getBeta(), xtalInfo2.getBeta(), 0.0001);
		assertEquals(xtalInfo.getGamma(), xtalInfo2.getGamma(), 0.0001);

		assertEquals(xtalInfo.getSpaceGroup(),xtalInfo2.getSpaceGroup());
	}

}

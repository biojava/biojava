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
package org.biojava.nbio.structure;

import org.biojava.nbio.structure.chem.ChemComp;
import org.biojava.nbio.structure.chem.ChemCompGroupFactory;
import org.biojava.nbio.structure.chem.ChemCompProvider;
import org.biojava.nbio.structure.chem.DownloadChemCompProvider;
import org.biojava.nbio.structure.chem.PolymerType;
import org.biojava.nbio.structure.chem.ReducedChemCompProvider;
import org.biojava.nbio.structure.chem.ResidueType;
import org.junit.Test;

import static org.junit.Assert.*;

public class ChemCompTest {

	@Test
	public void testALA(){
		String chemID = "ALA";

		ChemComp cc = ChemCompGroupFactory.getChemComp(chemID);

		assertNotNull(cc.getPolymerType());

		assertEquals(cc.getPolymerType(), PolymerType.peptide);

		assertNotNull("residue type has not been set!", cc.getResidueType());

		assertTrue (" is not amino ", cc.getResidueType().equals(ResidueType.lPeptideLinking));
	}

	@Test
	public void testMEA(){


		String chemID = "MEA";



		// also test replacing providers...
		ChemCompProvider oldProvider = ChemCompGroupFactory.getChemCompProvider();

		DownloadChemCompProvider all = new DownloadChemCompProvider();

		ChemCompGroupFactory.setChemCompProvider(all);

		ChemComp cc = ChemCompGroupFactory.getChemComp(chemID);

		assertNotNull(cc);

		assertTrue(" is not mea" , cc.getId().equals(chemID));

		assertEquals(" one letter code is not correct", "F", cc.getOneLetterCode());

		assertEquals("MEA",cc.getThreeLetterCode());

		assertNotNull(cc.getPolymerType());

		assertEquals( PolymerType.peptide, cc.getPolymerType());

		assertNotNull("residue type has not been set!", cc.getResidueType());

		assertTrue (" is not amino ", cc.getResidueType().equals(ResidueType.lPeptideLinking));

		Group g = ChemCompGroupFactory.getGroupFromChemCompDictionary(chemID);

		assertTrue( g.getType().equals(GroupType.AMINOACID));

		ChemCompGroupFactory.setChemCompProvider(oldProvider);
	}

	@Test
	public void testPRR(){

		ChemCompProvider oldProvider = ChemCompGroupFactory.getChemCompProvider();

		DownloadChemCompProvider all = new DownloadChemCompProvider();

		ChemCompGroupFactory.setChemCompProvider(all);

		String chemID = "PRR";

		Group g = ChemCompGroupFactory.getGroupFromChemCompDictionary(chemID);

		assertTrue("Got back group of wrong type! " + g.getClass().getName(),  g instanceof AminoAcid);

		AminoAcid aa = (AminoAcid) g;

		assertNotNull(aa.getAminoType());

		assertTrue(aa.getAminoType().equals('X'));

		ChemCompGroupFactory.setChemCompProvider(oldProvider);

	}

	@Test
	public void testChangingProviders(){

		// test for issue #145

	String chemID = "MEA";

	// first we test with reduced chem comp provider
	ChemCompGroupFactory.setChemCompProvider(new ReducedChemCompProvider());

	ChemComp cc = ChemCompGroupFactory.getChemComp(chemID);

	assertNotNull(cc);

	assertTrue(" is not mea" , cc.getId().equals(chemID));

	assertTrue(cc.isEmpty());

	// now we change to download chem comp provider

	ChemCompGroupFactory.setChemCompProvider(new DownloadChemCompProvider());

	cc = ChemCompGroupFactory.getChemComp(chemID);

	assertNotNull(cc);

	assertTrue(" is not mea" , cc.getId().equals(chemID));

	assertEquals("MEA",cc.getThreeLetterCode());



	// now testing in opposite order

	// first we test with download chem comp provider

	ChemCompGroupFactory.setChemCompProvider(new DownloadChemCompProvider());

	cc = ChemCompGroupFactory.getChemComp(chemID);

	assertNotNull(cc);

	assertTrue(" is not mea" , cc.getId().equals(chemID));

	assertEquals("MEA",cc.getThreeLetterCode());

	// now we change to reduced chem comp provider
	ChemCompGroupFactory.setChemCompProvider(new ReducedChemCompProvider());

	cc = ChemCompGroupFactory.getChemComp(chemID);

	assertNotNull(cc);

	assertTrue(" is not mea" , cc.getId().equals(chemID));

	//the cached description contains all information even with the ReducedProvider
	assertNotNull(cc.getThreeLetterCode());


	}

}

package org.biojava.bio.structure;

import org.biojava.bio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.bio.structure.io.mmcif.ChemCompProvider;
import org.biojava.bio.structure.io.mmcif.DownloadChemCompProvider;
import org.biojava.bio.structure.io.mmcif.chem.PolymerType;
import org.biojava.bio.structure.io.mmcif.chem.ResidueType;
import org.biojava.bio.structure.io.mmcif.model.ChemComp;

import junit.framework.TestCase;

public class ChemCompTest extends TestCase {

	public void testALA(){
		String chemID = "ALA";

		ChemComp cc = ChemCompGroupFactory.getChemComp(chemID);

		assertNotNull(cc.getPolymerType());

		assertEquals(cc.getPolymerType(), PolymerType.peptide);

		assertNotNull("residue type has not been set!", cc.getResidueType());

		assertTrue (" is not amino ", cc.getResidueType().equals(ResidueType.lPeptideLinking));
	}

	public void testMEA(){
		
		
		String chemID = "MEA";
		
		ChemCompProvider oldProvider = ChemCompGroupFactory.getChemCompProvider();
		
		DownloadChemCompProvider all = new DownloadChemCompProvider();
		//all.setDownloadAll(true);
		
		ChemCompGroupFactory.setChemCompProvider(all);

		ChemComp cc = ChemCompGroupFactory.getChemComp(chemID);

		assertNotNull(cc);
		
		assertTrue(" is not mea" , cc.getId().equals(chemID));

		assertEquals(" one letter code is not correct", cc.getOne_letter_code(),"F");

		assertEquals(cc.getThree_letter_code(),"MEA");

		assertNotNull(cc.getPolymerType());

		assertEquals(cc.getPolymerType(), PolymerType.peptide);

		assertNotNull("residue type has not been set!", cc.getResidueType());

		assertTrue (" is not amino ", cc.getResidueType().equals(ResidueType.lPeptideLinking));

		Group g = ChemCompGroupFactory.getGroupFromChemCompDictionary(chemID);

		assertTrue( g.getType().equals("amino"));

		ChemCompGroupFactory.setChemCompProvider(oldProvider);
	}
	
	public void testPRR(){
		
		String chemID = "PRR"; 
				
		Group g = ChemCompGroupFactory.getGroupFromChemCompDictionary(chemID);
		
		assertEquals(g.getType(), AminoAcidImpl.type);
		
		AminoAcid aa = (AminoAcid) g;
		
		assertNotNull(aa.getAminoType());
		
		assertTrue(aa.getAminoType().equals('X'));
		
	}

}

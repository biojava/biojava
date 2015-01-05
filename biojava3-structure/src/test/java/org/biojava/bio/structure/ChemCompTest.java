package org.biojava.bio.structure;

import org.biojava.bio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.bio.structure.io.mmcif.ChemCompProvider;
import org.biojava.bio.structure.io.mmcif.DownloadChemCompProvider;
import org.biojava.bio.structure.io.mmcif.ReducedChemCompProvider;
import org.biojava.bio.structure.io.mmcif.chem.PolymerType;
import org.biojava.bio.structure.io.mmcif.chem.ResidueType;
import org.biojava.bio.structure.io.mmcif.model.ChemComp;
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

		assertEquals(" one letter code is not correct", "F", cc.getOne_letter_code());

		assertEquals("MEA",cc.getThree_letter_code());

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

        // an empty description is returned as expected
        assertNull(cc.getThree_letter_code());


        // now we change to download chem comp provider

        ChemCompGroupFactory.setChemCompProvider(new DownloadChemCompProvider());

        cc = ChemCompGroupFactory.getChemComp(chemID);

        assertNotNull(cc);

        assertTrue(" is not mea" , cc.getId().equals(chemID));

        assertEquals("MEA",cc.getThree_letter_code());

        
        
        // now testing in opposite order
        

        // first we test with download chem comp provider

        ChemCompGroupFactory.setChemCompProvider(new DownloadChemCompProvider());

        cc = ChemCompGroupFactory.getChemComp(chemID);

        assertNotNull(cc);

        assertTrue(" is not mea" , cc.getId().equals(chemID));

        assertEquals("MEA",cc.getThree_letter_code());
        
        // now we change to reduced chem comp provider 
        ChemCompGroupFactory.setChemCompProvider(new ReducedChemCompProvider());

        cc = ChemCompGroupFactory.getChemComp(chemID);

        assertNotNull(cc);

        assertTrue(" is not mea" , cc.getId().equals(chemID));

        // an empty description is returned as expected
        assertNull(cc.getThree_letter_code());

        
    }

}

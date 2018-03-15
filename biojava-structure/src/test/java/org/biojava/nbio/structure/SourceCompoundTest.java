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
 * created at Apr 5, 2008
 */

package org.biojava.nbio.structure;

import org.biojava.nbio.structure.io.PDBFileParser;
import org.junit.Assert;
import org.junit.Test;

import java.io.IOException;
import java.io.InputStream;
import java.util.List;


public class SourceCompoundTest {

	private Structure getStructure(String fileName){

		InputStream inStream = this.getClass().getResourceAsStream(fileName);
		Assert.assertNotNull(inStream);

		PDBFileParser pdbpars = new PDBFileParser();
		Structure structure = null;
		try {
			structure = pdbpars.parsePDBFile(inStream) ;
		} catch (IOException e) {
			e.printStackTrace();
		}
		return structure;
	}


	@Test
	public void testCompoundSourceStructure(){

		Structure s2 = getStructure("/2gox.pdb");
		// note since biojava 5.0 we are finding entities for all molecules even if
		// the annotation is not present for them
		// 2 protein entities and 1 water entity
		Assert.assertEquals(3, s2.getEntityInfos().size());
		for (EntityInfo compound : s2.getEntityInfos()){
			if (compound.getMolId()==1) {
				Assert.assertEquals("COMPLEMENT C3", compound.getDescription());
				Assert.assertEquals("[A, C]", compound.getChainIds().toString());
				Assert.assertEquals("FRAGMENT OF ALPHA CHAIN: RESIDUES 996-1287", compound.getFragment());
				Assert.assertEquals("YES", compound.getEngineered());
				Assert.assertEquals("YES", compound.getMutation());
				Assert.assertEquals("HOMO SAPIENS", compound.getOrganismScientific());
				Assert.assertEquals("HUMAN", compound.getOrganismCommon());
				Assert.assertEquals("C3", compound.getGene());
				Assert.assertEquals("ESCHERICHIA COLI", compound.getExpressionSystem());
				Assert.assertEquals("BL21(DE3)", compound.getExpressionSystemStrain());
				Assert.assertEquals("PLASMID", compound.getExpressionSystemVectorType());
				Assert.assertEquals("PT7-", compound.getExpressionSystemPlasmid());
			}
			if (compound.getMolId()==2) {
				Assert.assertEquals("FIBRINOGEN-BINDING PROTEIN", compound.getDescription());
				Assert.assertEquals("[B, D]", compound.getChainIds().toString());
				Assert.assertEquals("C-TERMINAL DOMAIN: RESIDUES 101-165", compound.getFragment());
				Assert.assertEquals("YES", compound.getEngineered());
				Assert.assertEquals("STAPHYLOCOCCUS AUREUS", compound.getOrganismScientific());
				Assert.assertEquals("BACTERIA", compound.getOrganismCommon());
				Assert.assertEquals("MU50 / ATCC 700699", compound.getStrain());
				Assert.assertEquals("EFB", compound.getGene());
				Assert.assertEquals("ESCHERICHIA COLI", compound.getExpressionSystem());
				Assert.assertEquals("BL21(DE3)", compound.getExpressionSystemStrain());
				Assert.assertEquals("PLASMID", compound.getExpressionSystemVectorType());
				Assert.assertEquals("PT7HMT", compound.getExpressionSystemPlasmid());
			}
		}

	}

	@Test
	public void testCOMPNDsectionFRAGMENT(){
		Structure s2 = getStructure("/2gox.pdb");
		Structure s4 = getStructure("/3cfy.pdb");

		// this file has a CHAIN: string in the value for the FRAGMENT: filed which breaks the version 1.4 parser

		for (EntityInfo compound : s2.getEntityInfos()) {
			if (compound.getMolId()==1) {
				Assert.assertEquals("FRAGMENT OF ALPHA CHAIN: RESIDUES 996-1287", compound.getFragment());
			}

		}

		for (EntityInfo compound : s4.getEntityInfos()) {
			if (compound.getMolId()==1) {
				Assert.assertEquals("SIGNAL RECEIVER DOMAIN: RESIDUES 2-128", compound.getFragment());
			}

		}

	}

	@Test
	public void testCOMPNDsectionCHAINS(){
		Structure s3 = getStructure("/2pos.pdb");
		// note since biojava 5.0 we are finding entities for all molecules even if
		// the annotation is not present for them
		// thus for 2pos.pdb we have 1 protein entity, but 3 non-polymer entities and 1 water entity
		EntityInfo compound = s3.getEntityById(1);
		Assert.assertEquals(5, s3.getEntityInfos().size());
		Assert.assertEquals(1, compound.getMolId());
		Assert.assertEquals("SYLVATICIN", compound.getDescription());
		Assert.assertEquals("[A, B, C, D]", compound.getChainIds().toString());
		Assert.assertEquals("PYTHIUM SYLVATICUM", compound.getOrganismScientific());
		Assert.assertEquals("STRAIN 37", compound.getStrain());

	}

	@Test
	public void testSOURCEsectionSTRAIN(){
		Structure s4 = getStructure("/3cfy.pdb");
		for (EntityInfo compound : s4.getEntityInfos()){
			if (compound.getMolId()==1) {
				/*System.out.println(compound.getMolId());
				System.out.println(compound.getMolName());
				System.out.println(compound.getChainName().toString());
				System.out.println(compound.getFragment());
				System.out.println(compound.getEngineered());
				System.out.println(compound.getOrganismScientific());
				System.out.println(compound.getOrganismCommon());
				System.out.println(compound.getStrain());
				System.out.println(compound.getGene());
				System.out.println(compound.getExpressionSystem());
				System.out.println(compound.getExpressionSystemVectorType());
				System.out.println(compound.getExpressionSystemVector());
				System.out.println(compound.getExpressionSystemPlasmid());
				 */
				Assert.assertEquals(1, compound.getMolId());
				Assert.assertEquals("PUTATIVE LUXO REPRESSOR PROTEIN", compound.getDescription());
				Assert.assertEquals("[A]", compound.getChainIds().toString());
				Assert.assertEquals("SIGNAL RECEIVER DOMAIN: RESIDUES 2-128", compound.getFragment());
				Assert.assertEquals("YES", compound.getEngineered());
				Assert.assertEquals("VIBRIO PARAHAEMOLYTICUS RIMD 2210633", compound.getOrganismScientific());
				Assert.assertEquals("BACTERIA", compound.getOrganismCommon());
				Assert.assertEquals("RIMD 2210633 / SEROTYPE O3:K6", compound.getStrain());
				Assert.assertEquals("VP1469", compound.getGene());
				Assert.assertEquals("ESCHERICHIA COLI", compound.getExpressionSystem());
				Assert.assertEquals("PLASMID", compound.getExpressionSystemVectorType());
				Assert.assertEquals("PET", compound.getExpressionSystemVector());
				Assert.assertEquals("BC-PSGX3(BC)", compound.getExpressionSystemPlasmid());

			}
		}
	}

	@Test
	public void testSOURCEsectionORGSCI(){
		Structure s5 = getStructure("/3cdl.pdb");
		for (EntityInfo compound : s5.getEntityInfos()){
			if (compound.getMolId()==1) {
				//System.out.println(compound.getOrganismScientific());
				Assert.assertEquals("PSEUDOMONAS SYRINGAE PV. TOMATO STR. DC3000", compound.getOrganismScientific());
			}
		}
	}

	/**
	 * There is a file format change in v3.2 of the PDB file format, adding the
	 * tax id.
	 * This test makes sure that the tax id for the organism and expression
	 * systems is set correctly.
	 */
	@Test
	public void testSourceTaxIdVersion32File(){
		Structure structure = getStructure("/3dl7_v32.pdb");

		EntityInfo comp = structure.getEntityById(1);

		Assert.assertEquals("10090", comp.getOrganismTaxId());
		Assert.assertEquals("9606", comp.getExpressionSystemTaxId());

	}

	/**
	 * 3.2 format includes PMID and DOI in the JRNL section.
	 */
	@Test
	public void testJournalRefs(){
//        JRNL        AUTH   M.HAMMEL,G.SFYROERA,D.RICKLIN,P.MAGOTTI,
//        JRNL        AUTH 2 J.D.LAMBRIS,B.V.GEISBRECHT
//        JRNL        TITL   A STRUCTURAL BASIS FOR COMPLEMENT INHIBITION BY
//        JRNL        TITL 2 STAPHYLOCOCCUS AUREUS.
//        JRNL        REF    NAT.IMMUNOL.                  V.   8   430 2007
//        JRNL        REFN                   ISSN 1529-2908
//        JRNL        PMID   17351618
//        JRNL        DOI    10.1038/NI1450
		Structure structure = getStructure("/2gox_v315.pdb");
		//check that there really is an publication
		Assert.assertTrue(structure.hasJournalArticle());

		if (structure.hasJournalArticle()) {
			JournalArticle journal = structure.getJournalArticle();
			List<Author> authorList = journal.getAuthorList();
			Author firstAuthor = authorList.get(0);
			//check the authors
			Assert.assertEquals(6, authorList.size());
			Assert.assertEquals("HAMMEL", firstAuthor.getSurname());
			Assert.assertEquals("M.", firstAuthor.getInitials());
			//check the other publication details
			Assert.assertEquals("A STRUCTURAL BASIS FOR COMPLEMENT INHIBITION BY STAPHYLOCOCCUS AUREUS.", journal.getTitle());
			Assert.assertEquals("NAT.IMMUNOL.", journal.getJournalName());
			Assert.assertEquals(2007, journal.getPublicationDate());
			Assert.assertEquals("8", journal.getVolume());
			Assert.assertEquals("430", journal.getStartPage());
			Assert.assertEquals("ISSN 1529-2908", journal.getRefn());
			Assert.assertEquals("17351618", journal.getPmid());
			Assert.assertEquals("10.1038/NI1450", journal.getDoi());
		}
	}
}

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

package org.biojava.bio.structure;

import java.io.IOException;
import java.io.InputStream;

import java.util.List;

import org.biojava.bio.structure.io.PDBFileParser;


import junit.framework.TestCase;


public class SourceCompoundTest extends TestCase{

	private Structure getStructure(String fileName){

		InputStream inStream = this.getClass().getResourceAsStream(fileName);
		assertNotNull(inStream);

		PDBFileParser pdbpars = new PDBFileParser();
		Structure structure = null;
		try {
			structure = pdbpars.parsePDBFile(inStream) ;
		} catch (IOException e) {
			e.printStackTrace();
		}
		return structure;
	}


	public void testCompoundSourceStructure(){

		Structure s2 = getStructure("/2gox.pdb");
		assertEquals(2, s2.getCompounds().size());
		for (Compound compound : s2.getCompounds()){
			if (compound.getMolId().equals("1")) {
				assertEquals("COMPLEMENT C3", compound.getMolName());
				assertEquals("[A, C]", compound.getChainId().toString());
				assertEquals("FRAGMENT OF ALPHA CHAIN: RESIDUES 996-1287", compound.getFragment());
				assertEquals("YES", compound.getEngineered());
				assertEquals("YES", compound.getMutation());
				assertEquals("HOMO SAPIENS", compound.getOrganismScientific());
				assertEquals("HUMAN", compound.getOrganismCommon());
				assertEquals("C3", compound.getGene());
				assertEquals("ESCHERICHIA COLI", compound.getExpressionSystem());
				assertEquals("BL21(DE3)", compound.getExpressionSystemStrain());
				assertEquals("PLASMID", compound.getExpressionSystemVectorType());
				assertEquals("PT7-", compound.getExpressionSystemPlasmid());
			}
			if (compound.getMolId().equals("2")) {
				assertEquals("FIBRINOGEN-BINDING PROTEIN", compound.getMolName());
				assertEquals("[B, D]", compound.getChainId().toString());
				assertEquals("C-TERMINAL DOMAIN: RESIDUES 101-165", compound.getFragment());
				assertEquals("YES", compound.getEngineered());
				assertEquals("STAPHYLOCOCCUS AUREUS", compound.getOrganismScientific());
				assertEquals("BACTERIA", compound.getOrganismCommon());
				assertEquals("MU50 / ATCC 700699", compound.getStrain());
				assertEquals("EFB", compound.getGene());
				assertEquals("ESCHERICHIA COLI", compound.getExpressionSystem());
				assertEquals("BL21(DE3)", compound.getExpressionSystemStrain());
				assertEquals("PLASMID", compound.getExpressionSystemVectorType());
				assertEquals("PT7HMT", compound.getExpressionSystemPlasmid());
			}
		}

	}

	public void testCOMPNDsectionFRAGMENT(){
		Structure s2 = getStructure("/2gox.pdb");
		Structure s4 = getStructure("/3cfy.pdb");

		// this file has a CHAIN: string in the value for the FRAGMENT: filed which breaks the version 1.4 parser

		for (Compound compound : s2.getCompounds()) {
			if (compound.getMolId().equals("1")) {
				assertEquals("FRAGMENT OF ALPHA CHAIN: RESIDUES 996-1287", compound.getFragment());
			}

		}

		for (Compound compound : s4.getCompounds()) {
			if (compound.getMolId().equals("1")) {
				assertEquals("SIGNAL RECEIVER DOMAIN: RESIDUES 2-128", compound.getFragment());
			}

		}

	}

	public void testCOMPNDsectionCHAINS(){
		Structure s3 = getStructure("/2pos.pdb");
		assertEquals(1, s3.getCompounds().size());
		for (Compound compound : s3.getCompounds()){
			/*System.out.println(compound.getMolId());
			System.out.println(compound.getMolName());
			System.out.println(compound.getChainId().toString());
			System.out.println(compound.getOrganismScientific());
			System.out.println(compound.getStrain());
	*/
			assertEquals("1", compound.getMolId());
			assertEquals("SYLVATICIN", compound.getMolName());
			assertEquals("[A, B, C, D]", compound.getChainId().toString());
			assertEquals("PYTHIUM SYLVATICUM", compound.getOrganismScientific());
			assertEquals("STRAIN 37", compound.getStrain());

		}
	}

	public void testSOURCEsectionSTRAIN(){
		Structure s4 = getStructure("/3cfy.pdb");
		for (Compound compound : s4.getCompounds()){
			if (compound.getMolId().equals("1")) {
				/*System.out.println(compound.getMolId());
				System.out.println(compound.getMolName());
				System.out.println(compound.getChainId().toString());
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
				assertEquals("1", compound.getMolId());
				assertEquals("PUTATIVE LUXO REPRESSOR PROTEIN", compound.getMolName());
				assertEquals("[A]", compound.getChainId().toString());
				assertEquals("SIGNAL RECEIVER DOMAIN: RESIDUES 2-128", compound.getFragment());
				assertEquals("YES", compound.getEngineered());
				assertEquals("VIBRIO PARAHAEMOLYTICUS RIMD 2210633", compound.getOrganismScientific());
				assertEquals("BACTERIA", compound.getOrganismCommon());
				assertEquals("RIMD 2210633 / SEROTYPE O3:K6", compound.getStrain());
				assertEquals("VP1469", compound.getGene());
				assertEquals("ESCHERICHIA COLI", compound.getExpressionSystem());
				assertEquals("PLASMID", compound.getExpressionSystemVectorType());
				assertEquals("PET", compound.getExpressionSystemVector());
				assertEquals("BC-PSGX3(BC)", compound.getExpressionSystemPlasmid());

			}
		}
	}

	public void testSOURCEsectionORGSCI(){
		Structure s5 = getStructure("/3cdl.pdb");
		for (Compound compound : s5.getCompounds()){
			if (compound.getMolId().equals("1")) {
				//System.out.println(compound.getOrganismScientific());
				assertEquals("PSEUDOMONAS SYRINGAE PV. TOMATO STR. DC3000", compound.getOrganismScientific());
			}
		}
	}

   /**
     * There is a file format change in v3.2 of the PDB file format, adding the
     * tax id.
	 * This test makes sure that the tax id for the organism and expression
     * systems is set correctly.
	 */
	public void testSourceTaxIdVersion32File(){
		Structure structure = getStructure("/3dl7_v32.pdb");

        Compound comp = structure.getCompoundById("1");

        comp.showSource();

        assertEquals("10090", comp.getOrganismTaxId());
        assertEquals("9606", comp.getExpressionSystemTaxId());

	}

    /**
     * 3.2 format includes PMID and DOI in the JRNL section.
     */
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
        assertTrue(structure.hasJournalArticle());

        if (structure.hasJournalArticle()) {
            JournalArticle journal = structure.getJournalArticle();
            List<Author> authorList = journal.getAuthorList();
            Author firstAuthor = authorList.get(0);
            //check the authors
            assertEquals(6, authorList.size());
            assertEquals("HAMMEL", firstAuthor.getSurname());
            assertEquals("M.", firstAuthor.getInitials());
            //check the other publication details
            assertEquals("A STRUCTURAL BASIS FOR COMPLEMENT INHIBITION BY STAPHYLOCOCCUS AUREUS.", journal.getTitle());               
            assertEquals("NAT.IMMUNOL.", journal.getJournalName());
            assertEquals(2007, journal.getPublicationDate());
            assertEquals("8", journal.getVolume());
            assertEquals("430", journal.getStartPage());
            assertEquals("ISSN 1529-2908", journal.getRefn());
            assertEquals("17351618", journal.getPmid());
            assertEquals("10.1038/NI1450", journal.getDoi());
        }
    }
}

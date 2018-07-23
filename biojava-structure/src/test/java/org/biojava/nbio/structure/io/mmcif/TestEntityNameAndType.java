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
package org.biojava.nbio.structure.io.mmcif;

import static org.junit.Assert.assertArrayEquals;

import java.io.IOException;

import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.EntityInfo;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.junit.Test;
/**
 * Test to ensure the entity name and type
 * @author Anthony Bradley
 *
 */
public class TestEntityNameAndType {

	@Test	
	public void testEntityId() throws IOException, StructureException {
		
		// Set up the atom cache to parse on Internal chain id
		AtomCache cache = new AtomCache();
		cache.setUseMmCif(true);
		FileParsingParameters params = cache.getFileParsingParams();

		DownloadChemCompProvider cc = new DownloadChemCompProvider();
		ChemCompGroupFactory.setChemCompProvider(cc);
		cc.checkDoFirstInstall();
		cache.setFileParsingParams(params);
		StructureIO.setAtomCache(cache);
		// This is hte information we want to test against
		String[] typeInformation = new String[] {"POLYMER", "NONPOLYMER", "NONPOLYMER", "NONPOLYMER", "NONPOLYMER", "WATER"};
		String[] descriptionInformation = new String[] {"BROMODOMAIN ADJACENT TO ZINC FINGER DOMAIN PROTEIN 2B","4-Fluorobenzamidoxime",  "METHANOL", "METHANOL", "METHANOL", "water"};	
		
		// Now some other information fields to test this data is collated correctly
		String[] geneSourceSciName = new String[] {"HOMO SAPIENS", null, null, null, null, null};
		String[] geneSourceTaxId = new String[] {"9606", null, null, null, null, null};
		String[] hostOrganismSciName = new String[] {"ESCHERICHIA COLI", null, null, null, null, null};
		String[] hostOrganismTaxId = new String[] {"469008", null, null, null, null, null};		

		
		
		/// TODO GET ALL THE ENTITY INFORMATION REQUIRED FOR 4CUP 
		// Get this structure
		Structure bioJavaStruct = StructureIO.getStructure("4cup");
		String[] testTypeInfo = new String[6];
		String[] testDescInfo = new String[6];
		
		
		String[] testGeneSourceSciName = new String[6];
		String[] testGeneSourceTaxId = new String[6];
		String[] testHostOrganismSciName = new String[6];
		String[] testHostOrganismTaxId = new String[6];

		// Now loop through the structure
		int chainCounter = 0;
		for (Chain c: bioJavaStruct.getChains()) {
			// Now get the entity information we want to test
			EntityInfo thisCmpd = c.getEntityInfo();
			testTypeInfo[chainCounter] = thisCmpd.getType().toString();
			testDescInfo[chainCounter] = thisCmpd.getDescription();
			testGeneSourceSciName[chainCounter] =  thisCmpd.getOrganismScientific();
			testGeneSourceTaxId[chainCounter] = thisCmpd.getOrganismTaxId();
			testHostOrganismSciName[chainCounter] = thisCmpd.getExpressionSystem();
			testHostOrganismTaxId[chainCounter] = thisCmpd.getExpressionSystemTaxId();

			chainCounter++;
		}
		// Now check they're both the same
		assertArrayEquals(descriptionInformation, testDescInfo);
		assertArrayEquals(typeInformation, testTypeInfo);
		// Now check these work too
		assertArrayEquals(geneSourceSciName, testGeneSourceSciName);
		assertArrayEquals(geneSourceTaxId, testGeneSourceTaxId);
		assertArrayEquals(hostOrganismSciName, testHostOrganismSciName);
		assertArrayEquals(hostOrganismTaxId, testHostOrganismTaxId);

		
		
	}
}

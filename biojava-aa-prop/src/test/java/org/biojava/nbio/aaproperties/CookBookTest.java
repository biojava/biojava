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
package org.biojava.nbio.aaproperties;

import org.biojava.nbio.aaproperties.PeptideProperties;
import org.biojava.nbio.aaproperties.xml.AminoAcidCompositionTable;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.xml.bind.JAXBException;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.Map;

public class CookBookTest {

	// TODO this test doesn't assert anything, i.e. it is not a test! Should fix! - JD 2016-03-07

	private final static Logger logger = LoggerFactory.getLogger(CookBookTest.class);

	@Test
	public void shortExample1(){
		String sequence = "QIKDLLVSSSTDLDTTLVLVNAIYFKGMWKTAFNAEDTREMPFHVTKQESKPVQMMCMNNSFNVATLPAE";
		logger.debug("Molecular Weight: {}", PeptideProperties.getMolecularWeight(sequence));
	}

	@Test
	public void shortExample2() throws FileNotFoundException, JAXBException{
		String sequence = "QIKDLLVSSSTDLDTTLVLVNAIYFKGMWKTAFNAEDTREMPFHVTKQESKPVQMMCMNNSFNVATLPAE";
		File aminoAcidCompositionFile = new File("./src/main/resources/AminoAcidComposition.xml");
		logger.debug("Molecular Weight: {}", PeptideProperties.getMolecularWeight(sequence, aminoAcidCompositionFile));
	}

	@Test
	public void shortExample3() throws Exception{
		String[] sequences = new String[3];
		sequences[0] = "QIKDLLVSSSTDLDTTLVLVNAIYFKGMWKTAFNAEDTREMPFHVTKQESKPVQMMCMNNSFNVATLPAE";
		sequences[1] = "KMKILELPFASGDLSMLVLLPDEVSDLERIEKTINFEKLTEWTNPNTMEKRRVKVYLPQMKIEEKYNLTS";
		sequences[2] = "VLMALGMTDLFIPSANLTGISSAESLKISQAVHGAFMELSEDGIEMAGSTGVIEDIKHSPESEQFRADHP";

		File aminoAcidCompositionFile = new File("./src/main/resources/AminoAcidComposition.xml");
		AminoAcidCompositionTable table = PeptideProperties.obtainAminoAcidCompositionTable(aminoAcidCompositionFile);

		for(String sequence:sequences){
		    logger.debug("Molecular Weight: {}", PeptideProperties.getMolecularWeightBasedOnXML(sequence, table));
		}
	}

	@Test
	public void shortExample4(){
		String sequence = "QIKDLLVSSSTDLDTTLVLVNAIYFKGMWKTAFNAEDTREMPFHVTKQESKPVQMMCMNNSFNVATLPAE";

		//Enrichment of a specific amino acid type
		logger.debug("Composition of A: {}", PeptideProperties.getEnrichment(sequence, "A"));

		//Enrichment of a list of amino acid types
		Map<String, Double> composition = PeptideProperties.getAACompositionString(sequence);
		for(String aa:composition.keySet()){
			logger.debug("Composition of {}: {}", aa, composition.get(aa));
		}
	}


	@Test
	public void shortExample5(){
		String sequence = "QIKDLLVSSSTDLDTTLVLVNAIYFKGMWKTAFNAEDTRECMPFHVTKQESKPVQMMCMNNSFNVATLPAE";

		//Absorbance
		logger.debug("Absorbance (Cys Reduced): {}", PeptideProperties.getAbsorbance(sequence, true));
		logger.debug("Absorbance (Cys Not Reduced): {}", PeptideProperties.getAbsorbance(sequence, false));

		//Extinction Coefficient
		logger.debug("Extinction Coefficient (Cys Reduced): {}", PeptideProperties.getExtinctionCoefficient(sequence, true));
		logger.debug("Extinction Coefficient (Cys Not Reduced): {}", PeptideProperties.getExtinctionCoefficient(sequence, false));

		//Instability Index
		logger.debug("Instability Index: {}", PeptideProperties.getInstabilityIndex(sequence));

		//Apliphatic Index
		logger.debug("Apliphatic Index: {}", PeptideProperties.getApliphaticIndex(sequence));

		//Average Hydropathy Value
		logger.debug("Average Hydropathy Value: {}", PeptideProperties.getAvgHydropathy(sequence));

		//Isoelectric Point
		logger.debug("Isoelectric Point: {}", PeptideProperties.getIsoelectricPoint(sequence));

		//Net Charge
		logger.debug("Net Charge at pH 7: {}", PeptideProperties.getNetCharge(sequence));
	}
}

package org.biojava3.aaproperties;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Map;

import javax.xml.bind.JAXBException;

import org.biojava3.aaproperties.PeptideProperties;
import org.biojava3.aaproperties.xml.AminoAcidCompositionTable;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class CookBookTester {
	
	private final static Logger logger = LoggerFactory.getLogger(CookBookTester.class);

	@Test
	public void shortExample1(){
		String sequence = "QIKDLLVSSSTDLDTTLVLVNAIYFKGMWKTAFNAEDTREMPFHVTKQESKPVQMMCMNNSFNVATLPAE";
		logger.info("Molecular Weight: {}", PeptideProperties.getMolecularWeight(sequence));
	}
	
	@Test
	public void shortExample2() throws FileNotFoundException, JAXBException{
		String sequence = "QIKDLLVSSSTDLDTTLVLVNAIYFKGMWKTAFNAEDTREMPFHVTKQESKPVQMMCMNNSFNVATLPAE";
		File aminoAcidCompositionFile = new File("./src/main/resources/AminoAcidComposition.xml");
		logger.info("Molecular Weight: {}", PeptideProperties.getMolecularWeight(sequence, aminoAcidCompositionFile));
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
		    logger.info("Molecular Weight: {}", PeptideProperties.getMolecularWeightBasedOnXML(sequence, table));
		}
	}
	
	@Test
	public void shortExample4(){
		String sequence = "QIKDLLVSSSTDLDTTLVLVNAIYFKGMWKTAFNAEDTREMPFHVTKQESKPVQMMCMNNSFNVATLPAE";
		
		//Enrichment of a specific amino acid type
		logger.info("Composition of A: {}", PeptideProperties.getEnrichment(sequence, "A"));
		
		//Enrichment of a list of amino acid types
		Map<String, Double> composition = PeptideProperties.getAACompositionString(sequence);
		for(String aa:composition.keySet()){
			logger.info("Composition of {}: {}", aa, composition.get(aa));
		}
	}
	
	
	@Test
	public void shortExample5(){
		String sequence = "QIKDLLVSSSTDLDTTLVLVNAIYFKGMWKTAFNAEDTRECMPFHVTKQESKPVQMMCMNNSFNVATLPAE";
		
		//Absorbance
		logger.info("Absorbance (Cys Reduced): {}", PeptideProperties.getAbsorbance(sequence, true));
		logger.info("Absorbance (Cys Not Reduced): {}", PeptideProperties.getAbsorbance(sequence, false));
		
		//Extinction Coefficient
		logger.info("Extinction Coefficient (Cys Reduced): {}", PeptideProperties.getExtinctionCoefficient(sequence, true));
		logger.info("Extinction Coefficient (Cys Not Reduced): {}", PeptideProperties.getExtinctionCoefficient(sequence, false));
		
		//Instability Index
		logger.info("Instability Index: {}", PeptideProperties.getInstabilityIndex(sequence));
		
		//Apliphatic Index
		logger.info("Apliphatic Index: {}", PeptideProperties.getApliphaticIndex(sequence));
		
		//Average Hydropathy Value
		logger.info("Average Hydropathy Value: {}", PeptideProperties.getAvgHydropathy(sequence));
		
		//Isoelectric Point
		logger.info("Isoelectric Point: {}", PeptideProperties.getIsoelectricPoint(sequence));
		
		//Net Charge
		logger.info("Net Charge at pH 7: {}", PeptideProperties.getNetCharge(sequence));
	}
}

package org.biojava3.aaproperties.xml;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Map;

import javax.xml.bind.JAXBException;

import org.biojava3.aaproperties.PeptideProperties;
import org.junit.Test;

public class CookBookTester {
	private final boolean ignoreCase = true;
	
	@Test
	public void shortExample1(){
		String sequence = "QIKDLLVSSSTDLDTTLVLVNAIYFKGMWKTAFNAEDTREMPFHVTKQESKPVQMMCMNNSFNVATLPAE";
		System.out.println("Molecular Weight: " + PeptideProperties.getMolecularWeight(sequence, ignoreCase));
	}
	
	@Test
	public void shortExample2() throws FileNotFoundException, JAXBException{
		String sequence = "QIKDLLVSSSTDLDTTLVLVNAIYFKGMWKTAFNAEDTREMPFHVTKQESKPVQMMCMNNSFNVATLPAE";
		File aminoAcidCompositionFile = new File("./src/main/resources/AminoAcidComposition.xml");
		System.out.println("Molecular Weight: " + PeptideProperties.getMolecularWeight(sequence, aminoAcidCompositionFile, ignoreCase));
	}
	
	@Test
	public void shortExample3() throws Exception{
		String[] sequences = new String[3];
		sequences[0] = "QIKDLLVSSSTDLDTTLVLVNAIYFKGMWKTAFNAEDTREMPFHVTKQESKPVQMMCMNNSFNVATLPAE";
		sequences[1] = "KMKILELPFASGDLSMLVLLPDEVSDLERIEKTINFEKLTEWTNPNTMEKRRVKVYLPQMKIEEKYNLTS";
		sequences[2] = "VLMALGMTDLFIPSANLTGISSAESLKISQAVHGAFMELSEDGIEMAGSTGVIEDIKHSPESEQFRADHP";
		 
		File aminoAcidCompositionFile = new File("./src/main/resources/AminoAcidComposition.xml");
		AminoAcidCompositionTable table = PeptideProperties.obtainAminoAcidCompositionTable(aminoAcidCompositionFile, ignoreCase);
		 
		for(String sequence:sequences){
		    System.out.println("Molecular Weight: " + PeptideProperties.getMolecularWeightBasedOnXML(sequence, table, ignoreCase));
		}
	}
	
	@Test
	public void shortExample4(){
		String sequence = "QIKDLLVSSSTDLDTTLVLVNAIYFKGMWKTAFNAEDTREMPFHVTKQESKPVQMMCMNNSFNVATLPAE";
		
		//Enrichment of a specific amino acid type
		System.out.println("Composition of A: " + PeptideProperties.getEnrichment(sequence, "A", ignoreCase));
		
		//Enrichment of a list of amino acid types
		Map<String, Double> composition = PeptideProperties.getAACompositionString(sequence, ignoreCase);
		for(String aa:composition.keySet()){
			System.out.println("Composition of " + aa + ": " + composition.get(aa));
		}
	}
	
	
	@Test
	public void shortExample5(){
		String sequence = "QIKDLLVSSSTDLDTTLVLVNAIYFKGMWKTAFNAEDTRECMPFHVTKQESKPVQMMCMNNSFNVATLPAE";
		
		//Absorbance
		System.out.println("Absorbance (Cys Reduced): " + PeptideProperties.getAbsorbance(sequence, true, ignoreCase));
		System.out.println("Absorbance (Cys Not Reduced): " + PeptideProperties.getAbsorbance(sequence, false, ignoreCase));
		
		//Extinction Coefficient
		System.out.println("Extinction Coefficient (Cys Reduced): " + PeptideProperties.getExtinctionCoefficient(sequence, true, ignoreCase));
		System.out.println("Extinction Coefficient (Cys Not Reduced): " + PeptideProperties.getExtinctionCoefficient(sequence, false, ignoreCase));
		
		//Instability Index
		System.out.println("Instability Index: " + PeptideProperties.getInstabilityIndex(sequence, ignoreCase));
		
		//Apliphatic Index
		System.out.println("Apliphatic Index: " + PeptideProperties.getApliphaticIndex(sequence, ignoreCase));
		
		//Average Hydropathy Value
		System.out.println("Average Hydropathy Value: " + PeptideProperties.getAvgHydropathy(sequence, ignoreCase));
		
		//Isoelectric Point
		System.out.println("Isoelectric Point: " + PeptideProperties.getIsoelectricPoint(sequence, ignoreCase));
		
		//Net Charge
		System.out.println("Net Charge at pH 7: " + PeptideProperties.getNetCharge(sequence, ignoreCase));
	}
}

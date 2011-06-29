package org.biojava3.aaproperties.xml;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Map;

import javax.xml.bind.JAXBException;

import org.biojava3.aaproperties.PeptideProperties;
import org.junit.Test;

public class CookBookTester {
	@Test
	public void shortExample1(){
		String sequence = "QIKDLLVSSSTDLDTTLVLVNAIYFKGMWKTAFNAEDTREMPFHVTKQESKPVQMMCMNNSFNVATLPAE";
		System.out.println("Molecular Weight: " + PeptideProperties.getMolecularWeight(sequence));
	}
	
	@Test
	public void shortExample2() throws FileNotFoundException, JAXBException{
		String sequence = "QIKDLLVSSSTDLDTTLVLVNAIYFKGMWKTAFNAEDTREMPFHVTKQESKPVQMMCMNNSFNVATLPAE";
		File elementMassFile = new File("./src/main/resources/ElementMass.xml");
		File aminoAcidCompositionFile = new File("./src/main/resources/AminoAcidComposition.xml");
		System.out.println("Molecular Weight: " + PeptideProperties.getMolecularWeight(sequence, elementMassFile, aminoAcidCompositionFile));
	}
	
	@Test
	public void shortExample3() throws Exception{
		String[] sequences = new String[3];
		sequences[0] = "QIKDLLVSSSTDLDTTLVLVNAIYFKGMWKTAFNAEDTREMPFHVTKQESKPVQMMCMNNSFNVATLPAE";
		sequences[1] = "KMKILELPFASGDLSMLVLLPDEVSDLERIEKTINFEKLTEWTNPNTMEKRRVKVYLPQMKIEEKYNLTS";
		sequences[2] = "VLMALGMTDLFIPSANLTGISSAESLKISQAVHGAFMELSEDGIEMAGSTGVIEDIKHSPESEQFRADHP";
		 
		File elementMassFile = new File("./src/main/resources/ElementMass.xml");
		File aminoAcidCompositionFile = new File("./src/main/resources/AminoAcidComposition.xml");
		AminoAcidCompositionTable table = PeptideProperties.obtainAminoAcidCompositionTable(elementMassFile, aminoAcidCompositionFile);
		 
		for(String sequence:sequences){
		    System.out.println("Molecular Weight: " + PeptideProperties.getMolecularWeightBasedOnXML(sequence, table));
		}
	}
	
	@Test
	public void shortExample4(){
		String sequence = "QIKDLLVSSSTDLDTTLVLVNAIYFKGMWKTAFNAEDTREMPFHVTKQESKPVQMMCMNNSFNVATLPAE";
		System.out.println("Composition of A: " + PeptideProperties.getEnrichment(sequence, "A"));
	}
	
	@Test
	public void shortExample5(){
		String sequence = "QIKDLLVSSSTDLDTTLVLVNAIYFKGMWKTAFNAEDTREMPFHVTKQESKPVQMMCMNNSFNVATLPAE";
		Map<String, Double> composition = PeptideProperties.getAACompositionString(sequence);
		for(String aa:composition.keySet()){
			System.out.println("Composition of " + aa + ": " + composition.get(aa));
		}
	}
	
	@Test
	public void shortExample6(){
		String sequence = "QIKDLLVSSSTDLDTTLVLVNAIYFKGMWKTAFNAEDTRECMPFHVTKQESKPVQMMCMNNSFNVATLPAE";
		System.out.println("Absorbance (Cys Reduced): " + PeptideProperties.getAbsorbance(sequence, true));
		System.out.println("Absorbance (Cys Not Reduced): " + PeptideProperties.getAbsorbance(sequence, false));
		System.out.println("Extinction Coefficient (Cys Reduced): " + PeptideProperties.getExtinctionCoefficient(sequence, true));
		System.out.println("Extinction Coefficient (Cys Not Reduced): " + PeptideProperties.getExtinctionCoefficient(sequence, false));
		System.out.println("Instability Index: " + PeptideProperties.getInstabilityIndex(sequence));
		System.out.println("Apliphatic Index: " + PeptideProperties.getApliphaticIndex(sequence));
		System.out.println("Average Hydropathy Value: " + PeptideProperties.getAvgHydropathy(sequence));
		System.out.println("Isoelectric Point: " + PeptideProperties.getIsoelectricPoint(sequence));
		System.out.println("Net Charge at pH 7: " + PeptideProperties.getNetCharge(sequence));
	}
}

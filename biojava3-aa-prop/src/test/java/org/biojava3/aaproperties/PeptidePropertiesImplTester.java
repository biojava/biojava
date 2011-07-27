package org.biojava3.aaproperties;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Map;

import javax.xml.bind.JAXBException;

import org.biojava3.aaproperties.xml.AminoAcidCompositionTable;
import org.junit.Test;
import static junit.framework.Assert.*;

public class PeptidePropertiesImplTester {

	/**
	 * Test input 
	 */
	private final String sequence = "MTADGPCRELLCQLRAAVRHRWWC1";
	private final String fullInvalidSequence = "3176412372301230183--2310";
	private final boolean ignoreCase = true;
	
	@Test
	public void testAAComposition() { 
		//'W', 'C', 'M', 'H', 'Y', 'F', 'Q', 'N', 'I', 'R', 'D', 'P', 'T', 'K', 'E', 'V', 'S', 'G', 'A', 'L'
		Map<String, Double> composition = PeptideProperties.getAACompositionString(sequence, ignoreCase);
		int sequenceLength = sequence.length() - Utils.getNumberOfInvalidChar(sequence, null, ignoreCase);
		assertEquals(2.0/sequenceLength, composition.get("W"));
		assertEquals(3.0/sequenceLength, composition.get("C"));
		assertEquals(1.0/sequenceLength, composition.get("M"));
		assertEquals(1.0/sequenceLength, composition.get("H"));
		assertEquals(0.0/sequenceLength, composition.get("Y"));
		assertEquals(0.0/sequenceLength, composition.get("F"));
		assertEquals(1.0/sequenceLength, composition.get("Q"));
		assertEquals(0.0/sequenceLength, composition.get("N"));
		assertEquals(0.0/sequenceLength, composition.get("I"));
		assertEquals(4.0/sequenceLength, composition.get("R"));
		assertEquals(1.0/sequenceLength, composition.get("D"));
		assertEquals(1.0/sequenceLength, composition.get("P"));
		assertEquals(1.0/sequenceLength, composition.get("T"));
		assertEquals(0.0/sequenceLength, composition.get("K"));
		assertEquals(1.0/sequenceLength, composition.get("E"));
		assertEquals(1.0/sequenceLength, composition.get("V"));
		assertEquals(0.0/sequenceLength, composition.get("S"));
		assertEquals(1.0/sequenceLength, composition.get("G"));
		assertEquals(3.0/sequenceLength, composition.get("A"));
		assertEquals(3.0/sequenceLength, composition.get("L"));
		
		Map<String, Double> iComposition = PeptideProperties.getAACompositionString(fullInvalidSequence, ignoreCase);
		assertEquals(0.0, iComposition.get("W"));
		assertEquals(0.0, iComposition.get("C"));
		assertEquals(0.0, iComposition.get("M"));
		assertEquals(0.0, iComposition.get("H"));
		assertEquals(0.0, iComposition.get("Y"));
		assertEquals(0.0, iComposition.get("F"));
		assertEquals(0.0, iComposition.get("Q"));
		assertEquals(0.0, iComposition.get("N"));
		assertEquals(0.0, iComposition.get("I"));
		assertEquals(0.0, iComposition.get("R"));
		assertEquals(0.0, iComposition.get("D"));
		assertEquals(0.0, iComposition.get("P"));
		assertEquals(0.0, iComposition.get("T"));
		assertEquals(0.0, iComposition.get("K"));
		assertEquals(0.0, iComposition.get("E"));
		assertEquals(0.0, iComposition.get("V"));
		assertEquals(0.0, iComposition.get("S"));
		assertEquals(0.0, iComposition.get("G"));
		assertEquals(0.0, iComposition.get("A"));
		assertEquals(0.0, iComposition.get("L"));

		//Null would be returned for invalid character
		assertNotSame(0d, composition.get("Z"));
		assertNull(composition.get(null));
		assertNull(composition.get(""));
		assertNull(composition.get("1"));
	}
	
	@Test()
	public void testEnrichment() {
		//'W', 'C', 'M', 'H', 'Y', 'F', 'Q', 'N', 'I', 'R', 'D', 'P', 'T', 'K', 'E', 'V', 'S', 'G', 'A', 'L'
		int sequenceLength = sequence.length();
		assertEquals(2.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "W", ignoreCase));
		assertEquals(3.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "C", ignoreCase));
		assertEquals(1.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "M", ignoreCase));
		assertEquals(1.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "H", ignoreCase));
		assertEquals(0.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "Y", ignoreCase));
		assertEquals(0.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "F", ignoreCase));
		assertEquals(1.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "Q", ignoreCase));
		assertEquals(0.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "N", ignoreCase));
		assertEquals(0.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "I", ignoreCase));
		assertEquals(4.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "R", ignoreCase));
		assertEquals(1.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "D", ignoreCase));
		assertEquals(1.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "P", ignoreCase));
		assertEquals(1.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "T", ignoreCase));
		assertEquals(0.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "K", ignoreCase));
		assertEquals(1.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "E", ignoreCase));
		assertEquals(1.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "V", ignoreCase));
		assertEquals(0.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "S", ignoreCase));
		assertEquals(1.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "G", ignoreCase));
		assertEquals(3.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "A", ignoreCase));
		assertEquals(3.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "L", ignoreCase));
		
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "W", ignoreCase));
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "C", ignoreCase));
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "M", ignoreCase));
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "H", ignoreCase));
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "Y", ignoreCase));
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "F", ignoreCase));
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "Q", ignoreCase));
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "N", ignoreCase));
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "I", ignoreCase));
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "R", ignoreCase));
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "D", ignoreCase));
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "P", ignoreCase));
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "T", ignoreCase));
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "K", ignoreCase));
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "E", ignoreCase));
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "V", ignoreCase));
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "S", ignoreCase));
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "G", ignoreCase));
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "A", ignoreCase));
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "L", ignoreCase));
		assertEquals(0.0, PeptideProperties.getEnrichment(sequence, "X", ignoreCase));
	}
	
	@Test (expected = NullPointerException.class)
	public void testEnrichmentNull(){
		assertNull(PeptideProperties.getEnrichment(sequence, "1", ignoreCase));
		assertEquals(0.0, PeptideProperties.getEnrichment(sequence, "", ignoreCase));
	}
	
	@Test
	public void testMolecularWeight(){
		//http://www.innovagen.se/custom-peptide-synthesis/peptide-property-calculator/peptide-property-calculator.asp
		//http://au.expasy.org/cgi-bin/protparam
		//2872.4 is the value computed by the above two web tools
		assertEquals(2872.4, Utils.roundToDecimals(PeptideProperties.getMolecularWeight(sequence, ignoreCase), 1));
		assertEquals(0.0, PeptideProperties.getMolecularWeight("Z", ignoreCase));
		assertEquals(0.0, PeptideProperties.getMolecularWeight("1", ignoreCase));
		
		assertEquals(0.0, PeptideProperties.getMolecularWeight(fullInvalidSequence, ignoreCase));
	}
	
	@Test
	public void testMolecularWeightXML() throws FileNotFoundException, JAXBException{
		File elementMassFile = new File("./src/main/resources/ElementMass.xml");
		File aminoAcidCompositionFile = new File("./src/main/resources/AminoAcidComposition.xml");
		
		assertEquals(
			Utils.roundToDecimals(PeptideProperties.getMolecularWeight("A", elementMassFile, aminoAcidCompositionFile, ignoreCase) 
					* 5.0 -  4 * (17.0073 + 1.0079), 5), 
			Utils.roundToDecimals(PeptideProperties.getMolecularWeight("AAAAA", elementMassFile, aminoAcidCompositionFile, ignoreCase), 5));
	}
	
	@Test
	public void testMolecularWeightXMLSingleFile() throws FileNotFoundException, JAXBException{
		File aminoAcidCompositionFile = new File("./src/main/resources/AminoAcidComposition.xml");
		
		assertEquals(
			Utils.roundToDecimals(PeptideProperties.getMolecularWeight("A", aminoAcidCompositionFile, ignoreCase) * 5.0 -  4 * (17.0073 + 1.0079), 5), 
			Utils.roundToDecimals(PeptideProperties.getMolecularWeight("AAAAA", aminoAcidCompositionFile, ignoreCase), 5));
	}
	
	@Test
	public void testMolecularWeightBasedOnAminoAcidCompositionTable() throws Exception{
		File elementMassFile = new File("./src/main/resources/ElementMass.xml");
		File aminoAcidCompositionFile = new File("./src/main/resources/AminoAcidComposition.xml");
		AminoAcidCompositionTable table = PeptideProperties.obtainAminoAcidCompositionTable(elementMassFile, aminoAcidCompositionFile, ignoreCase);
		
		assertEquals(
			Utils.roundToDecimals(PeptideProperties.getMolecularWeightBasedOnXML("A", table, ignoreCase) * 5.0 -  4 * (17.0073 + 1.0079), 5), 
			Utils.roundToDecimals(PeptideProperties.getMolecularWeightBasedOnXML("AAAAA", table, ignoreCase), 5));
	}
	
	@Test (expected = NullPointerException.class)
	public void testMolecularWeightXMLNull() throws FileNotFoundException, JAXBException{
		PeptideProperties.getMolecularWeight(sequence, null, null, ignoreCase);
	}
	
	@Test
	public void testExtinctionCoefficient(){
		//http://au.expasy.org/cgi-bin/protparam
		assertEquals(11125.0, PeptideProperties.getExtinctionCoefficient(sequence, false, ignoreCase));
		assertEquals(11000.0, PeptideProperties.getExtinctionCoefficient(sequence, true, ignoreCase));
		
		assertEquals(0.0, PeptideProperties.getExtinctionCoefficient(fullInvalidSequence, true, ignoreCase));
		assertEquals(0.0, PeptideProperties.getExtinctionCoefficient(fullInvalidSequence, false, ignoreCase));
	}
	
	@Test (expected = NullPointerException.class)
	public void testExtinctionCoefficientNull(){
		assertEquals(11000.0, PeptideProperties.getExtinctionCoefficient(null, true, ignoreCase));
	}
	
	@Test
	public void testAbsorbance(){
		//http://au.expasy.org/cgi-bin/protparam
		assertEquals(3.830, Utils.roundToDecimals(PeptideProperties.getAbsorbance(sequence, true, ignoreCase), 3));
		assertEquals(3.873, Utils.roundToDecimals(PeptideProperties.getAbsorbance(sequence, false, ignoreCase), 3));
		
		assertEquals(0.0, Utils.roundToDecimals(PeptideProperties.getAbsorbance(fullInvalidSequence, true, ignoreCase), 3));
		assertEquals(0.0, Utils.roundToDecimals(PeptideProperties.getAbsorbance(fullInvalidSequence, false, ignoreCase), 3));
	}
	
	@Test (expected = NullPointerException.class)
	public void testAbsorbanceNull(){
		assertEquals(3.830, PeptideProperties.getAbsorbance(null, false, ignoreCase));
	}
	
	
	@Test
	public void testInstabilityIndex(){
		//http://au.expasy.org/cgi-bin/protparam
		assertEquals(38.48, Utils.roundToDecimals(PeptideProperties.getInstabilityIndex(sequence, ignoreCase), 2));
		assertEquals(0.0, Utils.roundToDecimals(PeptideProperties.getInstabilityIndex(fullInvalidSequence, ignoreCase), 2));
	}
	
	@Test (expected = NullPointerException.class)
	public void testInstabilityIndexNull(){
		assertEquals(38.48, PeptideProperties.getInstabilityIndex(null, ignoreCase));
	}
	
	@Test
	public void testApliphaticIndex(){
		//http://au.expasy.org/cgi-bin/protparam
		assertEquals(73.33, Utils.roundToDecimals(PeptideProperties.getApliphaticIndex(sequence, ignoreCase), 2));
		assertEquals(0.0, Utils.roundToDecimals(PeptideProperties.getApliphaticIndex(fullInvalidSequence, ignoreCase), 2));
	}
	
	@Test (expected = NullPointerException.class)
	public void testApliphaticIndexNull(){
		assertEquals(73.33, Utils.roundToDecimals(PeptideProperties.getApliphaticIndex(null, ignoreCase), 2));
	}
	
	@Test
	public void testAverageHydropathy(){
		//http://au.expasy.org/cgi-bin/protparam
		assertEquals(-0.242, Utils.roundToDecimals(PeptideProperties.getAvgHydropathy(sequence, ignoreCase), 3));
		assertEquals(0.0, Utils.roundToDecimals(PeptideProperties.getAvgHydropathy(fullInvalidSequence, ignoreCase), 3));
	}
	
	@Test (expected = NullPointerException.class)
	public void testAverageHydropathyNull(){
		assertEquals(-0.242, Utils.roundToDecimals(PeptideProperties.getAvgHydropathy(null, ignoreCase), 3));
	}
	
	@Test
	public void testIsoelectricPoint(){
		//http://www.innovagen.se/custom-peptide-synthesis/peptide-property-calculator/peptide-property-calculator.asp
		assertEquals(8.6, Utils.roundToDecimals(PeptideProperties.getIsoelectricPoint(sequence, ignoreCase), 1));
		assertEquals(7.0, Utils.roundToDecimals(PeptideProperties.getIsoelectricPoint(fullInvalidSequence, ignoreCase), 1));
	}
	
	@Test (expected = NullPointerException.class)
	public void testIsoelectricPointNull(){
		assertEquals(8.6, Utils.roundToDecimals(PeptideProperties.getIsoelectricPoint(null, ignoreCase), 1));
	}
	
	@Test
	public void testNetCharge(){
		//http://www.innovagen.se/custom-peptide-synthesis/peptide-property-calculator/peptide-property-calculator.asp
		assertEquals(2.0, Utils.roundToDecimals(PeptideProperties.getNetCharge(sequence, ignoreCase), 1));
		assertEquals(0.0, Utils.roundToDecimals(PeptideProperties.getNetCharge(fullInvalidSequence, ignoreCase), 1));
	}
	
	@Test (expected = NullPointerException.class)
	public void testNetChargeNull(){
		assertEquals(8.6, PeptideProperties.getNetCharge(null, ignoreCase));
	}
}

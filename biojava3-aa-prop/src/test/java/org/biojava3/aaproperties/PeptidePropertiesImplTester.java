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

	@Test
	public void testAAComposition() { 
		//'W', 'C', 'M', 'H', 'Y', 'F', 'Q', 'N', 'I', 'R', 'D', 'P', 'T', 'K', 'E', 'V', 'S', 'G', 'A', 'L'
		Map<String, Double> composition = PeptideProperties.getAACompositionString(sequence);
		int sequenceLength = sequence.length() - Utils.getNumberOfInvalidChar(sequence, null, true);
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

		Map<String, Double> iComposition = PeptideProperties.getAACompositionString(fullInvalidSequence);
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
		assertEquals(2.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "W"));
		assertEquals(3.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "C"));
		assertEquals(1.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "M"));
		assertEquals(1.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "H"));
		assertEquals(0.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "Y"));
		assertEquals(0.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "F"));
		assertEquals(1.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "Q"));
		assertEquals(0.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "N"));
		assertEquals(0.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "I"));
		assertEquals(4.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "R"));
		assertEquals(1.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "D"));
		assertEquals(1.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "P"));
		assertEquals(1.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "T"));
		assertEquals(0.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "K"));
		assertEquals(1.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "E"));
		assertEquals(1.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "V"));
		assertEquals(0.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "S"));
		assertEquals(1.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "G"));
		assertEquals(3.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "A"));
		assertEquals(3.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "L"));

		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "W"));
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "C"));
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "M"));
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "H"));
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "Y"));
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "F"));
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "Q"));
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "N"));
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "I"));
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "R"));
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "D"));
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "P"));
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "T"));
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "K"));
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "E"));
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "V"));
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "S"));
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "G"));
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "A"));
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "L"));
		assertEquals(0.0, PeptideProperties.getEnrichment(sequence, "X"));
	}

	@Test (expected = NullPointerException.class)
	public void testEnrichmentNull(){
		assertNull(PeptideProperties.getEnrichment(sequence, "1"));
		assertEquals(0.0, PeptideProperties.getEnrichment(sequence, ""));
	}

	@Test
	public void testMolecularWeight(){
		//http://www.innovagen.se/custom-peptide-synthesis/peptide-property-calculator/peptide-property-calculator.asp
		//http://au.expasy.org/cgi-bin/protparam
		//2872.4 is the value computed by the above two web tools
		assertEquals(2872.4, Utils.roundToDecimals(PeptideProperties.getMolecularWeight(sequence), 1));
		assertEquals(0.0, PeptideProperties.getMolecularWeight("Z"));
		assertEquals(0.0, PeptideProperties.getMolecularWeight("1"));

		assertEquals(0.0, PeptideProperties.getMolecularWeight(fullInvalidSequence));
	}

	@Test
	public void testMolecularWeightXML() throws FileNotFoundException, JAXBException{
		File elementMassFile = new File("./src/main/resources/ElementMass.xml");
		File aminoAcidCompositionFile = new File("./src/main/resources/AminoAcidComposition.xml");

		assertEquals(
				Utils.roundToDecimals(PeptideProperties.getMolecularWeight("A", elementMassFile, aminoAcidCompositionFile) 
						* 5.0 -  4 * (17.0073 + 1.0079), 5), 
						Utils.roundToDecimals(PeptideProperties.getMolecularWeight("AAAAA", elementMassFile, aminoAcidCompositionFile), 5));
	}

	@Test
	public void testMolecularWeightXMLSingleFile() throws FileNotFoundException, JAXBException{
		File aminoAcidCompositionFile = new File("./src/main/resources/AminoAcidComposition.xml");

		assertEquals(
				Utils.roundToDecimals(PeptideProperties.getMolecularWeight("A", aminoAcidCompositionFile) * 5.0 -  4 * (17.0073 + 1.0079), 5), 
				Utils.roundToDecimals(PeptideProperties.getMolecularWeight("AAAAA", aminoAcidCompositionFile), 5));
	}

	@Test
	public void testMolecularWeightBasedOnAminoAcidCompositionTable() throws Exception{
		File elementMassFile = new File("./src/main/resources/ElementMass.xml");
		File aminoAcidCompositionFile = new File("./src/main/resources/AminoAcidComposition.xml");
		AminoAcidCompositionTable table = PeptideProperties.obtainAminoAcidCompositionTable(elementMassFile, aminoAcidCompositionFile);

		assertEquals(
				Utils.roundToDecimals(PeptideProperties.getMolecularWeightBasedOnXML("A", table) * 5.0 -  4 * (17.0073 + 1.0079), 5), 
				Utils.roundToDecimals(PeptideProperties.getMolecularWeightBasedOnXML("AAAAA", table), 5));
	}

	@Test (expected = NullPointerException.class)
	public void testMolecularWeightXMLNull() throws FileNotFoundException, JAXBException{
		PeptideProperties.getMolecularWeight(sequence, null, null);
	}

	@Test
	public void testExtinctionCoefficient(){
		//http://au.expasy.org/cgi-bin/protparam
		assertEquals(11125.0, PeptideProperties.getExtinctionCoefficient(sequence, false));
		assertEquals(11000.0, PeptideProperties.getExtinctionCoefficient(sequence, true));

		assertEquals(0.0, PeptideProperties.getExtinctionCoefficient(fullInvalidSequence, true));
		assertEquals(0.0, PeptideProperties.getExtinctionCoefficient(fullInvalidSequence, false));
	}

	@Test (expected = NullPointerException.class)
	public void testExtinctionCoefficientNull(){
		assertEquals(11000.0, PeptideProperties.getExtinctionCoefficient(null, true));
	}

	@Test
	public void testAbsorbance(){
		//http://au.expasy.org/cgi-bin/protparam
		assertEquals(3.830, Utils.roundToDecimals(PeptideProperties.getAbsorbance(sequence, true), 3));
		assertEquals(3.873, Utils.roundToDecimals(PeptideProperties.getAbsorbance(sequence, false), 3));

		assertEquals(0.0, Utils.roundToDecimals(PeptideProperties.getAbsorbance(fullInvalidSequence, true), 3));
		assertEquals(0.0, Utils.roundToDecimals(PeptideProperties.getAbsorbance(fullInvalidSequence, false), 3));
	}

	@Test (expected = NullPointerException.class)
	public void testAbsorbanceNull(){
		assertEquals(3.830, PeptideProperties.getAbsorbance(null, false));
	}


	@Test
	public void testInstabilityIndex(){
		//http://au.expasy.org/cgi-bin/protparam
		assertEquals(38.48, Utils.roundToDecimals(PeptideProperties.getInstabilityIndex(sequence), 2));
		assertEquals(0.0, Utils.roundToDecimals(PeptideProperties.getInstabilityIndex(fullInvalidSequence), 2));
	}

	@Test (expected = NullPointerException.class)
	public void testInstabilityIndexNull(){
		assertEquals(38.48, PeptideProperties.getInstabilityIndex(null));
	}

	@Test
	public void testApliphaticIndex(){
		//http://au.expasy.org/cgi-bin/protparam
		assertEquals(73.33, Utils.roundToDecimals(PeptideProperties.getApliphaticIndex(sequence), 2));
		assertEquals(0.0, Utils.roundToDecimals(PeptideProperties.getApliphaticIndex(fullInvalidSequence), 2));
	}

	@Test (expected = NullPointerException.class)
	public void testApliphaticIndexNull(){
		assertEquals(73.33, Utils.roundToDecimals(PeptideProperties.getApliphaticIndex(null), 2));
	}

	@Test
	public void testAverageHydropathy(){
		//http://au.expasy.org/cgi-bin/protparam
		assertEquals(-0.242, Utils.roundToDecimals(PeptideProperties.getAvgHydropathy(sequence), 3));
		assertEquals(0.0, Utils.roundToDecimals(PeptideProperties.getAvgHydropathy(fullInvalidSequence), 3));
	}

	@Test (expected = NullPointerException.class)
	public void testAverageHydropathyNull(){
		assertEquals(-0.242, Utils.roundToDecimals(PeptideProperties.getAvgHydropathy(null), 3));
	}

	@Test
	public void testIsoelectricPointInnovagen(){
		/*
		 * Test for Innovagen
		 */
		//http://www.innovagen.se/custom-peptide-synthesis/peptide-property-calculator/peptide-property-calculator.asp
		assertEquals(9.01, Utils.roundToDecimals(PeptideProperties.getIsoelectricPoint(sequence, false), 2));
		assertEquals(7.00, Utils.roundToDecimals(PeptideProperties.getIsoelectricPoint(fullInvalidSequence, false), 2));

		assertEquals(2.70, Utils.roundToDecimals(PeptideProperties.getIsoelectricPoint("ACCACAAADADADACA", false), 2));
	}

	@Test
	public void testIsoelectricPointExpasy(){
		/*
		 * Test for Expasy
		 */
		assertEquals(3.42, Utils.roundToDecimals(PeptideProperties.getIsoelectricPoint("ACCACAAADADADACA"), 2));
		assertEquals(3.42, Utils.roundToDecimals(PeptideProperties.getIsoelectricPoint("ACCACAAADADADACM"), 2));
		//		
		assertEquals(3.37, Utils.roundToDecimals(PeptideProperties.getIsoelectricPoint("ECCACAAADADADACS", true), 2));


		assertEquals(3.24, Utils.roundToDecimals(PeptideProperties.getIsoelectricPoint("ADCCACAAADADADACDAAAAAAAAAAAA", true), 2));

		//3.32 at Expasy
		assertEquals(3.32, Utils.roundToDecimals(PeptideProperties.getIsoelectricPoint("DCCACAAADADADACS", true), 2));

		assertEquals(3.17, Utils.roundToDecimals(PeptideProperties.getIsoelectricPoint("DCCACAAADADADACDAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAD", true), 2));
		assertEquals(3.37, Utils.roundToDecimals(PeptideProperties.getIsoelectricPoint("ACCACAAADADADACE", true), 2));
		assertEquals(3.32, Utils.roundToDecimals(PeptideProperties.getIsoelectricPoint("ACCACAAADADADACAAAAAAAAAAAAAAD", true), 2));
		assertEquals(3.28, Utils.roundToDecimals(PeptideProperties.getIsoelectricPoint("DCCACAAADADADACE", true), 2));

		assertEquals(8.71, Utils.roundToDecimals(PeptideProperties.getIsoelectricPoint("MTADGPCRELLCQLRAAVRHRWWC", true), 2));
		assertEquals(8.71, Utils.roundToDecimals(PeptideProperties.getIsoelectricPoint(sequence, true), 2));
		assertEquals(0.0, Utils.roundToDecimals(PeptideProperties.getIsoelectricPoint(fullInvalidSequence, true), 1));
	}

	@Test (expected = NullPointerException.class)
	public void testIsoelectricPointNull(){
		assertEquals(8.6, Utils.roundToDecimals(PeptideProperties.getIsoelectricPoint(null), 1));
	}

	@Test
	public void testNetCharge(){
		/*
		 * Test for Innovagen
		 */
		//http://www.innovagen.se/custom-peptide-synthesis/peptide-property-calculator/peptide-property-calculator.asp
		assertEquals(2.9, Utils.roundToDecimals(PeptideProperties.getNetCharge(sequence, false), 1));
		assertEquals(0.0, Utils.roundToDecimals(PeptideProperties.getNetCharge(fullInvalidSequence, false), 1));

		assertEquals(-3.2, Utils.roundToDecimals(PeptideProperties.getNetCharge("ACCACAAADADADACA", false), 1));
		/*
		 * Did not test for Expasy because in their website, net charge is not given.
		 * However, since Isoelectric point is given which rely on getNetCharge values therefore, 
		 * 	we infer that if getIsoelectricPoint is correct, getNetCharge would be correct for Expasy.
		 */
		
		
		/*
		 * Provided by Steve Darnell to compare the difference between Innovagen and Expasy
		 */
		String[] alpha = {"A",/*"B",*/"C","D","E","F","G","H","I",/*"J",*/
				"K","L","M","N",/*"O",*/"P","Q","R","S","T",
				/*"U",*/"V","W",/*"X",*/"Y"/*,"Z"*/};
		for (String aa : alpha) {
			String p = String.format("AA%sAA", aa);
			System.out.println(p);
			System.out.println("pH\tInnovagen\tExpasy\tdiff");
			for ( int i = 1; i < 15; i++) {
				double phPoint = (new Double(i)).doubleValue();
				double chrgInnovagen = PeptideProperties.getNetCharge(p,false,phPoint);
				double chrgExpasy = PeptideProperties.getNetCharge(p,true,phPoint);
				System.out.println(String.format("%2.1f\t%2.2f\t%2.2f\t%2.2f", phPoint, chrgInnovagen, 
						chrgExpasy, chrgInnovagen - chrgExpasy));
			}
		}
	}

	@Test (expected = NullPointerException.class)
	public void testNetChargeNull(){
		assertEquals(8.6, PeptideProperties.getNetCharge(null));
	}
}

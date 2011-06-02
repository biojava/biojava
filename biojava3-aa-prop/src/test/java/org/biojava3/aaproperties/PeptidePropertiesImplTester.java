package org.biojava3.aaproperties;

import java.util.Map;

import org.junit.Test;
import static junit.framework.Assert.*;

public class PeptidePropertiesImplTester {

	/**
	 * Test input 
	 */
	private final String sequence = "MTADGPCRELLCQLRAAVRHRWWC";
	
	@Test
	public void testAAComposition() { 
		//'W', 'C', 'M', 'H', 'Y', 'F', 'Q', 'N', 'I', 'R', 'D', 'P', 'T', 'K', 'E', 'V', 'S', 'G', 'A', 'L'
		Map<String, Double> composition = PeptideProperties.getAACompositionString(sequence);
		int sequenceLength = sequence.length();
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
	}
	
	@Test
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
	}
	
	@Test
	public void testMolecularWeight(){
		//http://www.innovagen.se/custom-peptide-synthesis/peptide-property-calculator/peptide-property-calculator.asp
		//http://au.expasy.org/cgi-bin/protparam
		//2872.4 is the value computed by the above two web tools
		assertEquals(2872.4, PeptideProperties.getMolecularWeight(sequence));
	}
	
	@Test
	public void testExtinctionCoefficient(){
		//http://au.expasy.org/cgi-bin/protparam
		assertEquals(11125.0, PeptideProperties.getExtinctionCoefficient(sequence, true));
		assertEquals(11000.0, PeptideProperties.getExtinctionCoefficient(sequence, false));
	}
	
	@Test
	public void testAbsorbance(){
		//http://au.expasy.org/cgi-bin/protparam
		assertEquals(3.873, PeptideProperties.getAbsorbance(sequence, true));
		assertEquals(3.830, PeptideProperties.getAbsorbance(sequence, false));
	}
	
	@Test
	public void testInstabilityIndex(){
		//http://au.expasy.org/cgi-bin/protparam
		assertEquals(38.48, PeptideProperties.getInstabilityIndex(sequence));
	}
	
	@Test
	public void testApliphaticIndex(){
		//http://au.expasy.org/cgi-bin/protparam
		assertEquals(73.33, PeptideProperties.getApliphaticIndex(sequence));
	}
	
	@Test
	public void testAverageHydropathy(){
		//http://au.expasy.org/cgi-bin/protparam
		assertEquals(-0.242, PeptideProperties.getAvgHydropathy(sequence));
	}
	
	@Test
	public void testIsoelectricPoint(){
		//http://www.innovagen.se/custom-peptide-synthesis/peptide-property-calculator/peptide-property-calculator.asp
		assertEquals(8.6, PeptideProperties.getIsoelectricPoint(sequence));
	}
	
	@Test
	public void testNetCharge(){
		//http://www.innovagen.se/custom-peptide-synthesis/peptide-property-calculator/peptide-property-calculator.asp
		assertEquals(2.0, PeptideProperties.getNetCharge(sequence));
	}
}

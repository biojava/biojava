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
import org.biojava.nbio.aaproperties.Utils;
import org.biojava.nbio.aaproperties.xml.AminoAcidCompositionTable;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.xml.bind.JAXBException;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.Map;

import static org.junit.Assert.*;

public class PeptidePropertiesImplTest {

	private final static Logger logger = LoggerFactory.getLogger(PeptidePropertiesImplTest.class);

	private static final double delta = 0.00001;

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
		assertEquals(2.0/sequenceLength,composition.get("W"), delta);
		assertEquals(3.0/sequenceLength,composition.get("C"), delta);
		assertEquals(1.0/sequenceLength,composition.get("M"), delta);
		assertEquals(1.0/sequenceLength,composition.get("H"), delta);
		assertEquals(0.0/sequenceLength,composition.get("Y"), delta);
		assertEquals(0.0/sequenceLength,composition.get("F"), delta);
		assertEquals(1.0/sequenceLength,composition.get("Q"), delta);
		assertEquals(0.0/sequenceLength,composition.get("N"), delta);
		assertEquals(0.0/sequenceLength,composition.get("I"), delta);
		assertEquals(4.0/sequenceLength,composition.get("R"), delta);
		assertEquals(1.0/sequenceLength,composition.get("D"), delta);
		assertEquals(1.0/sequenceLength,composition.get("P"), delta);
		assertEquals(1.0/sequenceLength,composition.get("T"), delta);
		assertEquals(0.0/sequenceLength,composition.get("K"), delta);
		assertEquals(1.0/sequenceLength,composition.get("E"), delta);
		assertEquals(1.0/sequenceLength,composition.get("V"), delta);
		assertEquals(0.0/sequenceLength,composition.get("S"), delta);
		assertEquals(1.0/sequenceLength,composition.get("G"), delta);
		assertEquals(3.0/sequenceLength,composition.get("A"), delta);
		assertEquals(3.0/sequenceLength,composition.get("L"), delta);

		Map<String, Double> iComposition = PeptideProperties.getAACompositionString(fullInvalidSequence);
		assertEquals(0.0,iComposition.get("W"), delta);
		assertEquals(0.0,iComposition.get("C"), delta);
		assertEquals(0.0,iComposition.get("M"), delta);
		assertEquals(0.0,iComposition.get("H"), delta);
		assertEquals(0.0,iComposition.get("Y"), delta);
		assertEquals(0.0,iComposition.get("F"), delta);
		assertEquals(0.0,iComposition.get("Q"), delta);
		assertEquals(0.0,iComposition.get("N"), delta);
		assertEquals(0.0,iComposition.get("I"), delta);
		assertEquals(0.0,iComposition.get("R"), delta);
		assertEquals(0.0,iComposition.get("D"), delta);
		assertEquals(0.0,iComposition.get("P"), delta);
		assertEquals(0.0,iComposition.get("T"), delta);
		assertEquals(0.0,iComposition.get("K"), delta);
		assertEquals(0.0,iComposition.get("E"), delta);
		assertEquals(0.0,iComposition.get("V"), delta);
		assertEquals(0.0,iComposition.get("S"), delta);
		assertEquals(0.0,iComposition.get("G"), delta);
		assertEquals(0.0,iComposition.get("A"), delta);
		assertEquals(0.0,iComposition.get("L"), delta);

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
		assertEquals(2.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "W"), delta);
		assertEquals(3.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "C"), delta);
		assertEquals(1.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "M"), delta);
		assertEquals(1.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "H"), delta);
		assertEquals(0.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "Y"), delta);
		assertEquals(0.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "F"), delta);
		assertEquals(1.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "Q"), delta);
		assertEquals(0.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "N"), delta);
		assertEquals(0.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "I"), delta);
		assertEquals(4.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "R"), delta);
		assertEquals(1.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "D"), delta);
		assertEquals(1.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "P"), delta);
		assertEquals(1.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "T"), delta);
		assertEquals(0.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "K"), delta);
		assertEquals(1.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "E"), delta);
		assertEquals(1.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "V"), delta);
		assertEquals(0.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "S"), delta);
		assertEquals(1.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "G"), delta);
		assertEquals(3.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "A"), delta);
		assertEquals(3.0/sequenceLength, PeptideProperties.getEnrichment(sequence, "L"), delta);

		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "W"), delta);
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "C"), delta);
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "M"), delta);
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "H"), delta);
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "Y"), delta);
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "F"), delta);
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "Q"), delta);
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "N"), delta);
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "I"), delta);
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "R"), delta);
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "D"), delta);
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "P"), delta);
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "T"), delta);
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "K"), delta);
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "E"), delta);
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "V"), delta);
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "S"), delta);
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "G"), delta);
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "A"), delta);
		assertEquals(0.0, PeptideProperties.getEnrichment(fullInvalidSequence, "L"), delta);
		assertEquals(0.0, PeptideProperties.getEnrichment(sequence, "X"), delta);
	}

	@Test (expected = NullPointerException.class)
	public void testEnrichmentNull(){
		assertEquals(0.0, PeptideProperties.getEnrichment(sequence, ""), delta);
		assertNull(PeptideProperties.getEnrichment(sequence, "1"));
	}

	@Test
	public void testMolecularWeight(){
		//http://www.innovagen.se/custom-peptide-synthesis/peptide-property-calculator/peptide-property-calculator.asp
		//http://au.expasy.org/cgi-bin/protparam
		//2872.4 is the value computed by the above two web tools
		assertEquals(2872.4, Utils.roundToDecimals(PeptideProperties.getMolecularWeight(sequence), 1), delta);
		assertEquals(0.0, PeptideProperties.getMolecularWeight("Z"), delta);
		assertEquals(0.0, PeptideProperties.getMolecularWeight("1"), delta);

		assertEquals(0.0, PeptideProperties.getMolecularWeight(fullInvalidSequence), delta);
	}

	@Test
	public void testMolecularWeightXML() throws FileNotFoundException, JAXBException{
		File elementMassFile = new File("./src/main/resources/ElementMass.xml");
		File aminoAcidCompositionFile = new File("./src/main/resources/AminoAcidComposition.xml");

		assertEquals(
				PeptideProperties.getMolecularWeight("A", elementMassFile, aminoAcidCompositionFile)
						* 5.0 -  4 * (17.0073 + 1.0079),
				PeptideProperties.getMolecularWeight("AAAAA", elementMassFile, aminoAcidCompositionFile),
				delta);
	}

	@Test
	public void testMolecularWeightXMLSingleFile() throws FileNotFoundException, JAXBException{
		File aminoAcidCompositionFile = new File("./src/main/resources/AminoAcidComposition.xml");

		assertEquals(
				PeptideProperties.getMolecularWeight("A", aminoAcidCompositionFile) * 5.0 -  4 * (17.0073 + 1.0079),
				PeptideProperties.getMolecularWeight("AAAAA", aminoAcidCompositionFile),
				delta);
	}

	@Test
	public void testMolecularWeightBasedOnAminoAcidCompositionTable() throws Exception{
		File elementMassFile = new File("./src/main/resources/ElementMass.xml");
		File aminoAcidCompositionFile = new File("./src/main/resources/AminoAcidComposition.xml");
		AminoAcidCompositionTable table = PeptideProperties.obtainAminoAcidCompositionTable(elementMassFile, aminoAcidCompositionFile);

		assertEquals(
				PeptideProperties.getMolecularWeightBasedOnXML("A", table) * 5.0 -  4 * (17.0073 + 1.0079),
				PeptideProperties.getMolecularWeightBasedOnXML("AAAAA", table),
				delta);
	}

	@Test (expected = NullPointerException.class)
	public void testMolecularWeightXMLNull() throws FileNotFoundException, JAXBException{
		PeptideProperties.getMolecularWeight(sequence, null, null);
	}

	@Test
	public void testExtinctionCoefficient(){
		//http://au.expasy.org/cgi-bin/protparam
		assertEquals(11125.0, PeptideProperties.getExtinctionCoefficient(sequence, false), delta);
		assertEquals(11000.0, PeptideProperties.getExtinctionCoefficient(sequence, true), delta);

		assertEquals(0.0, PeptideProperties.getExtinctionCoefficient(fullInvalidSequence, true), delta);
		assertEquals(0.0, PeptideProperties.getExtinctionCoefficient(fullInvalidSequence, false), delta);
	}

	@Test (expected = NullPointerException.class)
	public void testExtinctionCoefficientNull(){
		assertEquals(11000.0, PeptideProperties.getExtinctionCoefficient(null, true), delta);
	}

	@Test
	public void testAbsorbance(){
		//http://au.expasy.org/cgi-bin/protparam
		assertEquals(3.830, PeptideProperties.getAbsorbance(sequence, true), 0.001);
		assertEquals(3.873, PeptideProperties.getAbsorbance(sequence, false), 0.001);

		assertEquals(0.0, PeptideProperties.getAbsorbance(fullInvalidSequence, true), 0.001);
		assertEquals(0.0, PeptideProperties.getAbsorbance(fullInvalidSequence, false), 0.001);
	}

	@Test (expected = NullPointerException.class)
	public void testAbsorbanceNull(){
		assertEquals(3.830, PeptideProperties.getAbsorbance(null, false), delta);
	}


	@Test
	public void testInstabilityIndex(){
		//http://au.expasy.org/cgi-bin/protparam
		assertEquals(38.48, PeptideProperties.getInstabilityIndex(sequence), 0.01);
		assertEquals(0.0, PeptideProperties.getInstabilityIndex(fullInvalidSequence), 0.01);
	}

	@Test (expected = NullPointerException.class)
	public void testInstabilityIndexNull(){
		assertEquals(38.48, PeptideProperties.getInstabilityIndex(null), delta);
	}

	@Test
	public void testApliphaticIndex(){
		//http://au.expasy.org/cgi-bin/protparam
		assertEquals(73.33, PeptideProperties.getApliphaticIndex(sequence), 0.01);
		assertEquals(0.0, PeptideProperties.getApliphaticIndex(fullInvalidSequence), 0.01);
	}

	@Test (expected = NullPointerException.class)
	public void testApliphaticIndexNull(){
		assertEquals(73.33, PeptideProperties.getApliphaticIndex(null), 0.01);
	}

	@Test
	public void testAverageHydropathy(){
		//http://au.expasy.org/cgi-bin/protparam
		assertEquals(-0.242, PeptideProperties.getAvgHydropathy(sequence), 0.001);
		assertEquals(0.0, PeptideProperties.getAvgHydropathy(fullInvalidSequence), 0.001);
	}

	@Test (expected = NullPointerException.class)
	public void testAverageHydropathyNull(){
		assertEquals(-0.242, PeptideProperties.getAvgHydropathy(null), 0.001);
	}

	@Test
	public void testIsoelectricPointInnovagen(){
		/*
		 * Test for Innovagen
		 */
		//http://www.innovagen.se/custom-peptide-synthesis/peptide-property-calculator/peptide-property-calculator.asp
		assertEquals(9.01, PeptideProperties.getIsoelectricPoint(sequence, false), 0.01);
		assertEquals(7.00, PeptideProperties.getIsoelectricPoint(fullInvalidSequence, false), 0.01);

		assertEquals(2.70, PeptideProperties.getIsoelectricPoint("ACCACAAADADADACA", false), 0.01);
	}

	@Test
	public void testIsoelectricPointExpasy(){
		/*
		 * Test for Expasy
		 */
		assertEquals(3.42, PeptideProperties.getIsoelectricPoint("ACCACAAADADADACA"), 0.01);
		assertEquals(3.42, PeptideProperties.getIsoelectricPoint("ACCACAAADADADACM"), 0.01);
		//
		assertEquals(3.37, PeptideProperties.getIsoelectricPoint("ECCACAAADADADACS", true), 0.01);


		assertEquals(3.24, PeptideProperties.getIsoelectricPoint("ADCCACAAADADADACDAAAAAAAAAAAA", true), 0.01);

		//3.32 at Expasy
		assertEquals(3.32, PeptideProperties.getIsoelectricPoint("DCCACAAADADADACS", true), 0.01);

		assertEquals(3.17, PeptideProperties.getIsoelectricPoint("DCCACAAADADADACDAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAD", true), 0.01);
		assertEquals(3.37, PeptideProperties.getIsoelectricPoint("ACCACAAADADADACE", true), 0.01);
		assertEquals(3.32, PeptideProperties.getIsoelectricPoint("ACCACAAADADADACAAAAAAAAAAAAAAD", true), 0.01);
		assertEquals(3.28, PeptideProperties.getIsoelectricPoint("DCCACAAADADADACE", true), 0.01);

		assertEquals(8.71, PeptideProperties.getIsoelectricPoint("MTADGPCRELLCQLRAAVRHRWWC", true), 0.01);
		assertEquals(8.71, PeptideProperties.getIsoelectricPoint(sequence, true), 0.01);
		assertEquals(0.0, PeptideProperties.getIsoelectricPoint(fullInvalidSequence, true), 0.01);
	}

	@Test (expected = NullPointerException.class)
	public void testIsoelectricPointNull(){
		assertEquals(8.6, PeptideProperties.getIsoelectricPoint(null), 0.1);
	}

	@Test
	public void testNetCharge(){
		/*
		 * Test for Innovagen
		 */
		//http://www.innovagen.se/custom-peptide-synthesis/peptide-property-calculator/peptide-property-calculator.asp
		assertEquals(2.9, PeptideProperties.getNetCharge(sequence, false), 0.1);
		assertEquals(0.0, PeptideProperties.getNetCharge(fullInvalidSequence, false), 0.1);

		assertEquals(-3.2, PeptideProperties.getNetCharge("ACCACAAADADADACA", false), 0.1);
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
			logger.debug(p);
			logger.debug("pH\tInnovagen\tExpasy\tdiff");
			for ( int i = 1; i < 15; i++) {
				double phPoint = (new Double(i)).doubleValue();
				double chrgInnovagen = PeptideProperties.getNetCharge(p,false,phPoint);
				double chrgExpasy = PeptideProperties.getNetCharge(p,true,phPoint);
				logger.debug(String.format("%2.1f\t%2.2f\t%2.2f\t%2.2f", phPoint, chrgInnovagen,
						chrgExpasy, chrgInnovagen - chrgExpasy));
			}
		}
	}

	@Test (expected = NullPointerException.class)
	public void testNetChargeNull(){
		assertEquals(8.6, PeptideProperties.getNetCharge(null), delta);
	}
}

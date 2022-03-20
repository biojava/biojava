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
package org.biojava.nbio.genome;

import junitx.framework.FileAssert;
import org.biojava.nbio.genome.parsers.gff.FeatureList;
import org.biojava.nbio.genome.parsers.gff.GFF3Reader;
import org.biojava.nbio.genome.parsers.gff.GFF3Writer;
import org.biojava.nbio.core.sequence.ChromosomeSequence;
import org.biojava.nbio.core.sequence.GeneSequence;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.io.FastaWriterHelper;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.FileOutputStream;
import java.util.Collection;
import java.util.LinkedHashMap;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class GeneFeatureHelperTest {

	private static final Logger logger = LoggerFactory.getLogger(GeneFeatureHelperTest.class);

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void testZeroLocation() throws Exception {

		@SuppressWarnings("unused")
		FeatureList listGenes = GFF3Reader.read("src/test/resources/amphimedon.gff3");
	}

	/**
	 * Test of loadFastaAddGeneFeaturesFromUpperCaseExonFastaFile method, of class GeneFeatureHelper.
	 *
	 * @throws Exception
	 */

	@Test
	public void testLoadFastaAddGeneFeaturesFromUpperCaseExonFastaFile() throws Exception {
		// logger.info("loadFastaAddGeneFeaturesFromUpperCaseExonFastaFile");
		File fastaSequenceFile = new File("src/test/resources/volvox_all.fna");
		File uppercaseFastaFile = new File("src/test/resources/volvox_all_genes_exon_uppercase.fna");
		boolean throwExceptionGeneNotFound = false;
		LinkedHashMap<String, ChromosomeSequence> chromosomeSequenceHashMap = GeneFeatureHelper
				.loadFastaAddGeneFeaturesFromUpperCaseExonFastaFile(fastaSequenceFile, uppercaseFastaFile,
						throwExceptionGeneNotFound);

		File tmp = File.createTempFile("volvox_all_genes_exon_uppercase", "gff3");
		tmp.deleteOnExit();
		FileOutputStream fo = new FileOutputStream(tmp);
		GFF3Writer gff3Writer = new GFF3Writer();
		gff3Writer.write(fo, chromosomeSequenceHashMap);
		fo.close();

	}

	/**
	 * Test of outputFastaSequenceLengthGFF3 method, of class GeneFeatureHelper.
	 */
	@Test
	public void testOutputFastaSequenceLengthGFF3() throws Exception {
		// logger.info("outputFastaSequenceLengthGFF3");

		File fastaSequenceFile = new File("src/test/resources/volvox_all.fna");
		File gffFile = File.createTempFile("volvox_length", "gff3");
		gffFile.deleteOnExit();
		GeneGFF3FeatureHelper.outputFastaSequenceLengthGFF3(fastaSequenceFile, gffFile);
		FileAssert.assertEquals("volvox_length.gff3 and volvox_length_output.gff3 are not equal", gffFile,
				new File("src/test/resources/volvox_length_reference.gff3"));

	}

	/**
	 * Test if the note from a gff3 file is added to the gene sequence
	 *
	 * @throws Exception
	 */

	@Test
	public void testAddGFF3Note() throws Exception {
		LinkedHashMap<String, ChromosomeSequence> chromosomeSequenceList = GeneGFF3FeatureHelper
				.loadFastaAddGeneFeaturesFromGmodGFF3(new File("src/test/resources/volvox_all.fna"), new File(
						"src/test/resources/volvox.gff3"), false);
		ChromosomeSequence ctgASequence = chromosomeSequenceList.get("ctgA");
		GeneSequence edenGeneSequence = ctgASequence.getGene("EDEN");
		logger.info("Note {}", edenGeneSequence.getNotesList());
	}

	/**
	 * Test of getProteinSequences method, of class GeneFeatureHelper. Used gff3 file that was modified from the volvox
	 * gff version. Do not have the reference protein that is generated from each CDS record so subject to being
	 * incorrect without a validated test case. Could not find anyone providing a gff3 test case with expected protein
	 * output.
	 */
	@Test
	public void testGetProteinSequences() throws Exception {
		LinkedHashMap<String, ChromosomeSequence> chromosomeSequenceList = GeneGFF3FeatureHelper
				.loadFastaAddGeneFeaturesFromGmodGFF3(new File("src/test/resources/volvox_all.fna"), new File(
						"src/test/resources/volvox.gff3"), false);
		LinkedHashMap<String, ProteinSequence> proteinSequenceList = GeneFeatureHelper
				.getProteinSequences(chromosomeSequenceList.values());
		// for(ProteinSequence proteinSequence : proteinSequenceList.values()){
		// logger.info("Output={}", proteinSequence.getSequenceAsString());
		// }
		File tmp = File.createTempFile("volvox_all", "faa");
		tmp.deleteOnExit();
		FastaWriterHelper.writeProteinSequence(tmp, proteinSequenceList.values());
		FileAssert.assertEquals("volvox_all_reference.faa and volvox_all.faa are not equal", new File(
				"src/test/resources/volvox_all_reference.faa"), tmp);
	}

	/**
	 * Test of getGeneSequences method, of class GeneFeatureHelper.
	 */
	@Test
	public void testGetGeneSequences() throws Exception {
		// logger.info("getGeneSequences");
		LinkedHashMap<String, ChromosomeSequence> chromosomeSequenceList = GeneGFF3FeatureHelper
				.loadFastaAddGeneFeaturesFromGmodGFF3(new File("src/test/resources/volvox_all.fna"), new File(
						"src/test/resources/volvox.gff3"), true);
		LinkedHashMap<String, GeneSequence> geneSequenceHashMap = GeneFeatureHelper
				.getGeneSequences(chromosomeSequenceList.values());
		Collection<GeneSequence> geneSequences = geneSequenceHashMap.values();

		File tmp = File.createTempFile("volvox_all_genes_exon_uppercase", "fna");
		tmp.deleteOnExit();
		FastaWriterHelper.writeGeneSequence(tmp, geneSequences, true);
	}

}

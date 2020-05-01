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
package org.biojava.nbio.core.sequence.loader;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.features.AbstractFeature;
import org.biojava.nbio.core.sequence.features.FeatureInterface;
import org.biojava.nbio.core.sequence.features.Qualifier;
import org.biojava.nbio.core.sequence.template.AbstractSequence;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.channels.Channels;
import java.nio.channels.ReadableByteChannel;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Map;

/**
 * Testing example for issue #834
 *
 * @author Jacek Grzebyta
 * @author Paolo Pavan
 */
@RunWith(Parameterized.class)
public class GenbankProxySequenceReaderTest {

	private String gi;
	private final static Logger logger = LoggerFactory.getLogger(GenbankProxySequenceReaderTest.class);

	public GenbankProxySequenceReaderTest(String gi) {
		this.gi = gi;
	}

	@Parameterized.Parameters
	public static Collection<String[]> getExamples() {
		String[][] accessorIds = new String[][]{
			{"399235158"},
			{"7525057"},
			{"379015144"},
			{"381353147"},
			{"381353148"},
			{"152970917"},
			{"381353149"},
			{"254839678"}
		};

		return Arrays.asList(accessorIds);
	}

	/**
	 * In {@link GenbankProxySequenceReader} there is a check to see if the requested files are already in the temp
	 * directory before attempting to retrieve them from the remote server. so simply copying the test files to the temp
	 * directory avoids calling out to the server and hitting a 429 status code from the server which fails the build.
	 * @throws IOException
	 */
	@Before
	public void copyTestFiles() throws IOException {
		Collection<String[]> accessorIds = getExamples();
		for (String[] arr: accessorIds) {
			copyTestFileToWorkingDirectory(arr[0]+".gb");
		}
	}

	/**
	 * Convenience method for {@link GenbankProxySequenceReaderTest#copyTestFiles()}
	 * @param filename name of the file to copy from the resource folder
	 * @throws IOException when something goes wrong with copying the files.
	 */
	private void copyTestFileToWorkingDirectory(String filename) throws IOException {
		String destRoot = System.getProperty("java.io.tmpdir");

		//if the directory does not end with a slash or backslash then add one
		if(!(destRoot.endsWith("/") || destRoot.endsWith("\\"))){
			destRoot += destRoot.contains("/")? "/" : "\\";
		}

		String dest =  destRoot + filename;
		String src = "org/biojava/nbio/core/sequence/loader/" + filename;

		//Remove any pre-existing files
		File d = new File(dest);
		d.delete();

		try(FileOutputStream destination = new FileOutputStream(d);
		InputStream is = this.getClass().getClassLoader().getResourceAsStream(src);
		ReadableByteChannel source = Channels.newChannel(is)) {

			destination.getChannel().transferFrom(source, 0, Long.MAX_VALUE);
		}
	}


	@Test
	public void testFeatures() throws IOException, InterruptedException, CompoundNotFoundException {
		logger.info("run test for protein: {}", gi);
		GenbankProxySequenceReader<AminoAcidCompound> genbankReader
				= new GenbankProxySequenceReader<>(System.getProperty("java.io.tmpdir"),
						this.gi,
						AminoAcidCompoundSet.getAminoAcidCompoundSet());

		// why only tests on protein sequences?
		ProteinSequence seq = new ProteinSequence(genbankReader);

		Assert.assertNotNull("protein sequence is null", seq);

		/*
		 parse description from header. There is no separate interface/abstract class for method getHeader()
		 so it should be done here (manualy).
		 */
		genbankReader.getHeaderParser().parseHeader(genbankReader.getHeader(), seq);

		// test description
		Assert.assertNotNull(seq.getDescription());

		// test accession Id
		logger.info("accession id: {}", seq.getAccession().getID());
		Assert.assertNotNull(seq.getAccession().getID());
		// test GID number
		if( seq.getAccession().getIdentifier() != null) { // GI: in header now optional. See #596
			Assert.assertEquals(gi, seq.getAccession().getIdentifier());
			logger.info("found identifier '{}'", seq.getAccession().getIdentifier());
		}
		// test taxonomy id
		logger.info("taxonomy id: {}", seq.getTaxonomy().getID());
		Assert.assertNotNull(seq.getTaxonomy().getID());
		Assert.assertNotNull(Integer.decode(seq.getTaxonomy().getID().split(":")[1]));

		// test taxonomy name
		String taxonName = seq.getFeaturesByType("source").get(0).getQualifiers().get("organism").get(0).getValue();
		logger.info("taxonomy name '{}'", taxonName);
		Assert.assertNotNull(taxonName);

		if (seq.getFeaturesByType("CDS").size() > 0) {
			FeatureInterface<AbstractSequence<AminoAcidCompound>, AminoAcidCompound> CDS = seq.getFeaturesByType("CDS").get(0);
			logger.info("CDS: {}", CDS);
			String codedBy = CDS.getQualifiers().get("coded_by").get(0).getValue();
			Assert.assertNotNull(codedBy);
			Assert.assertTrue(!codedBy.isEmpty());
			logger.info("\t\tcoded_by: {}", codedBy);
		}

		// genbank has limits on requests per second, we need to give it some time for next test or otherwise we get 429 http error codes - JD 2018-12-14
		// See https://github.com/biojava/biojava/issues/837
		Thread.sleep(500);
	}

	@Test
	public void testProteinSequenceFactoring() throws Exception {
		logger.info("create protein sequence test for target {}", gi);

		GenbankProxySequenceReader<AminoAcidCompound> genbankReader
				= new GenbankProxySequenceReader<>(System.getProperty("java.io.tmpdir"),
						this.gi,
						AminoAcidCompoundSet.getAminoAcidCompoundSet());

		ProteinSequence seq = new ProteinSequence(genbankReader);

		// if target protein contain CDS/coded_by than it should contain parent nucleotide seq
		List<AbstractFeature<AbstractSequence<AminoAcidCompound>, AminoAcidCompound>> CDSs = genbankReader.getFeatures().get("CDS");

		if (CDSs != null) {
			if (CDSs.size() == 1) {
				final Map<String, List<Qualifier>> qualifiers = CDSs.get(0).getQualifiers();
				List<Qualifier> codedByQualifiers = qualifiers.get("coded_by");
				Qualifier codedBy = codedByQualifiers.get(0);
				if (codedBy != null) {

					AbstractSequence<?> parentSeq = seq.getParentSequence();
					Assert.assertNotNull(parentSeq);

					/*
					 Sometimes protein might have many 'parents' with different accessions
					 so accession is not set.

					 That test is always failed
					 */
					//Assert.assertTrue(parentSeq.getAccession());
					Assert.assertTrue(!parentSeq.getSequenceAsString().isEmpty());
				}
			}
		} else {
			logger.info("target {} has no CDS", gi);
		}

		// genbank has limits on requests per second, we need to give it some time for next test or otherwise we get 429 http error codes - JD 2018-12-14
		// See https://github.com/biojava/biojava/issues/837
		Thread.sleep(500);

	}
}

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
import org.biojava.nbio.core.sequence.features.FeatureInterface;
import org.biojava.nbio.core.sequence.template.AbstractSequence;
import org.junit.Assert;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import org.biojava.nbio.core.sequence.features.AbstractFeature;
import org.biojava.nbio.core.sequence.features.Qualifier;

/**
 * Testing example for issue #834
 *
 * @author Jacek Grzebyta
 * @author Paolo Pavan
 * @see InfoTask
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
		String[][] out = new String[][]{
			{"399235158"},
			{"7525057"},
			{"379015144"},
			{"381353147"},
			{"381353148"},
			{"152970917"},
			{"381353149"},
			{"254839678"}
		};

		return Arrays.asList(out);
	}

	@Test
	public void testFeatures() throws IOException, InterruptedException, CompoundNotFoundException {
		logger.info("run test for protein: {}", gi);
		GenbankProxySequenceReader<AminoAcidCompound> genbankReader
				= new GenbankProxySequenceReader<AminoAcidCompound>(System.getProperty("java.io.tmpdir"),
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
		Assert.assertTrue(seq.getDescription() != null);

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
	}

	@Test
	public void testProteinSequenceFactoring() throws Exception {
		logger.info("create protein sequence test for target {}", gi);

		GenbankProxySequenceReader<AminoAcidCompound> genbankReader
				= new GenbankProxySequenceReader<AminoAcidCompound>(System.getProperty("java.io.tmpdir"),
						this.gi,
						AminoAcidCompoundSet.getAminoAcidCompoundSet());

		ProteinSequence seq = new ProteinSequence(genbankReader);

		// if target protein contain CDS/coded_by than it should contain parent nucleotide seq
		ArrayList<AbstractFeature> CDSs = genbankReader.getFeatures().get("CDS");

		if (CDSs != null) {
			if (CDSs.size() == 1) {
				ArrayList<Qualifier> qualifiers = (ArrayList)CDSs.get(0).getQualifiers().get("coded_by");
				Qualifier codedBy = qualifiers.get(0);
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

	}
}

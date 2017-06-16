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
package org.biojava.nbio.core.sequence.io;

import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.compound.DNACompoundSet;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.features.FeatureInterface;
import org.biojava.nbio.core.sequence.features.Qualifier;
import org.biojava.nbio.core.sequence.template.AbstractSequence;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Assert;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import static org.hamcrest.CoreMatchers.is;
import static org.junit.Assert.*;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 * @author Jacek Grzebyta
 * @author Philippe Soares
 */
public class GenbankReaderTest {

	private final static Logger logger = LoggerFactory.getLogger(GenbankReaderTest.class);

	public GenbankReaderTest() {
	}

	@BeforeClass
	public static void setUpClass() throws Exception {
	}

	@AfterClass
	public static void tearDownClass() throws Exception {
	}

	@Before
	public void setUp() {
	}

	@After
	public void tearDown() {
	}

	/**
	 * Test of process method, of class GenbankReader.
	 */
	@Test
	public void testProcess() throws Exception {

		logger.info("process protein");
		InputStream inStream = this.getClass().getResourceAsStream("/BondFeature.gb");
		assertNotNull(inStream);

		GenbankReader<ProteinSequence, AminoAcidCompound> genbankProtein
				= new GenbankReader<>(
						inStream,
						new GenericGenbankHeaderParser<>(),
						new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet())
				);

		LinkedHashMap<String, ProteinSequence> proteinSequences = genbankProtein.process();

		assertThat(proteinSequences.get("NP_000257").getComments().get(0),is(
				"VALIDATED REFSEQ: This record has undergone validation or\n" +
				"preliminary review. The reference sequence was derived from\n" +
				"AL034370.1, X65882.1 and BE139596.1.\n" +
				"Summary: NDP is the genetic locus identified as harboring mutations\n" +
				"that result in Norrie disease. Norrie disease is a rare genetic\n" +
				"disorder characterized by bilateral congenital blindness that is\n" +
				"caused by a vascularized mass behind each lens due to a\n" +
				"maldeveloped retina (pseudoglioma).\n" +
				"Publication Note:  This RefSeq record includes a subset of the\n" +
				"publications that are available for this gene. Please see the\n" +
				"Entrez Gene record to access additional publications."));

		assertThat(proteinSequences.get("NP_000257").getReferences().size(),is(11));
		assertThat(proteinSequences.get("NP_000257").getReferences().get(0).getAuthors(),
				is("Lev,D., Weigl,Y., Hasan,M., Gak,E., Davidovich,M., Vinkler,C.,\n" +
						"Leshinsky-Silver,E., Lerman-Sagie,T. and Watemberg,N."));
		assertThat(proteinSequences.get("NP_000257").getReferences().get(1).getTitle(),
				is("Novel mutations in Norrie disease gene in Japanese patients with\n" +
						"Norrie disease and familial exudative vitreoretinopathy"));
		assertThat(proteinSequences.get("NP_000257").getReferences().get(10).getJournal(),
				is("Nat. Genet. 1 (3), 199-203 (1992)"));

		assertNotNull(proteinSequences);
		assertEquals(1, proteinSequences.size());

		ProteinSequence proteinSequence = proteinSequences.get("NP_000257");
		assertNotNull(proteinSequences.get("NP_000257"));
		assertEquals("NP_000257", proteinSequence.getAccession().getID());
		assertEquals("4557789", proteinSequence.getAccession().getIdentifier());
		assertEquals("GENBANK", proteinSequence.getAccession().getDataSource().name());
		assertEquals(1, proteinSequence.getAccession().getVersion().intValue());
		assertTrue(genbankProtein.isClosed());

		logger.info("process DNA");
		inStream = this.getClass().getResourceAsStream("/NM_000266.gb");
		assertNotNull(inStream);

		GenbankReader<DNASequence, NucleotideCompound> genbankDNA
				= new GenbankReader<>(
						inStream,
						new GenericGenbankHeaderParser<>(),
						new DNASequenceCreator(DNACompoundSet.getDNACompoundSet())
				);
		LinkedHashMap<String, DNASequence> dnaSequences = genbankDNA.process();

		assertNotNull(dnaSequences);
		assertEquals(1, dnaSequences.size());

		DNASequence dnaSequence = dnaSequences.get("NM_000266");
		assertNotNull(dnaSequences.get("NM_000266"));
		assertEquals("NM_000266", dnaSequence.getAccession().getID());
		assertEquals("223671892", dnaSequence.getAccession().getIdentifier());
		assertEquals("GENBANK", dnaSequence.getAccession().getDataSource().name());
		assertEquals(3, dnaSequence.getAccession().getVersion().intValue());
		assertTrue(genbankDNA.isClosed());
	}

	/**
	 * Test the process method with a number of sequences to be read at each call.
	 * The underlying {@link InputStream} should remain open until the last call.
	 */
	@Test
	public void testPartialProcess() throws IOException, CompoundNotFoundException, NoSuchFieldException {
		InputStream inStream = this.getClass().getResourceAsStream("/two-dnaseqs.gb");

		GenbankReader<DNASequence, NucleotideCompound> genbankDNA
				= new GenbankReader<>(
				inStream,
				new GenericGenbankHeaderParser<>(),
				new DNASequenceCreator(DNACompoundSet.getDNACompoundSet())
		);

		// First call to process(1) returns the first sequence
		LinkedHashMap<String, DNASequence> dnaSequences = genbankDNA.process(1);

		assertNotNull(dnaSequences);
		assertEquals(1, dnaSequences.size());
		assertNotNull(dnaSequences.get("vPetite"));

		// Second call to process(1) returns the second sequence
		dnaSequences = genbankDNA.process(1);
		assertNotNull(dnaSequences);
		assertEquals(1, dnaSequences.size());
		assertNotNull(dnaSequences.get("sbFDR"));

		assertFalse(genbankDNA.isClosed());
		genbankDNA.close();
		assertTrue(genbankDNA.isClosed());

	}

	@Test
	public void CDStest() throws Exception {
		logger.info("CDS Test");

		InputStream inStream = this.getClass().getResourceAsStream("/BondFeature.gb");
		assertNotNull(inStream);

		GenbankReader<ProteinSequence, AminoAcidCompound> GenbankProtein
				= new GenbankReader<>(
						inStream,
						new GenericGenbankHeaderParser<>(),
						new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet())
				);
		LinkedHashMap<String, ProteinSequence> proteinSequences = GenbankProtein.process();
		inStream.close();


		Assert.assertTrue(proteinSequences.size() == 1);
		logger.debug("protein sequences: {}", proteinSequences);

		ProteinSequence protein = new ArrayList<>(proteinSequences.values()).get(0);

		FeatureInterface<AbstractSequence<AminoAcidCompound>, AminoAcidCompound> cdsFeature = protein.getFeaturesByType("CDS").get(0);
		String codedBy = cdsFeature.getQualifiers().get("coded_by").get(0).getValue();
		Map<String, List<Qualifier>> quals = cdsFeature.getQualifiers();
		List<Qualifier> dbrefs = quals.get("db_xref");

		Assert.assertNotNull(codedBy);
		Assert.assertTrue(!codedBy.isEmpty());
		assertEquals(codedBy, "NM_000266.2:503..904");
		assertEquals(5, dbrefs.size());

	}

}

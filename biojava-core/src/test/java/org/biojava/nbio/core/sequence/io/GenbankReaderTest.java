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

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.RNASequence;
import org.biojava.nbio.core.sequence.Strand;
import org.biojava.nbio.core.sequence.compound.*;
import org.biojava.nbio.core.sequence.features.FeatureInterface;
import org.biojava.nbio.core.sequence.features.Qualifier;
import org.biojava.nbio.core.sequence.location.template.AbstractLocation;
import org.biojava.nbio.core.sequence.template.AbstractSequence;
import org.biojava.nbio.core.sequence.template.Compound;
import org.junit.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.*;

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
	public void testPartialProcess() throws IOException, CompoundNotFoundException {
		CheckableInputStream inStream = new CheckableInputStream(this.getClass().getResourceAsStream("/two-dnaseqs.gb"));

		GenbankReader<DNASequence, NucleotideCompound> genbankDNA
				= new GenbankReader<>(
				inStream,
				new GenericGenbankHeaderParser<>(),
				new DNASequenceCreator(DNACompoundSet.getDNACompoundSet())
		);

		// First call to process(1) returns the first sequence
		LinkedHashMap<String, DNASequence> dnaSequences = genbankDNA.process(1);

		assertFalse(inStream.isclosed());
		assertNotNull(dnaSequences);
		assertEquals(1, dnaSequences.size());
		assertNotNull(dnaSequences.get("vPetite"));

		// Second call to process(1) returns the second sequence
		dnaSequences = genbankDNA.process(1);
		assertFalse(inStream.isclosed());
		assertNotNull(dnaSequences);
		assertEquals(1, dnaSequences.size());
		assertNotNull(dnaSequences.get("sbFDR"));

		assertFalse(genbankDNA.isClosed());
		genbankDNA.close();
		assertTrue(genbankDNA.isClosed());
		assertTrue(inStream.isclosed());
	}

	@Test
	public void CDStest() throws Exception {
		logger.info("CDS Test");

		CheckableInputStream inStream = new CheckableInputStream(this.getClass().getResourceAsStream("/BondFeature.gb"));
		assertNotNull(inStream);

		GenbankReader<ProteinSequence, AminoAcidCompound> GenbankProtein
				= new GenbankReader<>(
						inStream,
						new GenericGenbankHeaderParser<>(),
						new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet())
				);
		LinkedHashMap<String, ProteinSequence> proteinSequences = GenbankProtein.process();
		assertTrue(inStream.isclosed());


		Assert.assertEquals(1, proteinSequences.size());
		logger.debug("protein sequences: {}", proteinSequences);

		ProteinSequence protein = new ArrayList<>(proteinSequences.values()).get(0);

		FeatureInterface<AbstractSequence<AminoAcidCompound>, AminoAcidCompound> cdsFeature = protein.getFeaturesByType("CDS").get(0);
		String codedBy = cdsFeature.getQualifiers().get("coded_by").get(0).getValue();
		Map<String, List<Qualifier>> quals = cdsFeature.getQualifiers();
		List<Qualifier> dbrefs = quals.get("db_xref");

		Assert.assertNotNull(codedBy);
		Assert.assertFalse(codedBy.isEmpty());
		assertEquals("NM_000266.2:503..904", codedBy);
		assertEquals(5, dbrefs.size());

	}

	private DNASequence readGenbankResource(final String resource) throws IOException, CompoundNotFoundException {
		InputStream inputStream = getClass().getResourceAsStream(resource);
		GenbankReader<DNASequence, NucleotideCompound> genbankDNA
		= new GenbankReader<>(
				inputStream,
				new GenericGenbankHeaderParser<>(),
				new DNASequenceCreator(DNACompoundSet.getDNACompoundSet())
				);
		LinkedHashMap<String, DNASequence> dnaSequences = genbankDNA.process();
		return dnaSequences.values().iterator().next();
	}
	
	private RNASequence readGenbankRNAResource(final String resource) throws IOException, CompoundNotFoundException {
		InputStream inputStream = getClass().getResourceAsStream(resource);
		GenbankReader<RNASequence, NucleotideCompound> genbankRNA
		= new GenbankReader<>(
				inputStream,
				new GenericGenbankHeaderParser<>(),
				new RNASequenceCreator(RNACompoundSet.getRNACompoundSet())
				);
		LinkedHashMap<String, RNASequence> rnaSequences = genbankRNA.process();
		return rnaSequences.values().iterator().next();	
	}
	
	private ProteinSequence readGenbankProteinResource(final String resource) throws IOException, CompoundNotFoundException {
		InputStream inputStream = getClass().getResourceAsStream(resource);
		GenbankReader<ProteinSequence, AminoAcidCompound> genbankProtein
		= new GenbankReader<>(
				inputStream,
				new GenericGenbankHeaderParser<>(),
				new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet())
				);
		LinkedHashMap<String, ProteinSequence> proteinSequences = genbankProtein.process();
		return proteinSequences.values().iterator().next();	
	}
	
	private AbstractSequence<?> readUnknownGenbankResource(final String resource) throws IOException, CompoundNotFoundException {
		InputStream inputStream = getClass().getResourceAsStream(resource);
		GenbankSequenceParser<AbstractSequence<Compound>, Compound> genbankParser = new GenbankSequenceParser<>();
		BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(inputStream));
		String seqString = genbankParser.getSequence(bufferedReader, 0);
		String compoundSet = genbankParser.getCompoundType().getClass().getSimpleName();

		if (compoundSet.equals("AminoAcidCompoundSet")) {
			return readGenbankProteinResource(resource);
		} else if (compoundSet.equals("RNACompoundSet")) {
			return readGenbankRNAResource(resource);
		} else {
			return readGenbankResource(resource);
		}
	}

	@Test
	public void testNcbiExpandedAccessionFormats() throws IOException, CompoundNotFoundException {
		DNASequence header0 = readGenbankResource("/empty_header0.gb");
		assertEquals("CP032762             5868661 bp    DNA     circular BCT 15-OCT-2018", header0.getOriginalHeader());

		DNASequence header1 = readGenbankResource("/empty_header1.gb");
		assertEquals("AZZZAA02123456789 9999999999 bp    DNA     linear   PRI 15-OCT-2018", header1.getOriginalHeader());

		DNASequence header2 = readGenbankResource("/empty_header2.gb");
		assertEquals("AZZZAA02123456789 10000000000 bp    DNA     linear   PRI 15-OCT-2018", header2.getOriginalHeader());
	}
	
	@Test
	public void testLegacyLocusCompatable() throws IOException, CompoundNotFoundException {
		
		// Testing opening a genbank file with uppercase units, strand and topology
		AbstractSequence<? extends Compound> header0 = readUnknownGenbankResource("/org/biojava/nbio/core/sequence/io/uppercase_locus0.gb");
		assertEquals("ABC12.3_DE   7071 BP DS-DNA   CIRCULAR  SYN       22-JUL-1994", header0.getOriginalHeader());
		assertEquals("ABC12.3_DE", header0.getAccession().getID());
		assertEquals("DNACompoundSet", header0.getCompoundSet().getClass().getSimpleName());
		
		// Testing uppercase SS strand
		AbstractSequence<? extends Compound> header1 = readUnknownGenbankResource("/org/biojava/nbio/core/sequence/io//uppercase_locus1.gb");
		assertEquals("ABC12.3_DE   7071 BP SS-DNA   CIRCULAR  SYN       13-JUL-1994", header1.getOriginalHeader());
		assertEquals("ABC12.3_DE", header1.getAccession().getID());
		assertEquals("DNACompoundSet", header0.getCompoundSet().getClass().getSimpleName());
		
		// Testing uppercase MS strand
		AbstractSequence<? extends Compound> header2 = readUnknownGenbankResource("/org/biojava/nbio/core/sequence/io//uppercase_locus2.gb");
		assertEquals("ABC12.3_DE   7071 BP MS-DNA   CIRCULAR  SYN       13-JUL-1994", header2.getOriginalHeader());
		assertEquals("ABC12.3_DE", header2.getAccession().getID());
		assertEquals("DNACompoundSet", header0.getCompoundSet().getClass().getSimpleName());
		
		// Testing uppercase LINEAR topology
		AbstractSequence<? extends Compound> header3 = readUnknownGenbankResource("/org/biojava/nbio/core/sequence/io//uppercase_locus3.gb");
		assertEquals("ABC12.3_DE   7071 BP DNA   LINEAR    SYN       22-JUL-1994", header3.getOriginalHeader());
		assertEquals("ABC12.3_DE", header3.getAccession().getID());
		assertEquals("DNACompoundSet", header0.getCompoundSet().getClass().getSimpleName());
		
		// Testing uppercase units with no strand or topology
		AbstractSequence<? extends Compound> header4 = readUnknownGenbankResource("/org/biojava/nbio/core/sequence/io//uppercase_locus4.gb");
		assertEquals("ABC12.3_DE   7071 BP RNA             SYN       13-JUL-1994", header4.getOriginalHeader());
		assertEquals("ABC12.3_DE", header4.getAccession().getID());
		assertEquals("RNACompoundSet", header4.getCompoundSet().getClass().getSimpleName());
		
		// Testing uppercase units with no strand, topology, division or date
		AbstractSequence<? extends Compound> header5 = readUnknownGenbankResource("/org/biojava/nbio/core/sequence/io//uppercase_locus5.gb");
		assertEquals("ABC12.3_DE   7071 BP DNA", header5.getOriginalHeader());
		assertEquals("ABC12.3_DE", header5.getAccession().getID());
		
		// Testing uppercase units with no strand, molecule type, topology, division or date
		AbstractSequence<? extends Compound> header6 = readUnknownGenbankResource("/org/biojava/nbio/core/sequence/io//uppercase_locus6.gb");
		assertEquals("ABC12.3_DE   7071 BP", header6.getOriginalHeader());
		assertEquals("ABC12.3_DE", header6.getAccession().getID());
		assertEquals("DNACompoundSet", header0.getCompoundSet().getClass().getSimpleName());
		
		// Testing uppercase protein units
		AbstractSequence<? extends Compound> header7 = readUnknownGenbankResource("/org/biojava/nbio/core/sequence/io//uppercase_locus7.gb");
		assertEquals("ABC12.3_DE   7071 AA Protein", header7.getOriginalHeader());
		assertEquals("ABC12.3_DE", header7.getAccession().getID());
		assertEquals("AminoAcidCompoundSet", header7.getCompoundSet().getClass().getSimpleName());
		
	}

	@Test
	public void readSequenceWithZeroSpanFeature() throws IOException, CompoundNotFoundException {
		logger.info("make or read genbank file error when feature spans zero point of circular sequence (issue #855)");
		final DNASequence seq = readGenbankResource("/feature-spans-zero-point-circular-sequence.gb");

		assertNotNull(seq);

		final FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound> f = seq.getFeatures().get(33);
		final AbstractLocation fLocation = f.getLocations();

		assertTrue(fLocation.isCircular());
		assertEquals(7028, (int)fLocation.getStart().getPosition());
		assertEquals(286, (int)fLocation.getEnd().getPosition());
		assertEquals(Strand.NEGATIVE, fLocation.getStrand());
	}

	/**
	 * Biojava fails to parse anticodon and transl_except feature qualifiers when they line wrap.
	 * https://github.com/biojava/biojava/issues/843
	 */
	@Test
	public void testGithub843() throws Exception {
		CheckableInputStream inStream = new CheckableInputStream(this.getClass().getResourceAsStream("/NC_018080.gb"));
		assertNotNull(inStream);

		GenbankReader<DNASequence, NucleotideCompound> genbankDNA
				= new GenbankReader<>(
				inStream,
				new GenericGenbankHeaderParser<>(),
				new DNASequenceCreator(DNACompoundSet.getDNACompoundSet())
		);

		LinkedHashMap<String, DNASequence> dnaSequences = genbankDNA.process();
		assertNotNull(dnaSequences);

		DNASequence dna = new ArrayList<>(dnaSequences.values()).get(0);
		assertNotNull(dna);

		FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound> tRNAFeature = dna.getFeaturesByType("tRNA").get(0);
		String anticodon = tRNAFeature.getQualifiers().get("anticodon").get(0).getValue();
		assertEquals("(pos:complement(1123552..1123554),aa:Leu,seq:caa)", anticodon);
		String transl_except = tRNAFeature.getQualifiers().get("transl_except").get(0).getValue();
		assertEquals("(pos:complement(1123552..1123554),aa:Leu)",transl_except);
	}

	/**
	 * Helper class to be able to verify the closed state of the input stream.
	 */
	private static class CheckableInputStream extends BufferedInputStream {

		private boolean closed;

		CheckableInputStream(InputStream in) {
			super(in);
			closed = false;
		}

		@Override
		public void close() throws IOException {
			super.close();
			closed = true;
		}

		boolean isclosed() {
			return closed;
		}
	}
}

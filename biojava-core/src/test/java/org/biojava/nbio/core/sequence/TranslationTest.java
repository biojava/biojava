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
package org.biojava.nbio.core.sequence;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.compound.DNACompoundSet;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.io.*;
import org.biojava.nbio.core.sequence.io.util.ClasspathResource;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.biojava.nbio.core.sequence.transcription.Frame;
import org.biojava.nbio.core.sequence.transcription.TranscriptionEngine;
import org.biojava.nbio.core.sequence.transcription.TranscriptionEngine.Builder;
import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.InputStream;
import java.util.EnumMap;
import java.util.Map;
import java.util.Map.Entry;

import static org.biojava.nbio.core.sequence.io.util.IOUtils.close;
import static org.hamcrest.CoreMatchers.is;
import static org.junit.Assert.*;

public class TranslationTest {

	private final static Logger logger = LoggerFactory.getLogger(TranslationTest.class);

	private static DNACompoundSet dnaCs = DNACompoundSet.getDNACompoundSet();
	private static AminoAcidCompoundSet aaCs = AminoAcidCompoundSet.getAminoAcidCompoundSet();
	private static DNASequence brca2Dna;
	private static Sequence<AminoAcidCompound> brca2Pep;
	private static Sequence<NucleotideCompound> volvoxDna;
	private static Sequence<AminoAcidCompound> volvoxPep;

	@BeforeClass
	public static void parseSequences() {
		InputStream cdsIs = new ClasspathResource(
				"org/biojava/nbio/core/sequence/BRCA2-cds.fasta").getInputStream();
		InputStream pepIs = new ClasspathResource(
				"org/biojava/nbio/core/sequence/BRCA2-peptide.fasta").getInputStream();
		InputStream volDnaIs = new ClasspathResource(
				"org/biojava/nbio/core/sequence/volvox-cds.fasta").getInputStream();
		InputStream volPepIs = new ClasspathResource(
				"org/biojava/nbio/core/sequence/volvox-peptide.fasta").getInputStream();

		try {
			FastaReader<DNASequence, NucleotideCompound> dnaReader = new FastaReader<DNASequence, NucleotideCompound>(cdsIs,
					new GenericFastaHeaderParser<DNASequence, NucleotideCompound>(), new DNASequenceCreator(dnaCs));
			brca2Dna = dnaReader.process().values().iterator().next();
			FastaReader<ProteinSequence, AminoAcidCompound> pReader = new FastaReader<ProteinSequence, AminoAcidCompound>(
					pepIs, new GenericFastaHeaderParser<ProteinSequence, AminoAcidCompound>(), new ProteinSequenceCreator(
							aaCs));
			brca2Pep = pReader.process().values().iterator().next();

			FastaReader<DNASequence, NucleotideCompound> volvoxDnaReader = new FastaReader<DNASequence, NucleotideCompound>(volDnaIs,
					new GenericFastaHeaderParser<DNASequence, NucleotideCompound>(), new DNASequenceCreator(dnaCs));
			volvoxDna = volvoxDnaReader.process().values().iterator().next();
			FastaReader<ProteinSequence, AminoAcidCompound> volvoxPepReader = new FastaReader<ProteinSequence, AminoAcidCompound>(
					volPepIs, new GenericFastaHeaderParser<ProteinSequence, AminoAcidCompound>(), new ProteinSequenceCreator(
							aaCs));
			volvoxPep = volvoxPepReader.process().values().iterator().next();
		}
		catch (IOException e) {
			logger.error("Exception: ", e);
			Assert.fail("Encountered exception");
		}
		finally {
			close(cdsIs);
			close(pepIs);
			close(volDnaIs);
			close(volPepIs);
		}
	}

	@Test
	public void getUniversal() {
		IUPACParser.getInstance().getTable(1);
		IUPACParser.getInstance().getTable("UNIVERSAL");
	}

	@Test
	public void basicTranslation() throws CompoundNotFoundException {
		TranscriptionEngine e = TranscriptionEngine.getDefault();
		DNASequence dna = new DNASequence("ATG");
		RNASequence rna = dna.getRNASequence(e);
		ProteinSequence protein = rna.getProteinSequence(e);
		AminoAcidCompound initMet = protein.getCompoundAt(1);
		assertThat("Initator methionine wrong", initMet.toString(), is("M"));
	}

	@Test
	public void translateN() throws CompoundNotFoundException {
		TranscriptionEngine.Builder b = new TranscriptionEngine.Builder();
		b.translateNCodons(true).initMet(true);
		TranscriptionEngine e = b.build();
		DNASequence dna = new DNASequence("ATN");
		RNASequence rna = dna.getRNASequence(e);
		ProteinSequence protein = rna.getProteinSequence(e);
		assertThat("Ambiguous translation problem", protein.toString(), is("X"));
		DNASequence dna2 = new DNASequence("GTGGTNTAA");
		RNASequence rna2 = dna2.getRNASequence(e);
		ProteinSequence protein2 = rna2.getProteinSequence(e);
		assertThat("Ambiguous translation problem", protein2.toString(), is("VX"));
	}

	@SuppressWarnings("serial")
	@Test
	public void multiFrameTranslation() throws CompoundNotFoundException {
		TranscriptionEngine e = TranscriptionEngine.getDefault();
		DNASequence dna = new DNASequence("ATGGCGTGA");

		Map<Frame, String> expectedTranslations = new EnumMap<Frame, String>(Frame.class) {

			{
				put(Frame.ONE, "MA");
				put(Frame.TWO, "WR");
				put(Frame.THREE, "GV");
				put(Frame.REVERSED_ONE, "SRH");
				put(Frame.REVERSED_TWO, "HA");
				put(Frame.REVERSED_THREE, "TP");
			}
		};

		Map<Frame, Sequence<AminoAcidCompound>> translations =
				e.multipleFrameTranslation(dna, Frame.getAllFrames());

		for (Entry<Frame, Sequence<AminoAcidCompound>> entry : translations.entrySet()) {
			String expected = expectedTranslations.get(entry.getKey());
			Sequence<AminoAcidCompound> protein = entry.getValue();
			assertThat("Checking 6 frame translation", protein.toString(), is(expected));
		}
	}

	@Test
	public void lowerCases() throws CompoundNotFoundException {
		TranscriptionEngine e = TranscriptionEngine.getDefault();
		DNASequence dna = new DNASequence("atgcCt");
		RNASequence rna = dna.getRNASequence(e);
		Sequence<AminoAcidCompound> peptide = rna.getProteinSequence(e);
		assertThat("Checking lower casing is respected", peptide.getSequenceAsString(),
				is("MP"));
	}

	@Test
	public void translateBrca2ExonOne() throws CompoundNotFoundException {
		TranscriptionEngine e = TranscriptionEngine.getDefault();
		DNASequence dna = new DNASequence(
				"ATGCCTATTGGATCCAAAGAGAGGCCAACATTTTTTGAAATTTTTAAGACACGCTGCAACAAAGCA");
		RNASequence rna = dna.getRNASequence(e);
		Sequence<AminoAcidCompound> peptide = rna.getProteinSequence(e);
		assertThat("Initator methionine wrong", peptide.getSequenceAsString(),
				is("MPIGSKERPTFFEIFKTRCNKA"));
	}

	@Test(timeout=2000)
	public void translateBrca2() {
		TranscriptionEngine e =
				new TranscriptionEngine.Builder().decorateRna(true).build();
		for (int i = 0; i < 100; i++) {
			RNASequence rna = brca2Dna.getRNASequence(e);
			ProteinSequence protein = rna.getProteinSequence(e);
			assertThat("BRCA2 does not translate", protein.getSequenceAsString(),
					is(brca2Pep.getSequenceAsString()));
		}
	}

	@Test
	public void translateInternalStops() {
		TranscriptionEngine e = TranscriptionEngine.getDefault();
		Sequence<AminoAcidCompound> pep = e.translate(volvoxDna);
		assertThat("Ensure internal stops stay", pep.toString(), is(volvoxPep.toString()));
	}

	@Test
	public void translateStopAtInternalStops(){
		//This should stop translation at the first stop codon encountered
		TranscriptionEngine e = new TranscriptionEngine.Builder().stopAtStopCodons(true).build();
		RNASequence rna = ((DNASequence)volvoxDna).getRNASequence();
		String pep = rna.getProteinSequence(e).getSequenceAsString();
		String testpep = volvoxPep.getSequenceAsString().split("\\*")[0];
		assertThat("Translation stops at Stop", pep, is(testpep));
	}

	@Test
	public void waitForStartCodon() throws CompoundNotFoundException{
		//Should not start translation until a start codon is encountered
		TranscriptionEngine e = new TranscriptionEngine.Builder().waitForStartCodon(true).build();
		RNASequence rna = new RNASequence("UCCAUGAGC");
		String pep = rna.getProteinSequence(e).getSequenceAsString();
		assertThat("Translation starts at Start Codon",pep, is("MS"));

		//And should start at start of sequence (NB this is implied by success of all other tests)
		e = new TranscriptionEngine.Builder().waitForStartCodon(false).build();
		pep = rna.getProteinSequence(e).getSequenceAsString();
		assertThat("Translation starts at start of sequence", pep, is("SMS"));
	}

	@Test
	public void translateInitMet() throws CompoundNotFoundException {
		TranscriptionEngine e = TranscriptionEngine.getDefault();
		assertThat("Leucene (CTA) is not changed to init met", e.translate(new DNASequence("CTA")).toString(), is("L"));
		assertThat("Leucene (CTG) is changed to init met", e.translate(new DNASequence("CTG")).toString(), is("M"));

		e = new TranscriptionEngine.Builder().initMet(false).build();
		assertThat("Leucene (CTG) is not changed to init met", e.translate(new DNASequence("CTG")).toString(), is("L"));
	}

	/** test for https://github.com/biojava/biojava/issues/53  */
	@Test
	public void testHashCollision() throws CompoundNotFoundException{
		Builder builder = new TranscriptionEngine.Builder();
		builder.initMet(false);
		builder.translateNCodons(true);
		builder.trimStop(false);
		TranscriptionEngine engine = builder.build();
		Sequence<AminoAcidCompound> seq=engine.translate(new
				DNASequence("GTNTGTTAGTGT"));
		assertThat("XC*C", is(seq.toString()));
		Sequence<AminoAcidCompound> seq2=engine.translate(new
				DNASequence("ANAANG"));
		assertEquals("XX",seq2.toString());
		assertNotSame("HR",seq2.toString());
	}
}

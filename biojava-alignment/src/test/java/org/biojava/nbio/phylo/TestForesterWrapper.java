package org.biojava.nbio.phylo;

import java.io.ByteArrayOutputStream;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.LinkedHashMap;

import org.biojava.nbio.core.sequence.MultipleSequenceAlignment;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.io.FastaReader;
import org.biojava.nbio.core.sequence.io.FastaWriter;
import org.biojava.nbio.core.sequence.io.GenericFastaHeaderParser;
import org.biojava.nbio.core.sequence.io.ProteinSequenceCreator;
import org.biojava.nbio.core.sequence.io.template.FastaHeaderFormatInterface;
import org.forester.msa.Msa;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 * Test the BioJava-forester wrapper methods.
 * 
 * @author Aleix Lafita
 *
 */
public class TestForesterWrapper {

	@Test
	public void testMSAconversion() throws Exception {

		// Load the msa FASTA file into a BioJava MSA object
		InputStream inStream = TestForesterWrapper.class
				.getResourceAsStream("/1u6d_symm.fasta");

		FastaReader<ProteinSequence, AminoAcidCompound> fastaReader = 
				new FastaReader<ProteinSequence, AminoAcidCompound>(
				inStream,
				new GenericFastaHeaderParser<ProteinSequence, AminoAcidCompound>(),
				new ProteinSequenceCreator(AminoAcidCompoundSet
						.getAminoAcidCompoundSet()));

		LinkedHashMap<String, ProteinSequence> proteinSequences = fastaReader
				.process();

		inStream.close();

		MultipleSequenceAlignment<ProteinSequence, AminoAcidCompound> msa = 
				new MultipleSequenceAlignment<ProteinSequence, AminoAcidCompound>();

		String expected = "";
		for (ProteinSequence proteinSequence : proteinSequences.values()) {
			msa.addAlignedSequence(proteinSequence);
			expected += ">" + proteinSequence.getOriginalHeader() + "\n"
					+ proteinSequence.toString() + "\n";
		}

		// Convert the biojava MSA to a FASTA String
		OutputStream os = new ByteArrayOutputStream();
		FastaWriter<ProteinSequence, AminoAcidCompound> fastaW = 
				new FastaWriter<ProteinSequence, AminoAcidCompound>(os,
				msa.getAlignedSequences(),
				new FastaHeaderFormatInterface<ProteinSequence, AminoAcidCompound>() {
					@Override
					public String getHeader(ProteinSequence sequence) {
						return sequence.getAccession().toString();
					};
				});
		fastaW.process();
		String biojava = os.toString();

		// Convert the biojava MSA to a forester Msa
		Msa fMsa = ForesterWrapper.convert(msa);

		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < fMsa.getNumberOfSequences(); i++) {
			sb.append(">" + fMsa.getIdentifier(i) + "\n");
			sb.append(fMsa.getSequenceAsString(i) + "\n");
		}
		String forester = sb.toString();

		// Assert that all FASTA files are equal
		assertEquals(expected, biojava);
		assertEquals(expected, forester);

	}
}

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
/**
 *
 */
package org.biojava.nbio.core.sequence.io;

import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.io.template.GenbankHeaderFormatInterface;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.Sequence;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Collection;

/**
 * The class that should be used to write out genbank file of a sequence
 * collection
 *
 * @author mckeee1
 *
 */
public class GenbankWriterHelper {
	public static final String LINEAR_DNA = "linear";
	public static final String CIRCULAR_DNA = "circular";

	/**
	 * Write collection of protein sequences to a file
	 *
	 * @param file
	 * @param proteinSequences
	 * @throws Exception
	 */
	public static void writeProteinSequence(File file,
			Collection<ProteinSequence> proteinSequences) throws Exception {
		FileOutputStream outputStream = new FileOutputStream(file);
		BufferedOutputStream bo = new BufferedOutputStream(outputStream);
		writeProteinSequence(bo, proteinSequences);
		bo.close();
		outputStream.close();
	}

	/**
	 * Write collection of protein sequences to a stream
	 *
	 * @param outputStream
	 * @param proteinSequences
	 * @throws Exception
	 */

	public static void writeProteinSequence(OutputStream outputStream,
			Collection<ProteinSequence> proteinSequences) throws Exception {

		GenbankWriter<ProteinSequence, AminoAcidCompound> genbankWriter = new GenbankWriter<ProteinSequence, AminoAcidCompound>(
				outputStream,
				proteinSequences,
				new GenericGenbankHeaderFormat<ProteinSequence, AminoAcidCompound>());
		genbankWriter.process();

	}

	/**
	 * Write a collection of NucleotideSequences to a file
	 *
	 * @param file
	 * @param dnaSequences
	 * @throws Exception
	 */

	public static void writeNucleotideSequence(File file,
			Collection<DNASequence> dnaSequences) throws Exception {
		FileOutputStream outputStream = new FileOutputStream(file);
		BufferedOutputStream bo = new BufferedOutputStream(outputStream);
		writeNucleotideSequence(bo, dnaSequences);
		bo.close();
		outputStream.close();
	}

	/**
	 * Write a collection of NucleotideSequences to a file
	 *
	 * @param outputStream
	 * @param dnaSequences
	 * @throws Exception
	 */

	public static void writeNucleotideSequence(OutputStream outputStream,
			Collection<DNASequence> dnaSequences) throws Exception {
		writeNucleotideSequence(outputStream, dnaSequences, LINEAR_DNA);
	}

	/**
	 * Write a collection of NucleotideSequences to a file
	 *
	 * @param outputStream
	 * @param dnaSequences
	 * @param seqType
	 * @throws Exception
	 */

	public static void writeNucleotideSequence(OutputStream outputStream,
			Collection<DNASequence> dnaSequences, String seqType)
			throws Exception {
		GenericGenbankHeaderFormat<DNASequence, NucleotideCompound> genericGenbankHeaderFormat = new GenericGenbankHeaderFormat<DNASequence, NucleotideCompound>(
				seqType);
		// genericGenbankHeaderFormat.setLineSeparator(lineSep);
		GenbankWriter<DNASequence, NucleotideCompound> genbankWriter = new GenbankWriter<DNASequence, NucleotideCompound>(
				outputStream, dnaSequences, genericGenbankHeaderFormat);
		// genbankWriter.setLineSeparator(lineSep);
		genbankWriter.process();
	}

	/**
	 * Write a sequence to a file
	 *
	 * @param file
	 * @param sequence
	 * @throws Exception
	 */
	public static void writeSequence(File file, Sequence<?> sequence)
			throws Exception {
		FileOutputStream outputStream = new FileOutputStream(file);
		BufferedOutputStream bo = new BufferedOutputStream(outputStream);
		writeSequences(bo, singleSeqToCollection(sequence));
		bo.close();
		outputStream.close();
	}

	/**
	 * Write a sequence to OutputStream
	 *
	 * @param outputStream
	 * @param sequence
	 * @throws Exception
	 */
	public static void writeSequence(OutputStream outputStream,
			Sequence<?> sequence) throws Exception {
		writeSequences(outputStream, singleSeqToCollection(sequence));
	}

	/**
	 *
	 * @param sequence
	 * @return
	 */

	private static Collection<Sequence<?>> singleSeqToCollection(
			Sequence<?> sequence) {
		Collection<Sequence<?>> sequences = new ArrayList<Sequence<?>>();
		sequences.add(sequence);
		return sequences;
	}

	/**
	 * Method which will write your given Sequences to the specified
	 * {@link OutputStream}. This is a very generic method which writes just the
	 * AccessionID of the Sequence as the FASTA header.
	 *
	 * @param outputStream
	 *            Stream to write to; can be System.out
	 * @param sequences
	 *            The sequences to write out
	 * @throws Exception
	 *             Thrown normally thanks to IO problems
	 */
	public static void writeSequences(OutputStream outputStream,
			Collection<Sequence<?>> sequences) throws Exception {

		GenbankHeaderFormatInterface<Sequence<?>, Compound> fhfi = new GenbankHeaderFormatInterface<Sequence<?>, Compound>() {

			@Override
			public String getHeader(Sequence<?> sequence) {
				return sequence.getAccession().toString();
			}

			;
		};

		GenbankWriter<Sequence<?>, Compound> genbankWriter = new GenbankWriter<Sequence<?>, Compound>(
				outputStream, sequences, fhfi);

		genbankWriter.process();
	}
}

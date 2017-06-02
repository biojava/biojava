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
 * Created on 01-21-2010
 */
package org.biojava.nbio.core.sequence.io;

import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.RNASequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.compound.DNACompoundSet;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.compound.RNACompoundSet;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.LinkedHashMap;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class FastaReaderHelper {

	/**
	 * Selecting lazySequenceLoad=true will parse the FASTA file and figure out the accessionid and offsets and return sequence objects
	 * that can in the future read the sequence from the disk. This allows the loading of large fasta files where you are only interested
	 * in one sequence based on accession id.
	 * @param file
	 * @param lazySequenceLoad
	 * @return
	 * @throws IOException
	 */
	public static LinkedHashMap<String, DNASequence> readFastaDNASequence(File file, boolean lazySequenceLoad) throws IOException {
		if (!lazySequenceLoad) {
			return readFastaDNASequence(file);
		}

		FastaReader<DNASequence, NucleotideCompound> fastaProxyReader =
				new FastaReader<DNASequence, NucleotideCompound>(
						file,
						new GenericFastaHeaderParser<DNASequence, NucleotideCompound>(),
						new FileProxyDNASequenceCreator(
								file,
								DNACompoundSet.getDNACompoundSet(),
								new FastaSequenceParser()
							)
					);
		return fastaProxyReader.process();

	}

	/**
	 * Selecting lazySequenceLoad=true will parse the FASTA file and figure out the accessionid and offsets and return sequence objects
	 * that can in the future read the sequence from the disk. This allows the loading of large fasta files where you are only interested
	 * in one sequence based on accession id.
	 * @param file
	 * @param lazySequenceLoad
	 * @return
	 * @throws IOException
	 */
	public static LinkedHashMap<String, RNASequence> readFastaRNASequence(File file, boolean lazySequenceLoad) throws IOException {
		if (!lazySequenceLoad) {
			return readFastaRNASequence(file);
		}

		FastaReader<RNASequence, NucleotideCompound> fastaProxyReader =
				new FastaReader<RNASequence, NucleotideCompound>(
						file,
						new GenericFastaHeaderParser<RNASequence, NucleotideCompound>(),
						new FileProxyRNASequenceCreator(
								file,
								RNACompoundSet.getRNACompoundSet(),
								new FastaSequenceParser()
							)
					);
		return fastaProxyReader.process();

	}

	/**
	 * Read a fasta file containing amino acids with setup that would handle most
	 * cases.
	 *
	 * @param file
	 * @return
	 * @throws IOException
	 */
	public static LinkedHashMap<String, ProteinSequence> readFastaProteinSequence(
			File file) throws IOException {
		FileInputStream inStream = new FileInputStream(file);
		LinkedHashMap<String, ProteinSequence> proteinSequences = readFastaProteinSequence(inStream);
		inStream.close();
		return proteinSequences;
	}

	/**
	 * Read a fasta file containing amino acids with setup that would handle most
	 * cases. User is responsible for closing InputStream because you opened it
	 *
	 * @param inStream
	 * @return
	 * @throws IOException
	 */
	public static LinkedHashMap<String, ProteinSequence> readFastaProteinSequence(
			InputStream inStream) throws IOException {
		FastaReader<ProteinSequence, AminoAcidCompound> fastaReader = new FastaReader<ProteinSequence, AminoAcidCompound>(
				inStream,
				new GenericFastaHeaderParser<ProteinSequence, AminoAcidCompound>(),
				new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet()));
		return fastaReader.process();
	}

	/**
	 * Read a fasta DNA sequence
	 * @param inStream
	 * @return
	 * @throws IOException
	 */
	public static LinkedHashMap<String, DNASequence> readFastaDNASequence(
			InputStream inStream) throws IOException {
		FastaReader<DNASequence, NucleotideCompound> fastaReader = new FastaReader<DNASequence, NucleotideCompound>(
				inStream,
				new GenericFastaHeaderParser<DNASequence, NucleotideCompound>(),
				new DNASequenceCreator(DNACompoundSet.getDNACompoundSet()));
		return fastaReader.process();
	}

	/**
	 *
	 * @param file
	 * @return
	 * @throws IOException
	 */
	public static LinkedHashMap<String, DNASequence> readFastaDNASequence(
			File file) throws IOException {
		FileInputStream inStream = new FileInputStream(file);
		LinkedHashMap<String, DNASequence> dnaSequences = readFastaDNASequence(inStream);
		inStream.close();
		return dnaSequences;
	}

	/**
	 * Read a fasta RNA sequence
	 * @param inStream
	 * @return
	 * @throws IOException
	 */
	public static LinkedHashMap<String, RNASequence> readFastaRNASequence(
			InputStream inStream) throws IOException {
		FastaReader<RNASequence, NucleotideCompound> fastaReader = new FastaReader<RNASequence, NucleotideCompound>(
				inStream,
				new GenericFastaHeaderParser<RNASequence, NucleotideCompound>(),
				new RNASequenceCreator(RNACompoundSet.getRNACompoundSet()));
		return fastaReader.process();
	}

	/**
	 *
	 * @param file
	 * @return
	 * @throws IOException
	 */
	public static LinkedHashMap<String, RNASequence> readFastaRNASequence(
			File file) throws IOException {
		FileInputStream inStream = new FileInputStream(file);
		LinkedHashMap<String, RNASequence> rnaSequences = readFastaRNASequence(inStream);
		inStream.close();
		return rnaSequences;
	}

	public static void main(String[] args) throws Exception {

		LinkedHashMap<String, DNASequence> dnaSequences = FastaReaderHelper.readFastaDNASequence(new File("fasta.fna"));
		for (DNASequence sequence : dnaSequences.values()) {
			sequence.getRNASequence().getProteinSequence().getSequenceAsString();
		}
	}
}

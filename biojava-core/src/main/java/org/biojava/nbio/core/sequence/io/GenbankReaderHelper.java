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
import org.biojava.nbio.core.sequence.template.AbstractSequence;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.util.LinkedHashMap;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class GenbankReaderHelper {

	private final static Logger logger = LoggerFactory.getLogger(GenbankReaderHelper.class);

	/**
	 * Selecting lazySequenceLoad=true will parse the Genbank file and figure out the accessionid and offsets and return sequence objects
	 * that can in the future read the sequence from the disk. This allows the loading of large Genbank files where you are only interested
	 * in one sequence based on accession id.
	 * @param file
	 * @param lazySequenceLoad
	 * @return
	 * @throws Exception
	 */
	public static LinkedHashMap<String, DNASequence> readGenbankDNASequence(File file, boolean lazySequenceLoad) throws Exception {
		if (!lazySequenceLoad) {
			return readGenbankDNASequence(file);
		}

		GenbankReader<DNASequence, NucleotideCompound> GenbankProxyReader =
				new GenbankReader<DNASequence, NucleotideCompound>(
						file,
						new GenericGenbankHeaderParser<DNASequence, NucleotideCompound>(),
						new FileProxyDNASequenceCreator(
								file,
								DNACompoundSet.getDNACompoundSet(),
								new GenbankSequenceParser<AbstractSequence<NucleotideCompound>, NucleotideCompound>()
							)
					);
		return GenbankProxyReader.process();

	}

	/**
	 * Selecting lazySequenceLoad=true will parse the Genbank file and figure out the accessionid and offsets and return sequence objects
	 * that can in the future read the sequence from the disk. This allows the loading of large Genbank files where you are only interested
	 * in one sequence based on accession id.
	 * @param file
	 * @param lazySequenceLoad
	 * @return
	 * @throws Exception
	 */
	public static LinkedHashMap<String, ProteinSequence> readGenbankProteinSequence(File file, boolean lazySequenceLoad) throws Exception {
		if (!lazySequenceLoad) {
			return readGenbankProteinSequence(file);
		}

		GenbankReader<ProteinSequence, AminoAcidCompound> GenbankProxyReader =
				new GenbankReader<ProteinSequence, AminoAcidCompound>(
						file,
						new GenericGenbankHeaderParser<ProteinSequence, AminoAcidCompound>(),
						new FileProxyProteinSequenceCreator(
								file,
								AminoAcidCompoundSet.getAminoAcidCompoundSet(),
								new GenbankSequenceParser<AbstractSequence<AminoAcidCompound>, AminoAcidCompound>()
							)
					);
		return GenbankProxyReader.process();

	}

	/**
	 * Selecting lazySequenceLoad=true will parse the Genbank file and figure out the accessionid and offsets and return sequence objects
	 * that can in the future read the sequence from the disk. This allows the loading of large Genbank files where you are only interested
	 * in one sequence based on accession id.
	 * @param file
	 * @param lazySequenceLoad
	 * @return
	 * @throws Exception
	 */
	public static LinkedHashMap<String, RNASequence> readGenbankRNASequence(File file, boolean lazySequenceLoad) throws Exception {
		if (!lazySequenceLoad) {
			return readGenbankRNASequence(file);
		}

		GenbankReader<RNASequence, NucleotideCompound> GenbankProxyReader =
				new GenbankReader<RNASequence, NucleotideCompound>(
						file,
						new GenericGenbankHeaderParser<RNASequence, NucleotideCompound>(),
						new FileProxyRNASequenceCreator(
								file,
								RNACompoundSet.getRNACompoundSet(),
								new GenbankSequenceParser<AbstractSequence<NucleotideCompound>, NucleotideCompound>()
							)
					);
		return GenbankProxyReader.process();

	}

	/**
	 * Read a Genbank file containing amino acids with setup that would handle most
	 * cases.
	 *
	 * @param file
	 * @return
	 * @throws Exception
	 */
	public static LinkedHashMap<String, ProteinSequence> readGenbankProteinSequence(
			File file) throws Exception {
		FileInputStream inStream = new FileInputStream(file);
		LinkedHashMap<String, ProteinSequence> proteinSequences = readGenbankProteinSequence(inStream);
		inStream.close();
		return proteinSequences;
	}

	/**
	 * Read a Genbank file containing amino acids with setup that would handle most
	 * cases. User is responsible for closing InputStream because you opened it
	 *
	 * @param inStream
	 * @return
	 * @throws Exception
	 */
	public static LinkedHashMap<String, ProteinSequence> readGenbankProteinSequence(
			InputStream inStream) throws Exception {
		GenbankReader<ProteinSequence, AminoAcidCompound> GenbankReader = new GenbankReader<ProteinSequence, AminoAcidCompound>(
				inStream,
				new GenericGenbankHeaderParser<ProteinSequence, AminoAcidCompound>(),
				new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet()));
		return GenbankReader.process();
	}

	/**
	 * Read a Genbank DNA sequence
	 * @param inStream
	 * @return
	 * @throws Exception
	 */
	public static LinkedHashMap<String, DNASequence> readGenbankDNASequence(
			InputStream inStream) throws Exception {
		GenbankReader<DNASequence, NucleotideCompound> GenbankReader = new GenbankReader<DNASequence, NucleotideCompound>(
				inStream,
				new GenericGenbankHeaderParser<DNASequence, NucleotideCompound>(),
				new DNASequenceCreator(DNACompoundSet.getDNACompoundSet()));
		return GenbankReader.process();
	}

	/**
	 *
	 * @param file
	 * @return
	 * @throws Exception
	 */
	public static LinkedHashMap<String, DNASequence> readGenbankDNASequence(
			File file) throws Exception {
		FileInputStream inStream = new FileInputStream(file);
		LinkedHashMap<String, DNASequence> dnaSequences = readGenbankDNASequence(inStream);
		inStream.close();
		return dnaSequences;
	}
	/**
	 * Read a Genbank RNA sequence
	 * @param inStream
	 * @return
	 * @throws Exception
	 */
	public static LinkedHashMap<String, RNASequence> readGenbankRNASequence(
			InputStream inStream) throws Exception {
		GenbankReader<RNASequence, NucleotideCompound> GenbankReader = new GenbankReader<RNASequence, NucleotideCompound>(
				inStream,
				new GenericGenbankHeaderParser<RNASequence, NucleotideCompound>(),
				new RNASequenceCreator(RNACompoundSet.getRNACompoundSet()));
		return GenbankReader.process();
	}

	/**
	 *
	 * @param file
	 * @return
	 * @throws Exception
	 */
	public static LinkedHashMap<String, RNASequence> readGenbankRNASequence(
			File file) throws Exception {
		FileInputStream inStream = new FileInputStream(file);
		LinkedHashMap<String, RNASequence> rnaSequences = readGenbankRNASequence(inStream);
		inStream.close();
		return rnaSequences;
	}

	public static void main(String[] args) throws Exception {

		LinkedHashMap<String, DNASequence> dnaSequences = GenbankReaderHelper.readGenbankDNASequence(new File("src/test/resources/NM_000266.gb"), true);
		for (DNASequence sequence : dnaSequences.values()) {
			logger.info("DNA Sequence: {}", sequence.getRNASequence().getProteinSequence().getSequenceAsString());
		}

		LinkedHashMap<String, ProteinSequence> proteinSequences = GenbankReaderHelper.readGenbankProteinSequence(new File("src/test/resources/BondFeature.gb"), true);
		for (ProteinSequence sequence : proteinSequences.values()) {
			logger.info("Protein Sequence: {}", sequence.getSequenceAsString());
		}
	}
}

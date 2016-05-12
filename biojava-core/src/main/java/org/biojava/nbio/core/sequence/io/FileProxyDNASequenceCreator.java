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

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.io.template.SequenceCreatorInterface;
import org.biojava.nbio.core.sequence.io.template.SequenceParserInterface;
import org.biojava.nbio.core.sequence.loader.SequenceFileProxyLoader;
import org.biojava.nbio.core.sequence.template.AbstractSequence;
import org.biojava.nbio.core.sequence.template.CompoundSet;
import org.biojava.nbio.core.sequence.template.ProxySequenceReader;

import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 * This class is a good example of using the SequenceCreatorInterface where during parsing of the stream
 * the sequence and the offset index are passed to create a Protein sequence that will be loaded in lazily.
 * This way you can load very large fasta files and store accession id and delay loading the sequence to save
 * memory. The index is the file stream offset so when a DNASequence has a call to getSequence() the
 * SequenceFileProxyLoader will open the file and offset to the index and retrieve the sequence.
 *
 * Same approach can be used for genome sequence data stored in a local fasta file, in a database or via http
 * interface to a remote server
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class FileProxyDNASequenceCreator implements
		SequenceCreatorInterface<NucleotideCompound> {

	CompoundSet<NucleotideCompound> compoundSet = null;
	File file = null;
	SequenceParserInterface sequenceParser;

	/**
	 * Need File so that we can store full path name in SequenceFileProxyLoader for Random File access as a quick read
	 * @param fastaFile
	 * @param compoundSet
	 */
	public FileProxyDNASequenceCreator(File file,
			CompoundSet<NucleotideCompound> compoundSet,
			SequenceParserInterface sequenceParser) {
		this.compoundSet = compoundSet;
		this.file = file;
		this.sequenceParser = sequenceParser;
	}

	/**
	 * Even though we are passing in the sequence we really only care about the length of the sequence and the offset
	 * index in the fasta file.
	 * @param sequence
	 * @param index
	 * @return
	 * @throws CompoundNotFoundException
	 * @throws IOException
	 */
	@Override
	public AbstractSequence<NucleotideCompound> getSequence(String sequence, long index ) throws CompoundNotFoundException, IOException {
		SequenceFileProxyLoader<NucleotideCompound> sequenceFileProxyLoader = new SequenceFileProxyLoader<NucleotideCompound>(
				file,
				sequenceParser,
				index,
				sequence.length(),
				compoundSet);
		return new DNASequence(sequenceFileProxyLoader, compoundSet);
	}

	/**
	 * Should be able to extend the same concept to a remote URL call or database connection. Not supported yet
	 * @param proxyLoader
	 * @param index
	 * @return
	 */
	@Override
	public AbstractSequence<NucleotideCompound> getSequence(
			ProxySequenceReader<NucleotideCompound> proxyLoader, long index) {
		throw new UnsupportedOperationException("Not supported yet.");
	}

	/**
	 * Not sure of use case and currently not supported
	 * @param list
	 * @return
	 */
	@Override
	public AbstractSequence<NucleotideCompound> getSequence(
			List<NucleotideCompound> list) {
		throw new UnsupportedOperationException("Not supported yet.");
	}
}

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
package org.biojava3.core.sequence.io;

import java.io.File;
import java.util.List;

import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.io.template.SequenceCreatorInterface;
import org.biojava3.core.sequence.loader.SequenceFileProxyLoader;
import org.biojava3.core.sequence.template.AbstractSequence;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.ProxySequenceReader;

/**
 * This class is a good example of using the SequenceCreatorInterface where during parsing of the stream
 * the sequence and the offset index are passed to create a Protein sequence that will be loaded in lazily.
 * This way you can load very large fasta files and store accession id and delay loading the sequence to save
 * memory. The index is the file stream offset so when a ProteinSequence has a call to getSequence() the
 * SequenceFileProxyLoader will open the file and offset to the index and retrieve the sequence.
 *
 * Same approach can be used for genome sequence data stored in a local fasta file, in a database or via http
 * interface to a remote server
 *
 * @author Scooter Willis &lt;willishf at gmail dot com&gt;
 */
public class FileProxyProteinSequenceCreator implements
        SequenceCreatorInterface<AminoAcidCompound> {

    CompoundSet<AminoAcidCompound> compoundSet = null;
    File fastaFile = null;

    /**
     * Need File so that we can store full path name in SequenceFileProxyLoader for Random File access as a quick read
     * @param fastaFile
     * @param compoundSet
     */
    public FileProxyProteinSequenceCreator(File fastaFile,
            CompoundSet<AminoAcidCompound> compoundSet) {
        this.compoundSet = compoundSet;
        this.fastaFile = fastaFile;
    }

    /**
     * Even though we are passing in the sequence we really only care about the length of the sequence and the offset
     * index in the fasta file.
     * @param sequence
     * @param index
     * @return
     */

    public AbstractSequence<AminoAcidCompound> getSequence(String sequence,
            long index) {
        SequenceFileProxyLoader<AminoAcidCompound> sequenceFileProxyLoader = new SequenceFileProxyLoader<AminoAcidCompound>(
                fastaFile, new FastaSequenceParser(), index, sequence.length(),
                compoundSet);
        return new ProteinSequence(sequenceFileProxyLoader, compoundSet);
    }

    /**
     * Should be able to extend the same concept to a remote URL call or database connection. Not supported yet
     * @param proxyLoader
     * @param index
     * @return
     */
    public AbstractSequence<AminoAcidCompound> getSequence(
            ProxySequenceReader<AminoAcidCompound> proxyLoader, long index) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    /**
     * Not sure of use case and currently not supported
     * @param list
     * @return
     */
    public AbstractSequence<AminoAcidCompound> getSequence(
            List<AminoAcidCompound> list) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}

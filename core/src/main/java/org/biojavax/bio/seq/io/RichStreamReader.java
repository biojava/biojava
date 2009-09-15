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

package org.biojavax.bio.seq.io;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.NoSuchElementException;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojavax.Namespace;
import org.biojavax.bio.BioEntry;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.bio.seq.RichSequenceIterator;

/**
 * Parses a stream into sequences.
 * This object implements SequenceIterator, so you can loop over each sequence
 * produced. It consumes a stream, and uses a SequenceFormat to extract each
 * sequence from the stream.
 * It is assumed that the stream contains sequences that can be handled by the
 * one format, and that they are not seperated other than by delimiters that the
 * format can handle.
 * Sequences are instantiated when they are requested by nextSequence, not
 * before, so it is safe to use this object to parse a gigabyte fasta file, and
 * do sequence-by-sequence processing, while being guaranteed that RichStreamReader
 * will not require you to keep any of the sequences in memory.
 * @author Matthew Pocock
 * @author Thomas Down
 * @author Richard Holland
 * @since 1.5
 */

public class RichStreamReader implements RichSequenceIterator {
    
    /**
     * The symbol parser.
     */
    private Namespace ns;
    
    /**
     * The symbol parser.
     */
    private SymbolTokenization symParser;
    
    /**
     * The sequence format.
     */
    private RichSequenceFormat format;
    
    /**
     * The sequence-builder factory.
     */
    private RichSequenceBuilderFactory sf;
    
    /**
     * The stream of data to parse.
     */
    
    private BufferedReader reader;
    
    /**
     * Flag indicating if more sequences are available.
     */
    private boolean moreSequenceAvailable = true;
    
    /**
     * {@inheritDoc}
     */
    public Sequence nextSequence() throws NoSuchElementException, BioException {
        return this.nextRichSequence();
    }
    
    /**
     * {@inheritDoc}
     */
    public BioEntry nextBioEntry() throws NoSuchElementException, BioException {
        return this.nextRichSequence();
    }
    
    /**
     * {@inheritDoc}
     */
    public RichSequence nextRichSequence() throws NoSuchElementException, BioException {
        if(!moreSequenceAvailable)
            throw new NoSuchElementException("Stream is empty");
        try {
            RichSequenceBuilder builder = (RichSequenceBuilder)sf.makeSequenceBuilder();
            moreSequenceAvailable = format.readRichSequence(reader, symParser, builder, ns);
            return builder.makeRichSequence();
        } catch (Exception e) {
            throw new BioException("Could not read sequence",e);
        }
    }
    
    /**
     * {@inheritDoc}
     */
    public boolean hasNext() {
        return moreSequenceAvailable;
    }
    
    /**
     * Creates a new stream reader on the given input stream, which will attempt to read
     * sequences in the given format, having symbols from the given tokenization, and
     * pass them to the given factory to be transformed into RichSequence objects in
     * the given namespace.
     * @param is the input stream to read from
     * @param format the input file format
     * @param symParser the tokenizer that understands the sequence symbols in the file
     * @param sf the factory that will build the sequences
     * @param ns the namespace the sequences will be loaded into.
     */
    public RichStreamReader(InputStream is,
            RichSequenceFormat format,
            SymbolTokenization symParser,
            RichSequenceBuilderFactory sf,
            Namespace ns)  {
        this(new BufferedReader(new InputStreamReader(is)), format,symParser,sf,ns);
    }
    
    /**
     * Creates a new stream reader on the given reader, which will attempt to read
     * sequences in the given format, having symbols from the given tokenization, and
     * pass them to the given factory to be transformed into RichSequence objects in
     * the given namespace.
     * @param reader the reader to read from
     * @param format the input file format
     * @param symParser the tokenizer that understands the sequence symbols in the file
     * @param sf the factory that will build the sequences
     * @param ns the namespace the sequences will be loaded into.
     */
    public RichStreamReader(BufferedReader reader,
            RichSequenceFormat format,
            SymbolTokenization symParser,
            RichSequenceBuilderFactory sf,
            Namespace ns)  {
        this.reader = reader;
        this.format = format;
        this.symParser = symParser;
        this.sf = sf;
        this.ns = ns;
    }
}

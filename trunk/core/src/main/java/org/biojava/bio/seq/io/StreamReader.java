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


package org.biojava.bio.seq.io;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.NoSuchElementException;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;

/**
 * Parses a stream into sequences.
 * <p>
 * This object implements SequenceIterator, so you can loop over each sequence
 * produced. It consumes a stream, and uses a SequenceFormat to extract each
 * sequence from the stream.
 * <p>
 * It is assumed that the stream contains sequences that can be handled by the
 * one format, and that they are not seperated other than by delimiters that the
 * format can handle.
 * <p>
 * Sequences are instantiated when they are requested by nextSequence, not
 * before, so it is safe to use this object to parse a gigabyte fasta file, and
 * do sequence-by-sequence processing, while being guaranteed that StreamReader
 * will not require you to keep any of the sequences in memory.
 * 
 * <p>More functionality is offered by {@link org.biojavax.bio.seq.io.RichStreamReader RichStreamReader},
 * Use of this interface is prefered.</p>
 *
 * @author Matthew Pocock
 * @author Thomas Down
 * @see org.biojavax.bio.seq.io.RichStreamReader
 */

public class StreamReader implements SequenceIterator, org.biojava.utils.ParseErrorListener
{
    /**
     * The symbol parser.
     */
    private SymbolTokenization symParser;

    /**
     * The sequence format.
     */
    private SequenceFormat format;

    /**
     * The sequence-builder factory.
     */
    private SequenceBuilderFactory sf;

    /**
     * The stream of data to parse.
     */

    private BufferedReader reader;

    /**
     * Flag indicating if more sequences are available.
     */

    private boolean moreSequenceAvailable = true;

    /**
     * Pull the next sequence out of the stream.
     * <p>
     * This method will delegate parsing from the stream to a SequenceFormat
     * object, and then return the resulting sequence.
     *
     * @return the next Sequence
     * @throws NoSuchElementException if the end of the stream has been hit
     * @throws BioException if for any reason the next sequence could not be read
     */

    public Sequence nextSequence()
        throws NoSuchElementException, BioException
    {
        if(!moreSequenceAvailable)
            throw new NoSuchElementException("Stream is empty");
        try {
            SequenceBuilder builder = sf.makeSequenceBuilder();
            moreSequenceAvailable = format.readSequence(reader, symParser, builder);
            return builder.makeSequence();
        } catch (Exception e) {
            throw new BioException("Could not read sequence",e);
        }
    }

    public boolean hasNext() {
        return moreSequenceAvailable;
    }

    public StreamReader(InputStream is,
                        SequenceFormat format,
                        SymbolTokenization symParser,
                        SequenceBuilderFactory sf)  {
        this.reader = new BufferedReader(new InputStreamReader(is));
        this.format = format;
        this.symParser = symParser;
        this.sf = sf;
    }

    public StreamReader(BufferedReader reader,
                        SequenceFormat format,
                        SymbolTokenization symParser,
                        SequenceBuilderFactory sf)  {
        this.reader = reader;
        this.format = format;
        this.symParser = symParser;
        this.sf = sf;
        ((org.biojava.utils.ParseErrorSource)(this.format)).addParseErrorListener(this);
    }

        /**
         * This method determines the behaviour when a bad line is processed.
         * Some options are to log the error, throw an exception, ignore it
         * completely, or pass the event through.
         * <p>
         * This method should be overwritten when different behavior is desired.
         *
         * @param theEvent The event that contains the bad line and token.
         */
        public void BadLineParsed(org.biojava.utils.ParseErrorEvent theEvent)
        {
                System.err.println(theEvent.getMessage());
        }
}

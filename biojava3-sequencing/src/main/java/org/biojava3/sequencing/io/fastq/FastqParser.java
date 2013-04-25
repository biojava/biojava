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
package org.biojava3.sequencing.io.fastq;

import java.io.Closeable;
import java.io.IOException;

import com.google.common.io.CharStreams;
import com.google.common.io.InputSupplier;
import com.google.common.io.LineProcessor;

/**
 * Low-level event based parser for FASTQ formatted sequences.
 *
 * @since 3.0.3
 */
final class FastqParser
{

    /**
     * Parse the specified input supplier.
     *
     * @param supplier input supplier, must not be null
     * @param listener low-level event based parser callback, must not be null
     * @throws IOException if an I/O error occurs
     */
    static <R extends Readable & Closeable> void parse(final InputSupplier<R> supplier, final ParseListener listener)
        throws IOException
    {
        if (supplier == null)
        {
            throw new IllegalArgumentException("supplier must not be null");
        }
        FastqParserLineProcessor lineProcessor = new FastqParserLineProcessor(listener);
        CharStreams.readLines(supplier, lineProcessor);
        if (lineProcessor.getState() == State.COMPLETE)
        {
            listener.complete();
            lineProcessor.setState(State.DESCRIPTION);
        }
        if (lineProcessor.getState() != State.DESCRIPTION)
        {
            throw new IOException("truncated sequence"); // at line " + lineNumber);
        }
    }

    /**
     * FASTQ formatted sequence parser line processor.
     */
    private static class FastqParserLineProcessor implements LineProcessor<Object>
    {
        /** Parser state. */
        private State state = State.DESCRIPTION;

        /** Sequence length. */
        private int sequenceLength = 0;

        /** Quality length. */
        private int qualityLength = 0;

        /** Parse listener. */
        private final ParseListener listener;


        /**
         * Create a new FASTQ formatted sequence parser line processor with the specified parse listener.
         *
         * @param listener parse listener, must not be null
         */
        private FastqParserLineProcessor(final ParseListener listener)
        {
            if (listener == null)
            {
                throw new IllegalArgumentException("listener must not be null");
            }
            this.listener = listener;
        }


        /**
         * Return the parser state.
         *
         * @return the parser state
         */
        private State getState()
        {
            return state;
        }

        /**
         * Set the parser state to <code>state</code>.
         *
         * @param state parser state
         */
        private void setState(final State state)
        {
            this.state = state;
        }

        @Override
        public Object getResult()
        {
            return null;
        }

        @Override
        public boolean processLine(final String line) throws IOException
        {
            String sequence = null;
            String quality = null;
            switch (state)
            {
            case DESCRIPTION:
                if (line.startsWith("@"))
                {
                    listener.description(line.substring(1).trim());
                    state = State.SEQUENCE;
                }
                else
                {
                    throw new IOException("description must begin with a '@' character");
                }
                break;
            case SEQUENCE:
                sequence = line.trim();
                listener.sequence(sequence);
                sequenceLength = sequence.length();
                state = State.REPEAT_DESCRIPTION;
                break;
            case REPEAT_DESCRIPTION:
                if (line.startsWith("+"))
                {
                    listener.repeatDescription(line.substring(1).trim());
                    state = State.QUALITY;
                }
                else
                {
                    sequence = line.trim();
                    listener.appendSequence(sequence);
                    sequenceLength += sequence.length();
                }
                break;
            case QUALITY:
                quality = line.trim();
                listener.quality(quality);
                qualityLength = quality.length();
                state = State.COMPLETE;
                break;
            case COMPLETE:
                if (sequenceLength == qualityLength)
                {
                    listener.complete();

                    if (line.startsWith("@"))
                    {
                        listener.description(line.substring(1).trim());
                        state = State.SEQUENCE;
                    }
                    else
                    {
                        throw new IOException("description must begin with a '@' character");
                    }
                }
                else
                {
                    quality = line.trim();
                    listener.appendQuality(quality);
                    qualityLength += quality.length();
                }
                break;
            default:
                break;
            }
            return true;
        }
    }

    /** Parser state. */
    private static enum State
    {
        /** Description parser state. */
        DESCRIPTION,

        /** Sequence parser state. */
        SEQUENCE,

        /** Repeat description parser state. */
        REPEAT_DESCRIPTION,

        /** Quality score parser state. */
        QUALITY,

        /** Complete parser state. */
        COMPLETE;
    };
}
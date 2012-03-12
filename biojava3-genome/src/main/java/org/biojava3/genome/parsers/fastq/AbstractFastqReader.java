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
package org.biojava3.genome.parsers.fastq;

import java.io.*;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;

/**
 * Abstract reader implementation for FASTQ formatted sequences.
 *
 * @since 3.0.3
 */
abstract class AbstractFastqReader implements FastqReader {

    /**
     * Parser state.
     */
    private static enum State {

        /**
         * Description parser state.
         */
        DESCRIPTION,
        /**
         * Sequence parser state.
         */
        SEQUENCE,
        /**
         * Repeat description parser state.
         */
        REPEAT_DESCRIPTION,
        /**
         * Quality score parser state.
         */
        QUALITY,
        /**
         * Complete parser state.
         */
        COMPLETE;
    };

    /**
     * Return the FASTQ sequence format variant for this reader.
     *
     * @return the FASTQ sequence format variant for this reader
     */
    protected abstract FastqVariant getVariant();

    /**
     * Validate the specified description.
     *
     * @param builder FASTQ formatted sequence builder, will not be null
     * @param description description to validate, will not be null
     * @param lineNumber current line number in input stream
     * @throws IOException if the specified description is not valid
     */
    protected abstract void validateDescription(FastqBuilder builder, String description, int lineNumber) throws IOException;

    /**
     * Validate the specified sequence.
     *
     * @param builder FASTQ formatted sequence builder, will not be null
     * @param sequence sequence to validate, will not be null
     * @param lineNumber current line number in input stream
     * @throws IOException if the specified sequence is not valid
     */
    protected abstract void validateSequence(FastqBuilder builder, String sequence, int lineNumber) throws IOException;

    /**
     * Validate the specified repeat description.
     *
     * @param builder FASTQ formatted sequence builder, will not be null
     * @param repeatDescription repeat description to validate, will not be null
     * @param lineNumber current line number in input stream
     * @throws IOException if the specified repeat description is not valid
     */
    protected abstract void validateRepeatDescription(FastqBuilder builder, String repeatDescription, int lineNumber) throws IOException;

    /**
     * Validate the specified quality scores.
     *
     * @param builder FASTQ formatted sequence builder, will not be null
     * @param quality quality scores to validate, will not be null
     * @param lineNumber current line number in input stream
     * @throws IOException if the specified quality scores are not valid
     */
    protected abstract void validateQuality(FastqBuilder builder, String quality, int lineNumber) throws IOException;

    /**
     * {@inheritDoc}
     */
    @Override
    public final Iterable<Fastq> read(final File file) throws IOException {
        if (file == null) {
            throw new IllegalArgumentException("file must not be null");
        }
        InputStream inputStream = null;
        try {
            inputStream = new FileInputStream(file);
            return read(inputStream);
        } catch (IOException e) {
            throw e;
        } finally {
            if (inputStream != null) {
                try {
                    inputStream.close();
                } catch (IOException e) {
                    // ignore
                }
            }
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public final Iterable<Fastq> read(final URL url) throws IOException {
        if (url == null) {
            throw new IllegalArgumentException("url must not be null");
        }
        InputStream inputStream = null;
        try {
            inputStream = url.openStream();
            return read(inputStream);
        } catch (IOException e) {
            throw e;
        } finally {
            if (inputStream != null) {
                try {
                    inputStream.close();
                } catch (IOException e) {
                    // ignore
                }
            }
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public final Iterable<Fastq> read(final InputStream inputStream) throws IOException {
        if (inputStream == null) {
            throw new IllegalArgumentException("inputStream must not be null");
        }
        int lineNumber = 0;
        State state = State.DESCRIPTION;
        List<Fastq> result = new ArrayList<Fastq>();
        FastqBuilder builder = new FastqBuilder().withVariant(getVariant());
        BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream));
        while (reader.ready()) {
            String line = reader.readLine();
            switch (state) {
                case DESCRIPTION:
                    validateDescription(builder, line, lineNumber);
                    builder.withDescription(line.substring(1).trim());
                    state = State.SEQUENCE;
                    break;
                case SEQUENCE:
                    validateSequence(builder, line, lineNumber);
                    builder.withSequence(line.trim());
                    state = State.REPEAT_DESCRIPTION;
                    break;
                case REPEAT_DESCRIPTION:
                    if (!line.startsWith("+")) {
                        builder.appendSequence(line.trim());
                    } else {
                        validateRepeatDescription(builder, line, lineNumber);
                        state = State.QUALITY;
                    }
                    break;
                case QUALITY:
                    validateQuality(builder, line, lineNumber);
                    builder.withQuality(line.trim());
                    state = State.COMPLETE;
                    break;
                case COMPLETE:
                    if (!builder.sequenceAndQualityLengthsMatch()) {
                        builder.appendQuality(line.trim());
                    } else {
                        try {
                            result.add(builder.build());
                        } catch (IllegalStateException e) {
                            throw new IOException("caught an IllegalStateException at line " + lineNumber + " " + e.getMessage());
                            //throw new IOException("caught an IllegalStateException at line " + lineNumber, e);  jdk 1.6+
                        }
                        validateDescription(builder, line, lineNumber);
                        builder.withDescription(line.substring(1).trim());
                        state = State.SEQUENCE;
                    }
                    break;
                default:
                    break;
            }
            lineNumber++;
        }
        if (state == State.COMPLETE) {
            try {
                result.add(builder.build());
                state = State.DESCRIPTION;
            } catch (IllegalStateException e) {
                //throw new IOException("caught an IllegalStateException at line " + lineNumber + " " + e.getMessage());
                throw new IOException("caught an IllegalStateException at line " + lineNumber, e);
            }
        }
        if (state != State.DESCRIPTION) {
            throw new IOException("truncated sequence at line " + lineNumber);
        }
        return result;
    }
}
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

import java.io.IOException;

/**
 * Reader for {@link FastqVariant#FASTQ_SANGER} formatted sequences.
 *
 * @since 3.0.3
 */
public final class SangerFastqReader extends AbstractFastqReader {

    /**
     * {@inheritDoc}
     */
    @Override
    protected FastqVariant getVariant() {
        return FastqVariant.FASTQ_SANGER;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateDescription(final FastqBuilder builder, final String description, final int lineNumber) throws IOException {
        if (!description.startsWith("@")) {
            throw new IOException("description must begin with a '@' character at line " + lineNumber);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSequence(final FastqBuilder builder, final String sequence, final int lineNumber) throws IOException {
        // empty
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateRepeatDescription(final FastqBuilder builder, final String repeatDescription, final int lineNumber) throws IOException {
        String description = builder.getDescription();
        if ((description != null) && description.length() > 0 && repeatDescription.length() > 1) {
            if (!description.equals(repeatDescription.substring(1))) {
                throw new IOException("repeat description must match description at line " + lineNumber);
            }
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateQuality(final FastqBuilder builder, final String quality, final int lineNumber) throws IOException {
        for (int i = 0; i < quality.length(); i++) {
            int c = (int) quality.charAt(i);
            if (c < 33 || c > 126) {
                throw new IOException("quality scores must contain ASCII codes 33 to 126, found " + c
                        + " at line " + lineNumber);
            }
        }
    }
}
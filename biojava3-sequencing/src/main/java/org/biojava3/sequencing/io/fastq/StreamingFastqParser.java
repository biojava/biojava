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

import com.google.common.io.InputSupplier;

/**
 * Event based parser for FASTQ formatted sequences.
 *
 * @since 3.0.3
 */
final class StreamingFastqParser
{

    /**
     * Stream the specified input supplier.
     *
     * @param supplier input supplier, must not be null
     * @param variant FASTQ variant, must not be null
     * @param listener event based reader callback, must not be null
     * @throws IOException if an I/O error occurs
     */
    static <R extends Readable & Closeable> void stream(final InputSupplier<R> supplier,
                                                        final FastqVariant variant,
                                                        final StreamListener listener)
        throws IOException
    {
        if (supplier == null)
        {
            throw new IllegalArgumentException("supplier must not be null");
        }
        if (variant == null)
        {
            throw new IllegalArgumentException("variant must not be null");
        }
        if (listener == null)
        {
            throw new IllegalArgumentException("listener must not be null");
        }

        final FastqBuilder builder = new FastqBuilder().withVariant(variant);
        FastqParser.parse(supplier, new ParseListener()
            {
                @Override
                public void description(final String description) throws IOException
                {
                    builder.withDescription(description);
                }

                @Override
                public void sequence(final String sequence) throws IOException
                {
                    builder.withSequence(sequence);
                }

                @Override
                public void appendSequence(final String sequence) throws IOException
                {
                    builder.appendSequence(sequence);
                }

                @Override
                public void repeatDescription(final String repeatDescription) throws IOException
                {
                    String description = builder.getDescription();
                    if ((description != null) && (description.length() > 0) && (repeatDescription.length() > 0))
                    {
                        if (!description.equals(repeatDescription))
                        {
                            throw new IOException("repeat description must match description");
                        }
                    }
                }

                /**
                 * Validate the specified quality line.
                 *
                 * @param quality quality line to validate
                 * @throws IOException if an I/O error occurs
                 */
                private void validateQuality(final String quality) throws IOException
                {
                    for (int i = 0, size = quality.length(); i < size; i++)
                    {
                        char c = quality.charAt(i);
                        int qualityScore = variant.qualityScore(c);
                        if (qualityScore < variant.minimumQualityScore() || qualityScore > variant.maximumQualityScore())
                        {
                            throw new IOException("quality score must be between " + variant.minimumQualityScore() +
                                                  " and " + variant.maximumQualityScore() + ", was " + qualityScore +
                                                  " for ASCII char '" + c + "'");
                        }
                    }
                }

                @Override
                public void quality(final String quality) throws IOException
                {
                    validateQuality(quality);
                    builder.withQuality(quality);
                }

                @Override
                public void appendQuality(final String quality) throws IOException
                {
                    validateQuality(quality);
                    builder.appendQuality(quality);
                }

                @Override
                public void complete() throws IOException
                {
                    try
                    {
                        listener.fastq(builder.build());
                    }
                    catch (IllegalStateException e)
                    {
                        throw new IOException("caught an IllegalStateException", e);
                    }
                }
            });
    }
}

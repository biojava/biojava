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

/**
 * Fluent builder API for creating FASTQ formatted sequences.
 *
 * @since 3.0.3
 */
public final class FastqBuilder
{
    /** Description for this FASTQ formatted sequence builder. */
    private String description;

    /** Sequence for this FASTQ formatted sequence builder. */
    private StringBuilder sequence;

    /** Quality scores for this FASTQ formatted sequence builder. */
    private StringBuilder quality;

    /** FASTQ sequence format variant for this FASTQ formatted sequence builder. */
    private FastqVariant variant = DEFAULT_VARIANT;

    /** Default FASTQ sequence format variant, <code>FastqVariant.FASTQ_SANGER</code>. */
    public static final FastqVariant DEFAULT_VARIANT = FastqVariant.FASTQ_SANGER;


    /**
     * Create a new FASTQ formatted sequence builder.
     */
    public FastqBuilder()
    {
        // empty
    }


    /**
     * Return the description for this FASTQ formatted sequence builder.
     *
     * @return the description for this FASTQ formatted sequence builder
     */
    public String getDescription()
    {
        return description;
    }

    /**
     * Return this FASTQ formatted sequence builder configured with the specified description.
     *
     * @param description description for this FASTQ formatted sequence builder, must not be null
     * @return this FASTQ formatted sequence builder configured with the specified description
     */
    public FastqBuilder withDescription(final String description)
    {
        if (description == null)
        {
            throw new IllegalArgumentException("description must not be null");
        }
        this.description = description;
        return this;
    }

    /**
     * Return this FASTQ formatted sequence builder configured with the specified sequence.
     *
     * @param sequence sequence for this FASTQ formatted sequence builder, must not be null
     * @return this FASTQ formatted sequence builder configured with the specified sequence
     */
    public FastqBuilder withSequence(final String sequence)
    {
        if (sequence == null)
        {
            throw new IllegalArgumentException("sequence must not be null");
        }
        if (this.sequence == null)
        {
            this.sequence = new StringBuilder(sequence.length());
        }
        this.sequence.replace(0, this.sequence.length(), sequence);
        return this;
    }

    /**
     * Return this FASTQ formatted sequence builder configured with the specified sequence
     * appended to its current sequence.
     *
     * @param sequence sequence to append to the sequence for this FASTQ formatted sequence builder, must not be null
     * @return this FASTQ formatted sequence builder configured with the specified sequence
     *    appended to its current sequence
     */
    public FastqBuilder appendSequence(final String sequence)
    {
        if (sequence == null)
        {
            throw new IllegalArgumentException("sequence must not be null");
        }
        if (this.sequence == null)
        {
            this.sequence = new StringBuilder(sequence.length());
        }
        this.sequence.append(sequence);
        return this;
    }

    /**
     * Return this FASTQ formatted sequence builder configured with the specified quality scores.
     *
     * @param quality quality scores for this FASTQ formatted sequence builder, must not be null
     * @return this FASTQ formatted sequence builder configured with the specified quality scores
     */
    public FastqBuilder withQuality(final String quality)
    {
        if (quality == null)
        {
            throw new IllegalArgumentException("quality must not be null");
        }
        if (this.quality == null)
        {
            this.quality = new StringBuilder(quality.length());
        }
        this.quality.replace(0, this.quality.length(), quality);
        return this;
    }

    /**
     * Return this FASTQ formatted sequence builder configured with the specified quality scores
     * appended to its current quality scores.
     *
     * @param quality quality scores to append to the quality scores for this FASTQ formatted sequence
     *    builder, must not be null
     * @return this FASTQ formatted sequence builder configured with the specified quality scores
     *    appended to its current quality scores
     */
    public FastqBuilder appendQuality(final String quality)
    {
        if (quality == null)
        {
            throw new IllegalArgumentException("quality must not be null");
        }
        if (this.quality == null)
        {
            this.quality = new StringBuilder(quality.length());
        }
        this.quality.append(quality);
        return this;
    }

    /**
     * Return true if the sequence and quality scores for this FASTQ formatted sequence builder are equal in length.
     *
     * @return true if the sequence and quality scores for this FASTQ formatted sequence builder are equal in length
     */
    public boolean sequenceAndQualityLengthsMatch()
    {
        if (sequence == null && quality == null)
        {
            return true;
        }
        if ((sequence != null && quality == null) || (sequence == null && quality != null))
        {
            return false;
        }
        return sequence.length() == quality.length();
    }

    /**
     * Return this FASTQ formatted sequence builder configured with the specified FASTQ sequence format variant.
     *
     * @param variant FASTQ sequence format variant for this FASTQ formatted sequence builder, must not be null
     * @return this FASTQ formatted sequence builder configured with the specified FASTQ sequence format variant
     */
    public FastqBuilder withVariant(final FastqVariant variant)
    {
        if (variant == null)
        {
            throw new IllegalArgumentException("variant must not be null");
        }
        this.variant = variant;
        return this;
    }

    /**
     * Build and return a new FASTQ formatted sequence configured from the properties of this builder.
     *
     * @return a new FASTQ formatted sequence configured from the properties of this builder
     * @throws IllegalStateException if the configuration of this builder results in an illegal state
     */
    public Fastq build()
    {
        if (description == null)
        {
            throw new IllegalStateException("description must not be null");
        }
        if (sequence == null)
        {
            throw new IllegalStateException("sequence must not be null");
        }
        if (quality == null)
        {
            throw new IllegalStateException("quality must not be null");
        }
        if (!sequenceAndQualityLengthsMatch())
        {
            throw new IllegalStateException("sequence and quality scores must be the same length");
        }
        Fastq fastq = new Fastq(description, sequence.toString(), quality.toString(), variant);
        return fastq;
    }
}
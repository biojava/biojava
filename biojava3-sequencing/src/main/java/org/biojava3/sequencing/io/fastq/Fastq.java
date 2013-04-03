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
 * FASTQ formatted sequence.
 *
 * @since 3.0.3
 */
public final class Fastq
{
    /** Description of this FASTQ formatted sequence. */
    private final String description;

    /** Sequence for this FASTQ formatted sequence. */
    private final String sequence;

    /** Quality scores for this FASTQ formatted sequence. */
    private final String quality;

    /** FASTQ sequence format variant for this FASTQ formatted sequence. */
    private final FastqVariant variant;


    /**
     * Create a new FASTQ formatted sequence from the specified description, sequence, quality scores,
     * and sequence format variant.
     *
     * @param description description of this FASTQ formatted sequence, must not be null
     * @param sequence sequence for this FASTQ formatted sequence, must not be null
     * @param quality quality scores for this FASTQ formatted sequence, must not be null
     * @param variant FASTQ sequence format variant for this FASTQ formatted sequence, must not be null
     */
    Fastq(final String description,
          final String sequence,
          final String quality,
          final FastqVariant variant)
    {
        if (description == null)
        {
            throw new IllegalArgumentException("description must not be null");
        }
        if (sequence == null)
        {
            throw new IllegalArgumentException("sequence must not be null");
        }
        if (quality == null)
        {
            throw new IllegalArgumentException("quality must not be null");
        }
        if (variant == null)
        {
            throw new IllegalArgumentException("variant must not be null");
        }
        this.description = description;
        this.sequence = sequence;
        this.quality = quality;
        this.variant = variant;
    }


    /**
     * Return the description of this FASTQ formatted sequence.
     * The description will not be null.
     *
     * @return the description of this FASTQ formatted sequence
     */
    public String getDescription()
    {
        return description;
    }

    /**
     * Return the sequence for this FASTQ formatted sequence.
     * The sequence will not be null.
     *
     * @return the sequence for this FASTQ formatted sequence
     */
    public String getSequence()
    {
        return sequence;
    }

    /**
     * Return the quality scores for this FASTQ formatted sequence.
     * The quality scores will not be null.
     *
     * @return the quality scores for this FASTQ formatted sequence
     */
    public String getQuality()
    {
        return quality;
    }

    /**
     * Return the FASTQ sequence format variant for this FASTQ formatted sequence.
     * The FASTQ sequence format variant will not be null.
     *
     * @return the FASTQ sequence format variant for this FASTQ formatted sequence
     */
    public FastqVariant getVariant()
    {
        return variant;
    }
}
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

import java.util.HashMap;
import java.util.Map;

/**
 * FASTQ sequence format variant.
 *
 * @since 3.0.3
 */
public enum FastqVariant
{
    /** Sanger FASTQ sequence format variant. */
    FASTQ_SANGER("Original or Sanger format")
    {
        @Override
        public int minimumQualityScore()
        {
            return 0;
        }

        @Override
        public int maximumQualityScore()
        {
            return 93;
        }

        @Override
        public int qualityScore(final char c)
        {
            return ((int) c) - 33;
        }

        @Override
        public char quality(final int qualityScore)
        {
            if (qualityScore < minimumQualityScore())
            {
                throw new IllegalArgumentException("qualityScore must be greater than or equal to minimumQualityScore()");
            }
            if (qualityScore > maximumQualityScore())
            {
                throw new IllegalArgumentException("qualityScore must be less than or equal to maximumQualityScore()");
            }
            return (char) (qualityScore + 33);
        }

        @Override
        public double errorProbability(final int qualityScore)
        {
            return Math.pow(10.0d, ((double) qualityScore) / -10.0d);
        }
    },

    /** Solexa FASTQ sequence format variant. */
    FASTQ_SOLEXA("Solexa and early Illumina format")
    {
        @Override
        public int minimumQualityScore()
        {
            return -5;
        }

        @Override
        public int maximumQualityScore()
        {
            return 62;
        }

        @Override
        public int qualityScore(final char c)
        {
            return ((int) c) - 64;
        }

        @Override
        public char quality(final int qualityScore)
        {
            if (qualityScore < minimumQualityScore())
            {
                throw new IllegalArgumentException("qualityScore must be greater than or equal to minimumQualityScore()");
            }
            if (qualityScore > maximumQualityScore())
            {
                throw new IllegalArgumentException("qualityScore must be less than or equal to maximumQualityScore()");
            }
            return (char) (qualityScore + 64);
        }

        @Override
        public double errorProbability(final int qualityScore)
        {
            double q = Math.pow(10.0d, ((double) qualityScore) / -10.0d);
            return q / (1 + q);
        }
    },

    /** Illumina FASTQ sequence format variant. */
    FASTQ_ILLUMINA("Illumina 1.3+ format")
    {
        @Override
        public int minimumQualityScore()
        {
            return 0;
        }

        @Override
        public int maximumQualityScore()
        {
            return 62;
        }

        @Override
        public int qualityScore(final char c)
        {
            return ((int) c) - 64;
        }

        @Override
        public char quality(final int qualityScore)
        {
            if (qualityScore < minimumQualityScore())
            {
                throw new IllegalArgumentException("qualityScore must be greater than or equal to minimumQualityScore()");
            }
            if (qualityScore > maximumQualityScore())
            {
                throw new IllegalArgumentException("qualityScore must be less than or equal to maximumQualityScore()");
            }
            return (char) (qualityScore + 64);
        }

        @Override
        public double errorProbability(final int qualityScore)
        {
            return Math.pow(10.0d, ((double) qualityScore) / -10.0d);
        }
    };


    /** Map of FASTQ sequence format variants keyed by name and lowercase-with-dashes name. */
    private static final Map<String, FastqVariant> FASTQ_VARIANTS = new HashMap<String, FastqVariant>(6);

    static
    {
        for (FastqVariant fastqVariant : values())
        {
            FASTQ_VARIANTS.put(fastqVariant.name(), fastqVariant);
            FASTQ_VARIANTS.put(fastqVariant.lowercaseName(), fastqVariant);
        }
    }

    /** Description of this FASTQ sequence format variant. */
    private final String description;


    /**
     * Create a new FASTQ sequence format variant with the specified description.
     *
     * @param description description of this FASTQ sequence format variant, must not be null
     */
    private FastqVariant(final String description)
    {
        if (description == null)
        {
            throw new IllegalArgumentException("description must not be null");
        }
        this.description = description;
    }


    /**
     * Return the description of this FASTQ sequence format variant.
     * The description will not be null.
     *
     * @return the description of this FASTQ sequence format variant
     */
    public String getDescription()
    {
        return description;
    }

    /**
     * Return true if this FASTQ sequence format variant is {@link #FASTQ_SANGER}.
     *
     * @return true if this FASTQ sequence format variant is {@link #FASTQ_SANGER}
     */
    public boolean isSanger()
    {
        return (this == FASTQ_SANGER);
    }

    /**
     * Return true if this FASTQ sequence format variant is {@link #FASTQ_SOLEXA}.
     *
     * @return true if this FASTQ sequence format variant is {@link #FASTQ_SOLEXA}
     */
    public boolean isSolexa()
    {
        return (this == FASTQ_SOLEXA);
    }

    /**
     * Return true if this FASTQ sequence format variant is {@link #FASTQ_ILLUMINA}.
     *
     * @return true if this FASTQ sequence format variant is {@link #FASTQ_ILLUMINA}
     */
    public boolean isIllumina()
    {
        return (this == FASTQ_ILLUMINA);
    }

    /**
     * Return the minimum quality score for this FASTQ sequence format variant.
     *
     * @return the minimum quality score for this FASTQ sequence format variant.
     */
    public abstract int minimumQualityScore();

    /**
     * Return the maximum quality score for this FASTQ sequence format variant.
     *
     * @return the maximum quality score for this FASTQ sequence format variant.
     */
    public abstract int maximumQualityScore();

    /**
     * Convert the specified quality in ASCII format to a quality score.
     *
     * @param c quality in ASCII format
     * @return the specified quality in ASCII format converted to a quality score
     */
    public abstract int qualityScore(char c);

    /**
     * Convert the specified quality score to a quality in ASCII format.
     *
     * @since 3.0.6
     * @param qualityScore quality score, must be <code>&gt;= minimumQualityScore()</code>
     *    and <code>&lt;= maximumQualityScore()</code>
     * @return the quality in ASCII format converted from the specified quality score
     */
    public abstract char quality(int qualityScore);

    /**
     * Convert the specified quality in ASCII format to an error probability.
     *
     * @param c quality in ASCII format
     * @return the specified quality in ASCII format converted to an error probability
     */
    public double errorProbability(char c)
    {
        return errorProbability(qualityScore(c));
    }

    /**
     * Calculate the error probability given the specified quality score.
     *
     * @param qualityScore quality score
     * @return the error probability given the specified quality score
     */
    public abstract double errorProbability(int qualityScore);

    /**
     * Return the name of this FASTQ sequence format variant in <code>lowercase-with-dashes</code> style.
     *
     * @return the name of this FASTQ sequence format variant in <code>lowercase-with-dashes</code> style
     */
    public String lowercaseName()
    {
        return name().toLowerCase().replace('_', '-');
    }


    /**
     * Return the FASTQ sequence format variant with the specified name, if any.  The name may
     * be specified in either <code>UPPERCASE_WITH_UNDERSCORES</code>
     * or <code>lowercase-with-dashes</code> style.
     *
     * @param name name
     * @return the FASTQ sequence format variant with the specified name, or <code>null</code>
     *    if no such FASTQ sequence format variant exists
     */
    public static FastqVariant parseFastqVariant(final String name)
    {
        return FASTQ_VARIANTS.get(name);
    }
}
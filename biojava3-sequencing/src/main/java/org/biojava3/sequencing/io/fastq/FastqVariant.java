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
 * @since 1.7.1
 */
public enum FastqVariant
{
    /** Sanger FASTQ sequence format variant. */
    FASTQ_SANGER("Original or Sanger format"),

    /** Solexa FASTQ sequence format variant. */
    FASTQ_SOLEXA("Solexa and early Illumina format"),

    /** Illumina FASTQ sequence format variant. */
    FASTQ_ILLUMINA("Illumina 1.3+ format");


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
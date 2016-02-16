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
package org.biojava.nbio.sequencing.io.fastq;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.features.QualityFeature;
import org.biojava.nbio.core.sequence.features.QuantityFeature;
import org.biojava.nbio.core.sequence.template.AbstractSequence;

import java.util.List;

/**
 * Utility methods for FASTQ formatted sequences.
 *
 * @since 3.0.3
 */
public final class FastqTools
{

    /**
     * Private no-arg constructor.
     */
    private FastqTools()
    {
        // empty
    }


    /**
     * Create and return a new {@link DNASequence} from the specified FASTQ formatted sequence.
     *
     * @param fastq FASTQ formatted sequence, must not be null
     * @return a new {@link DNASequence} from the specified FASTQ formatted sequence
     * @throws CompoundNotFoundException if DNA sequence in fastq contains unrecognised compounds
     */
    public static DNASequence createDNASequence(final Fastq fastq) throws CompoundNotFoundException
    {
        if (fastq == null)
        {
            throw new IllegalArgumentException("fastq must not be null");
        }
        DNASequence sequence = new DNASequence(fastq.getSequence());
        sequence.setOriginalHeader(fastq.getDescription());
        return sequence;
    }

    /**
     * Create and return a new {@link DNASequence} with quality scores from the specified
     * FASTQ formatted sequence.  The quality scores are stored in a {@link QualityFeature}
     * with a type <code>"qualityScores"</code> the same length as the sequence.
     *
     * @param fastq FASTQ formatted sequence, must not be null
     * @return a new {@link DNASequence} with quality scores from the specified FASTQ formatted sequence
     * @throws CompoundNotFoundException if DNA sequence in fastq contains unrecognised compounds 
     */
    public static DNASequence createDNASequenceWithQualityScores(final Fastq fastq) throws CompoundNotFoundException
    {
        DNASequence sequence = createDNASequence(fastq);
        sequence.addFeature(1, sequence.getLength(), createQualityScores(fastq));
        return sequence;
    }

    /**
     * Create and return a new {@link DNASequence} with error probabilities from the specified
     * FASTQ formatted sequence.  The error probabilities are stored in a {@link QuantityFeature}
     * with a type <code>"errorProbabilities"</code> the same length as the sequence.
     *
     * @param fastq FASTQ formatted sequence, must not be null
     * @return a new {@link DNASequence} with error probabilities from the specified FASTQ formatted sequence
     * @throws CompoundNotFoundException if DNA sequence in fastq contains unrecognised compounds 
     */
    public static DNASequence createDNASequenceWithErrorProbabilities(final Fastq fastq) throws CompoundNotFoundException
    {
        DNASequence sequence = createDNASequence(fastq);
        sequence.addFeature(1, sequence.getLength(), createErrorProbabilities(fastq));
        return sequence;
    }

    /**
     * Create and return a new {@link DNASequence} with quality scores and error probabilities from the
     * specified FASTQ formatted sequence.  The quality scores are stored in a {@link QualityFeature}
     * with a type <code>"qualityScores"</code> the same length as the sequence and the error
     * probabilities are stored in a {@link QuantityFeature} with a type <code>"errorProbabilities"</code>
     * the same length as the sequence.
     *
     * @param fastq FASTQ formatted sequence, must not be null
     * @return a new {@link DNASequence} with quality scores and error probabilities from the specified
     *    FASTQ formatted sequence
     * @throws CompoundNotFoundException if DNA sequence in fastq contains unrecognised compounds     
     */
    public static DNASequence createDNASequenceWithQualityScoresAndErrorProbabilities(final Fastq fastq) throws CompoundNotFoundException
    {
        DNASequence sequence = createDNASequence(fastq);
        sequence.addFeature(1, sequence.getLength(), createQualityScores(fastq));
        sequence.addFeature(1, sequence.getLength(), createErrorProbabilities(fastq));
        return sequence;
    }

    /**
     * Create and return a new {@link QualityFeature} from the quality scores of the specified
     * FASTQ formatted sequence.  The quality scores feature has a type <code>"qualityScores"</code>
     * and will be the same length as the sequence.
     *
     * @param fastq FASTQ formatted sequence, must not be null
     * @return a new {@link QualityFeature} from the quality scores of the specified FASTQ
     *    formatted sequence
     */
    public static QualityFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> createQualityScores(final Fastq fastq)
    {
        if (fastq == null)
        {
            throw new IllegalArgumentException("fastq must not be null");
        }
        QualityFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> qualityScores = new QualityFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound>("qualityScores", "sequencing");
        qualityScores.setQualities(toList(qualityScores(fastq)));
        return qualityScores;
    }

    /**
     * Create and return a new {@link QuantityFeature} from the error probabilities of the specified
     * FASTQ formatted sequence.  The error probabilities feature has a type <code>"errorProbabilities"</code>
     * and will be the same length as the sequence.
     *
     * @param fastq FASTQ formatted sequence, must not be null
     * @return a new {@link QualityFeature} from the error probabilities of the specified FASTQ
     *    formatted sequence
     */
    public static QuantityFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> createErrorProbabilities(final Fastq fastq)
    {
        if (fastq == null)
        {
            throw new IllegalArgumentException("fastq must not be null");
        }
        QuantityFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> errorProbabilities = new QuantityFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound>("errorProbabilities", "sequencing");
        errorProbabilities.setQuantities(toList(errorProbabilities(fastq)));
        return errorProbabilities;
    }

    /**
     * Return the quality scores from the specified FASTQ formatted sequence.
     *
     * @param fastq FASTQ formatted sequence, must not be null
     * @return the quality scores from the specified FASTQ formatted sequence
     */
    public static Iterable<Number> qualityScores(final Fastq fastq)
    {
        if (fastq == null)
        {
            throw new IllegalArgumentException("fastq must not be null");
        }
        int size = fastq.getQuality().length();
        List<Number> qualityScores = Lists.newArrayListWithExpectedSize(size);
        FastqVariant variant = fastq.getVariant();
        for (int i = 0; i < size; i++)
        {
            char c = fastq.getQuality().charAt(i);
            qualityScores.add(variant.qualityScore(c));
        }
        return ImmutableList.copyOf(qualityScores);
    }

    /**
     * Copy the quality scores from the specified FASTQ formatted sequence into the specified int array.
     *
     * @param fastq FASTQ formatted sequence, must not be null
     * @param qualityScores int array of quality scores, must not be null and must be the same
     *    length as the FASTQ formatted sequence quality
     * @return the specified int array of quality scores
     */
    public static int[] qualityScores(final Fastq fastq, final int[] qualityScores)
    {
        if (fastq == null)
        {
            throw new IllegalArgumentException("fastq must not be null");
        }
        if (qualityScores == null)
        {
            throw new IllegalArgumentException("qualityScores must not be null");
        }
        int size = fastq.getQuality().length();
        if (qualityScores.length != size)
        {
            throw new IllegalArgumentException("qualityScores must be the same length as the FASTQ formatted sequence quality");
        }
        FastqVariant variant = fastq.getVariant();
        for (int i = 0; i < size; i++)
        {
            char c = fastq.getQuality().charAt(i);
            qualityScores[i] = variant.qualityScore(c);
        }
        return qualityScores;
    }

    /**
     * Return the error probabilities from the specified FASTQ formatted sequence.
     *
     * @param fastq FASTQ formatted sequence, must not be null
     * @return the error probabilities from the specified FASTQ formatted sequence
     */
    public static Iterable<Number> errorProbabilities(final Fastq fastq)
    {
        if (fastq == null)
        {
            throw new IllegalArgumentException("fastq must not be null");
        }
        int size = fastq.getQuality().length();
        List<Number> errorProbabilities = Lists.newArrayListWithExpectedSize(size);
        FastqVariant variant = fastq.getVariant();
        for (int i = 0; i < size; i++)
        {
            char c = fastq.getQuality().charAt(i);
            errorProbabilities.add(variant.errorProbability(c));
        }
        return ImmutableList.copyOf(errorProbabilities);
    }

    /**
     * Copy the error probabilities from the specified FASTQ formatted sequence into the specified double array.
     *
     * @param fastq FASTQ formatted sequence, must not be null
     * @param errorProbabilities double array of error probabilities, must not be null and must be the same
     *    length as the FASTQ formatted sequence quality
     * @return the specified double array of error probabilities
     */
    public static double[] errorProbabilities(final Fastq fastq, final double[] errorProbabilities)
    {
        if (fastq == null)
        {
            throw new IllegalArgumentException("fastq must not be null");
        }
        if (errorProbabilities == null)
        {
            throw new IllegalArgumentException("errorProbabilities must not be null");
        }
        int size = fastq.getQuality().length();
        if (errorProbabilities.length != size)
        {
            throw new IllegalArgumentException("errorProbabilities must be the same length as the FASTQ formatted sequence quality");
        }
        FastqVariant variant = fastq.getVariant();
        for (int i = 0; i < size; i++)
        {
            char c = fastq.getQuality().charAt(i);
            errorProbabilities[i] = variant.errorProbability(c);
        }
        return errorProbabilities;
    }

    /**
     * Convert the specified FASTQ formatted sequence to the
     * specified FASTQ sequence format variant.
     *
     * @since 4.2
     * @param fastq FASTQ formatted sequence, must not be null
     * @param variant FASTQ sequence format variant, must not be null
     * @return the specified FASTQ formatted sequence converted to the
     *    specified FASTQ sequence format variant
     */
    public static Fastq convert(final Fastq fastq, final FastqVariant variant)
    {
        if (fastq == null)
        {
            throw new IllegalArgumentException("fastq must not be null");
        }
        if (variant == null)
        {
            throw new IllegalArgumentException("variant must not be null");
        }
        if (fastq.getVariant().equals(variant))
        {
            return fastq;
        }
        return new Fastq(fastq.getDescription(), fastq.getSequence(), convertQualities(fastq, variant), variant);
    }

    /**
     * Convert the qualities in the specified FASTQ formatted sequence to the
     * specified FASTQ sequence format variant.
     *
     * @since 4.2
     * @param fastq FASTQ formatted sequence, must not be null
     * @param variant FASTQ sequence format variant, must not be null
     * @return the qualities in the specified FASTQ formatted sequence converted to the
     *    specified FASTQ sequence format variant
     */
    static String convertQualities(final Fastq fastq, final FastqVariant variant)
    {
        if (fastq == null)
        {
            throw new IllegalArgumentException("fastq must not be null");
        }
        if (variant == null)
        {
            throw new IllegalArgumentException("variant must not be null");
        }
        if (fastq.getVariant().equals(variant))
        {
            return fastq.getQuality();
        }
        int size = fastq.getQuality().length();
        double[] errorProbabilities = errorProbabilities(fastq, new double[size]);
        StringBuilder sb = new StringBuilder(size);
        for (int i = 0; i < size; i++)
        {
            sb.append(variant.quality(variant.qualityScore(errorProbabilities[i])));
        }
        return sb.toString();
    }

    /**
     * Return the specified iterable as a list.
     *
     * @paam <T> element type
     * @param iterable iterable
     * @return the specified iterable as a list
     */
    @SuppressWarnings("unchecked")
    static <T> List<T> toList(final Iterable<? extends T> iterable)
    {
        if (iterable instanceof List)
        {
            return (List<T>) iterable;
        }
        return ImmutableList.copyOf(iterable);
    }
}

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

import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.ImmutableList;

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
     * Return the quality scores from the specified FASTQ formatted sequence.
     *
     * @param fastq FASTQ formatted sequence, must not be null
     * @return the quality scores from the specified FASTQ formatted sequence
     */
    public static Iterable<Integer> qualityScores(final Fastq fastq)
    {
        if (fastq == null)
        {
            throw new IllegalArgumentException("fastq must not be null");
        }
        int size = fastq.getQuality().length();
        List<Integer> qualityScores = Lists.newArrayListWithExpectedSize(size);
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
    public static Iterable<Double> errorProbabilities(final Fastq fastq)
    {
        if (fastq == null)
        {
            throw new IllegalArgumentException("fastq must not be null");
        }
        int size = fastq.getQuality().length();
        List<Double> errorProbabilities = Lists.newArrayListWithExpectedSize(size);
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
}
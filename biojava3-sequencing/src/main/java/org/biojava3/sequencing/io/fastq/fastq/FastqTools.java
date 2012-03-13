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
package org.biojava.bio.program.fastq;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.ImmutableList;

import org.biojava.bio.Annotation;

import org.biojava.bio.dist.Distribution;

import org.biojava.bio.program.phred.PhredSequence;
import org.biojava.bio.program.phred.PhredTools;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Sequence;

import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.IntegerAlphabet;
import org.biojava.bio.symbol.IntegerAlphabet.SubIntegerAlphabet;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.bio.symbol.SimpleSymbolList;

/**
 * Utility methods for FASTQ formatted sequences.
 *
 * @since 1.8.2
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
     * Create and return a new DNA {@link SymbolList} from the specified FASTQ formatted sequence.
     *
     * @param fastq FASTQ formatted sequence, must not be null
     * @return a new DNA {@link SymbolList} from the specified FASTQ formatted sequence
     * @throws IllegalSymbolException if an illegal symbol is found
     */
    public static SymbolList createDNA(final Fastq fastq) throws IllegalSymbolException
    {
        if (fastq == null)
        {
            throw new IllegalArgumentException("fastq must not be null");
        }
        return DNATools.createDNA(fastq.getSequence());
    }

    /**
     * Create and return a new {@link SymbolList} of quality scores from the specified FASTQ formatted sequence.
     *
     * @param fastq FASTQ formatted sequence, must not be null
     * @return a new {@link SymbolList} of quality scores from the specified FASTQ formatted sequence
     * @throws IllegalSymbolException if an illegal symbol is found
     */
    public static SymbolList createQuality(final Fastq fastq) throws IllegalSymbolException
    {
        if (fastq == null)
        {
            throw new IllegalArgumentException("fastq must not be null");
        }
        FastqVariant variant = fastq.getVariant();
        SubIntegerAlphabet alphabet = IntegerAlphabet.getSubAlphabet(variant.minimumQualityScore(), variant.maximumQualityScore());
        SimpleSymbolList qualitySymbols = new SimpleSymbolList(alphabet);
        for (int i = 0, size = fastq.getQuality().length(); i < size; i++)
        {
            char c = fastq.getQuality().charAt(i);
            qualitySymbols.addSymbol(alphabet.getSymbol(variant.qualityScore(c)));
        }
        return qualitySymbols;
    }

    /**
     * Create and return a new DNA {@link Sequence} from the specified FASTQ formatted sequence.
     *
     * @param fastq FASTQ formatted sequence, must not be null
     * @return a new {@link Sequence} from the specified FASTQ formatted sequence
     * @throws IllegalSymbolException if an illegal symbol is found
     */
    public static Sequence createDNASequence(final Fastq fastq) throws IllegalSymbolException
    {
        if (fastq == null)
        {
            throw new IllegalArgumentException("fastq must not be null");
        }
        return DNATools.createDNASequence(fastq.getSequence(), fastq.getDescription());
    }

    /**
     * Create and return a new {@link PhredSequence} from the specified FASTQ formatted sequence.
     * Only Sanger variant FASTQ formatted sequences are supported.
     *
     * @param fastq FASTQ formatted sequence, must not be null and must be Sanger variant
     * @return a new {@link PhredSequence} from the specified FASTQ formatted sequence
     * @throws IllegalAlphabetException if an illegal alphabet is used
     * @throws IllegalSymbolException if an illegal symbol is found
     */
    public static PhredSequence createPhredSequence(final Fastq fastq) throws IllegalAlphabetException, IllegalSymbolException
    {
        if (fastq == null)
        {
            throw new IllegalArgumentException("fastq must not be null");
        }
        if (!fastq.getVariant().isSanger())
        {
            throw new IllegalArgumentException("fastq must be sanger variant, was " + fastq.getVariant());
        }
        SymbolList dnaSymbols = createDNA(fastq);

        // 0-99 subinteger alphabet required by PhredSequence, thus only Sanger variant is supported
        SubIntegerAlphabet alphabet = IntegerAlphabet.getSubAlphabet(0, 99);
        SimpleSymbolList qualitySymbols = new SimpleSymbolList(alphabet);
        for (int i = 0, size = fastq.getQuality().length(); i < size; i++)
        {
            char c = fastq.getQuality().charAt(i);
            qualitySymbols.addSymbol(alphabet.getSymbol(FastqVariant.FASTQ_SANGER.qualityScore(c)));
        }

        SymbolList phredSymbols = PhredTools.createPhred(dnaSymbols, qualitySymbols);
        return new PhredSequence(phredSymbols, fastq.getDescription(), null, Annotation.EMPTY_ANNOTATION);
    }

    /**
     * Create and return a new array of symbol {@link Distribution}s from the specified FASTQ formatted sequence.
     * Only Sanger variant FASTQ formatted sequences are supported.
     *
     * @param fastq FASTQ formatted sequence, must not be null and must be Sanger variant
     * @return a new array of symbol {@link Distribution}s from the specified FASTQ formatted sequence
     * @throws IllegalAlphabetException if an illegal alphabet is used
     * @throws IllegalSymbolException if an illegal symbol is found
     */
    public static Distribution[] createSymbolDistribution(final Fastq fastq) throws IllegalAlphabetException, IllegalSymbolException
    {
        PhredSequence phredSequence = createPhredSequence(fastq);
        return PhredTools.phredToDistArray(phredSequence);
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
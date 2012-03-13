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

import junit.framework.TestCase;

import org.biojava.bio.dist.Distribution;

import org.biojava.bio.program.phred.PhredSequence;
import org.biojava.bio.seq.Sequence;

import org.biojava.bio.symbol.SymbolList;

/**
 * Unit test for FastqTools.
 */
public final class FastqToolsTest extends TestCase
{
    private final FastqBuilder builder = new FastqBuilder().withDescription("foo").withSequence("ACTG").withQuality("ZZZZ");

    public void testCreateDNA() throws Exception
    {
        SymbolList dna = FastqTools.createDNA(builder.build());
        assertNotNull(dna);
        assertEquals(4, dna.length());
    }

    public void testCreateDNANullFastq() throws Exception
    {
        try
        {
            FastqTools.createDNA(null);
            fail("createDNA(null) expected IllegalArgumentException");
        }
        catch (IllegalArgumentException e)
        {
            // expected
        }
    }

    public void testCreateQuality() throws Exception
    {
        SymbolList quality = FastqTools.createQuality(builder.build());
        assertNotNull(quality);
        assertEquals(4, quality.length());
    }

    public void testCreateQualityNullFastq() throws Exception
    {
        try
        {
            FastqTools.createQuality(null);
            fail("createQuality(null) expected IllegalArgumentException");
        }
        catch (IllegalArgumentException e)
        {
            // expected
        }
    }

    public void testCreateDNASequence() throws Exception
    {
        Sequence sequence = FastqTools.createDNASequence(builder.build());
        assertNotNull(sequence);
        assertEquals(4, sequence.length());
    }

    public void testCreateDNASequenceNullFastq() throws Exception
    {
        try
        {
            FastqTools.createDNASequence(null);
            fail("createDNASequence(null) expected IllegalArgumentException");
        }
        catch (IllegalArgumentException e)
        {
            // expected
        }
    }

    public void testPhredSequence() throws Exception
    {
        PhredSequence sequence = FastqTools.createPhredSequence(builder.build());
        assertNotNull(sequence);
        assertEquals(4, sequence.length());
    }

    public void testCreatePhredSequenceNullFastq() throws Exception
    {
        try
        {
            FastqTools.createPhredSequence(null);
            fail("createPhredSequence(null) expected IllegalArgumentException");
        }
        catch (IllegalArgumentException e)
        {
            // expected
        }
    }

    public void testPhredSequenceNotSangerVariant() throws Exception
    {
        try
        {
            FastqTools.createPhredSequence(builder.withVariant(FastqVariant.FASTQ_ILLUMINA).build());
            fail("createPhredSequence(not sanger variant) expected IllegalArgumentException");
        }
        catch (IllegalArgumentException e)
        {
            // expected
        }
    }

    public void testCreateSymbolDistribution() throws Exception
    {
        Distribution[] symbolDistribution = FastqTools.createSymbolDistribution(builder.build());
        assertNotNull(symbolDistribution);
        assertEquals(4, symbolDistribution.length);
    }

    public void testCreateSymbolDistributionNullFastq() throws Exception
    {
        try
        {
            FastqTools.createSymbolDistribution(null);
            fail("createSymbolDistribution(null) expected IllegalArgumentException");
        }
        catch (IllegalArgumentException e)
        {
            // expected
        }
    }

    public void testCreateSymbolDistributionNotSangerVariant() throws Exception
    {
        try
        {
            FastqTools.createSymbolDistribution(builder.withVariant(FastqVariant.FASTQ_ILLUMINA).build());
            fail("createSymbolDistribution(not sanger variant) expected IllegalArgumentException");
        }
        catch (IllegalArgumentException e)
        {
            // expected
        }
    }

    public void testQualityScores()
    {
        Iterable<Integer> qualityScores = FastqTools.qualityScores(builder.build());
        assertNotNull(qualityScores);
        int count = 0;
        for (Integer qualityScore : qualityScores)
        {
            assertNotNull(qualityScore);
            count++;
        }
        assertEquals(4, count);
    }

    public void testQualityScoresNullFastq()
    {
        try
        {
            FastqTools.qualityScores(null);
            fail("qualityScores(null) expected IllegalArgumentException");
        }
        catch (IllegalArgumentException e)
        {
            // expected
        }
    }

    public void testQualityScoresIntArray()
    {
        int[] qualityScores = new int[4];
        FastqTools.qualityScores(builder.build(), qualityScores);
        for (int i = 0; i < 4; i++)
        {
            assertTrue(qualityScores[i] != 0);
        }
    }

    public void testQualityScoresIntArrayNullFastq()
    {
        try
        {
            FastqTools.qualityScores(null, new int[0]);
            fail("qualityScores(null, int[]) expected IllegalArgumentException");
        }
        catch (IllegalArgumentException e)
        {
            // expected
        }
    }

    public void testQualityScoresNullIntArray()
    {
        try
        {
            FastqTools.qualityScores(builder.build(), null);
            fail("qualityScores(fastq, null) expected IllegalArgumentException");
        }
        catch (IllegalArgumentException e)
        {
            // expected
        }
    }

    public void testQualityScoresQualityScoresTooSmall()
    {
        try
        {
            FastqTools.qualityScores(builder.build(), new int[3]);
            fail("expected IllegalArgumentException");
        }
        catch (IllegalArgumentException e)
        {
            // expected
        }
    }

    public void testQualityScoresQualityScoresTooLarge()
    {
        try
        {
            FastqTools.qualityScores(builder.build(), new int[5]);
            fail("expected IllegalArgumentException");
        }
        catch (IllegalArgumentException e)
        {
            // expected
        }
    }

    public void testErrorProbabilities()
    {
        Iterable<Double> errorProbabilities = FastqTools.errorProbabilities(builder.build());
        assertNotNull(errorProbabilities);
        int count = 0;
        for (Double errorProbability : errorProbabilities)
        {
            assertNotNull(errorProbability);
            count++;
        }
        assertEquals(4, count);
    }

    public void testErrorProbabilitiesNullFastq()
    {
        try
        {
            FastqTools.errorProbabilities(null);
            fail("errorProbabilities(null) expected IllegalArgumentException");
        }
        catch (IllegalArgumentException e)
        {
            // expected
        }
    }

    public void testErrorProbabilitiesDoubleArray()
    {
        double[] errorProbabilities = new double[4];
        FastqTools.errorProbabilities(builder.build(), errorProbabilities);
        for (int i = 0; i < 0; i++)
        {
            assertTrue(errorProbabilities[i] > 0.0d);
        }
    }

    public void testErrorProbabilitiesDoubleArrayNullFastq()
    {
        try
        {
            FastqTools.errorProbabilities(null, new double[0]);
            fail("errorProbabilities(null, double[]) expected IllegalArgumentException");
        }
        catch (IllegalArgumentException e)
        {
            // expected
        }
    }

    public void testErrorProbabilitiesNullErrorProbabilities()
    {
        try
        {
            FastqTools.errorProbabilities(builder.build(), null);
            fail("errorProbabilities(fastq, null) expected IllegalArgumentException");
        }
        catch (IllegalArgumentException e)
        {
            // expected
        }
    }

    public void testErrorProbabilitiesErrorProbabilitiesTooSmall()
    {
        try
        {
            FastqTools.errorProbabilities(builder.build(), new double[3]);
            fail("expected IllegalArgumentException");
        }
        catch (IllegalArgumentException e)
        {
            // expected
        }
    }

    public void testErrorProbabilitiesErrorProbabilitiesTooLarge()
    {
        try
        {
            FastqTools.errorProbabilities(builder.build(), new double[5]);
            fail("expected IllegalArgumentException");
        }
        catch (IllegalArgumentException e)
        {
            // expected
        }
    }
}

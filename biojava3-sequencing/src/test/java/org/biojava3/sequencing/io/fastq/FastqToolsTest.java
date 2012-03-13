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

import junit.framework.TestCase;

/**
 * Unit test for FastqTools.
 */
public final class FastqToolsTest extends TestCase
{
    private final FastqBuilder builder = new FastqBuilder().withDescription("foo").withSequence("ACTG").withQuality("ZZZZ");

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

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

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;

import junit.framework.TestCase;

import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.features.QualityFeature;
import org.biojava3.core.sequence.features.QuantityFeature;

/**
 * Unit test for FastqTools.
 */
public final class FastqToolsTest extends TestCase
{
    private final FastqBuilder builder = new FastqBuilder().withDescription("foo").withSequence("ACTG").withQuality("ZZZZ");

    public void testCreateDNASequence()
    {
        DNASequence sequence = FastqTools.createDNASequence(builder.build());
        assertNotNull(sequence);
    }

    public void testCreateDNASequenceNullFastq()
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

    public void testCreateDNASequenceWithQualityScores()
    {
        DNASequence sequence = FastqTools.createDNASequenceWithQualityScores(builder.build());
        assertNotNull(sequence);

        List features = sequence.getFeaturesByType("qualityScores");
        assertNotNull(features);
        assertEquals(1, features.size());
        QualityFeature qualityScores = (QualityFeature) features.get(0);
        assertEquals(sequence.getLength(), qualityScores.getQualities().size());
        assertEquals(sequence.getLength(), qualityScores.getLocations().getLength());
    }

    public void testCreateDNASequenceWithQualityScoresNullFastq()
    {
        try
        {
            FastqTools.createDNASequenceWithQualityScores(null);
            fail("createDNASequenceWithQualityScores(null) expected IllegalArgumentException");
        }
        catch (IllegalArgumentException e)
        {
            // expected
        }
    }

    public void testCreateDNASequenceWithErrorProbabilies()
    {
        DNASequence sequence = FastqTools.createDNASequenceWithErrorProbabilities(builder.build());
        assertNotNull(sequence);

        List features = sequence.getFeaturesByType("errorProbabilities");
        assertNotNull(features);
        assertEquals(1, features.size());
        QuantityFeature errorProbabilities = (QuantityFeature) features.get(0);
        assertEquals(sequence.getLength(), errorProbabilities.getQuantities().size());
        assertEquals(sequence.getLength(), errorProbabilities.getLocations().getLength());
    }

    public void testCreateDNASequenceWithErrorProbabilitiesNullFastq()
    {
        try
        {
            FastqTools.createDNASequenceWithErrorProbabilities(null);
            fail("createDNASequenceWithErrorProbabilities(null) expected IllegalArgumentException");
        }
        catch (IllegalArgumentException e)
        {
            // expected
        }
    }

    public void testCreateDNASequenceWithQualityScoresAndErrorProbabilities()
    {
        DNASequence sequence = FastqTools.createDNASequenceWithQualityScoresAndErrorProbabilities(builder.build());
        assertNotNull(sequence);

        List qualityScoresFeatures = sequence.getFeaturesByType("qualityScores");
        assertNotNull(qualityScoresFeatures);
        assertEquals(1, qualityScoresFeatures.size());
        QualityFeature qualityScores = (QualityFeature) qualityScoresFeatures.get(0);
        assertEquals(sequence.getLength(), qualityScores.getQualities().size());
        assertEquals(sequence.getLength(), qualityScores.getLocations().getLength());

        List errorProbabilitiesFeatures = sequence.getFeaturesByType("errorProbabilities");
        assertNotNull(errorProbabilitiesFeatures);
        assertEquals(1, errorProbabilitiesFeatures.size());
        QuantityFeature errorProbabilities = (QuantityFeature) errorProbabilitiesFeatures.get(0);
        assertEquals(sequence.getLength(), errorProbabilities.getQuantities().size());
        assertEquals(sequence.getLength(), errorProbabilities.getLocations().getLength());
    }

    public void testCreateDNASequenceWithQualityScoresAndErrorProbabilitiesNullFastq()
    {
        try
        {
            FastqTools.createDNASequenceWithQualityScoresAndErrorProbabilities(null);
            fail("createDNASequenceWithQualityScoresAndErrorProbabilities(null) expected IllegalArgumentException");
        }
        catch (IllegalArgumentException e)
        {
            // expected
        }
    }

    public void testCreateQualityScores()
    {
        Fastq fastq = builder.build();
        QualityFeature qualityScores = FastqTools.createQualityScores(fastq);
        assertNotNull(qualityScores);
        assertEquals(fastq.getSequence().length(), qualityScores.getQualities().size());
    }

    public void testCreateQualityScoresNullFastq()
    {
        try
        {
            FastqTools.createQualityScores(null);
            fail("createQualityScores(null) expected IllegalArgumentException");
        }
        catch (IllegalArgumentException e)
        {
            // expected
        }
    }

    public void testCreateErrorProbabilities()
    {
        Fastq fastq = builder.build();
        QuantityFeature errorProbabilities = FastqTools.createErrorProbabilities(fastq);
        assertNotNull(errorProbabilities);
        assertEquals(fastq.getSequence().length(), errorProbabilities.getQuantities().size());
    }

    public void testCreateErrorProbabilitiesNullFastq()
    {
        try
        {
            FastqTools.createErrorProbabilities(null);
            fail("createErrorProbabilities(null) expected IllegalArgumentException");
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

    public void testToList()
    {
        List<String> list = new ArrayList<String>();
        assertSame(list, FastqTools.toList(list));
    }

    public void testToListNotAList()
    {
        Collection<String> collection = new HashSet<String>();
        assertTrue(FastqTools.toList(collection) instanceof List);
        assertNotSame(collection, FastqTools.toList(collection));
    }
}

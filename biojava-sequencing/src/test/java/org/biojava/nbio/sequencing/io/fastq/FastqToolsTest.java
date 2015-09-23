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

import junit.framework.TestCase;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.features.FeatureInterface;
import org.biojava.nbio.core.sequence.features.QualityFeature;
import org.biojava.nbio.core.sequence.features.QuantityFeature;
import org.biojava.nbio.core.sequence.template.AbstractSequence;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;

/**
 * Unit test for FastqTools.
 */
public final class FastqToolsTest extends TestCase
{
    private final FastqBuilder builder = new FastqBuilder().withDescription("foo").withSequence("ACTG").withQuality("ZZZZ");

    public void testCreateDNASequence() throws CompoundNotFoundException
    {
        DNASequence sequence = FastqTools.createDNASequence(builder.build());
        assertNotNull(sequence);
    }

    public void testCreateDNASequenceNullFastq() throws CompoundNotFoundException
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

    public void testCreateDNASequenceWithQualityScores() throws CompoundNotFoundException
    {
        DNASequence sequence = FastqTools.createDNASequenceWithQualityScores(builder.build());
        assertNotNull(sequence);

        List<FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound>> features = sequence.getFeaturesByType("qualityScores");
        assertNotNull(features);
        assertEquals(1, features.size());
        QualityFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> qualityScores = (QualityFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound>) features.get(0);
        assertEquals(sequence.getLength(), qualityScores.getQualities().size());
        assertEquals(sequence.getLength(), qualityScores.getLocations().getLength());
    }

    public void testCreateDNASequenceWithQualityScoresNullFastq() throws CompoundNotFoundException
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

    public void testCreateDNASequenceWithErrorProbabilies() throws CompoundNotFoundException
    {
        DNASequence sequence = FastqTools.createDNASequenceWithErrorProbabilities(builder.build());
        assertNotNull(sequence);

        List<FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound>> features = sequence.getFeaturesByType("errorProbabilities");
        assertNotNull(features);
        assertEquals(1, features.size());
        QuantityFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> errorProbabilities = (QuantityFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound>) features.get(0);
        assertEquals(sequence.getLength(), errorProbabilities.getQuantities().size());
        assertEquals(sequence.getLength(), errorProbabilities.getLocations().getLength());
    }

    public void testCreateDNASequenceWithErrorProbabilitiesNullFastq() throws CompoundNotFoundException
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

    public void testCreateDNASequenceWithQualityScoresAndErrorProbabilities() throws CompoundNotFoundException
    {
        DNASequence sequence = FastqTools.createDNASequenceWithQualityScoresAndErrorProbabilities(builder.build());
        assertNotNull(sequence);

        List<FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound>> qualityScoresFeatures = sequence.getFeaturesByType("qualityScores");
        assertNotNull(qualityScoresFeatures);
        assertEquals(1, qualityScoresFeatures.size());
        QualityFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> qualityScores = (QualityFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound>) qualityScoresFeatures.get(0);
        assertEquals(sequence.getLength(), qualityScores.getQualities().size());
        assertEquals(sequence.getLength(), qualityScores.getLocations().getLength());

        List<FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound>> errorProbabilitiesFeatures = sequence.getFeaturesByType("errorProbabilities");
        assertNotNull(errorProbabilitiesFeatures);
        assertEquals(1, errorProbabilitiesFeatures.size());
        QuantityFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> errorProbabilities = (QuantityFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound>) errorProbabilitiesFeatures.get(0);
        assertEquals(sequence.getLength(), errorProbabilities.getQuantities().size());
        assertEquals(sequence.getLength(), errorProbabilities.getLocations().getLength());
    }

    public void testCreateDNASequenceWithQualityScoresAndErrorProbabilitiesNullFastq() throws CompoundNotFoundException
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
        QualityFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> qualityScores = FastqTools.createQualityScores(fastq);
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
        QuantityFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> errorProbabilities = FastqTools.createErrorProbabilities(fastq);
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
        Iterable<Number> qualityScores = FastqTools.qualityScores(builder.build());
        assertNotNull(qualityScores);
        int count = 0;
        for (Number qualityScore : qualityScores)
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
        Iterable<Number> errorProbabilities = FastqTools.errorProbabilities(builder.build());
        assertNotNull(errorProbabilities);
        int count = 0;
        for (Number errorProbability : errorProbabilities)
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

    public void testConvertNullFastq()
    {
        try
        {
            FastqTools.convert(null, FastqVariant.FASTQ_SANGER);
            fail("expected IllegalArgumentException");
        }
        catch (IllegalArgumentException e)
        {
            // expected
        }
    }

    public void testConvertNullVariant()
    {
        try
        {
            FastqTools.convert(builder.build(), null);
            fail("expected IllegalArgumentException");
        }
        catch (IllegalArgumentException e)
        {
            // expected
        }
    }

    public void testConvertSameVariant()
    {
        Fastq fastq = builder.build();
        assertEquals(fastq, FastqTools.convert(fastq, fastq.getVariant()));
    }

    public void testConvertQualitiesNullFastq()
    {
        try
        {
            FastqTools.convertQualities(null, FastqVariant.FASTQ_SANGER);
            fail("expected IllegalArgumentException");
        }
        catch (IllegalArgumentException e)
        {
            // expected
        }
    }

    public void testConvertQualitiesNullVariant()
    {
        try
        {
            FastqTools.convertQualities(builder.build(), null);
            fail("expected IllegalArgumentException");
        }
        catch (IllegalArgumentException e)
        {
            // expected
        }
    }

    public void testConvertQualitiesSameVariant()
    {
        Fastq fastq = builder.build();
        assertEquals(fastq.getQuality(), FastqTools.convertQualities(fastq, fastq.getVariant()));
    }

    public void testConvertQualitiesSangerToSolexa()
    {
        Fastq fastq = builder.build();
        assertEquals("yyyy", FastqTools.convertQualities(fastq, FastqVariant.FASTQ_SOLEXA));
    }

    public void testConvertQualitiesSangerToIllumina()
    {
        Fastq fastq = builder.build();
        assertEquals("yyyy", FastqTools.convertQualities(fastq, FastqVariant.FASTQ_ILLUMINA));
    }

    public void testConvertQualitiesSolexaToSanger()
    {
        Fastq fastq = builder.withVariant(FastqVariant.FASTQ_SOLEXA).build();
        assertEquals(";;;;", FastqTools.convertQualities(fastq, FastqVariant.FASTQ_SANGER));
    }

    public void testConvertQualitiesIlluminaToSanger()
    {
        Fastq fastq = builder.withVariant(FastqVariant.FASTQ_ILLUMINA).build();
        assertEquals(";;;;", FastqTools.convertQualities(fastq, FastqVariant.FASTQ_SANGER));
    }

    public void testConvertQualitiesSolexaToIllumina()
    {
        Fastq fastq = builder.withVariant(FastqVariant.FASTQ_SOLEXA).build();
        assertEquals("ZZZZ", FastqTools.convertQualities(fastq, FastqVariant.FASTQ_ILLUMINA));
    }

    public void testConvertQualitiesIlluminaToSolexa()
    {
        Fastq fastq = builder.withVariant(FastqVariant.FASTQ_ILLUMINA).build();
        assertEquals("ZZZZ", FastqTools.convertQualities(fastq, FastqVariant.FASTQ_SOLEXA));
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

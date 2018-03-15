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
package org.biojava.nbio.genome.io.fastq;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.features.FeatureInterface;
import org.biojava.nbio.core.sequence.features.QualityFeature;
import org.biojava.nbio.core.sequence.features.QuantityFeature;
import org.biojava.nbio.core.sequence.template.AbstractSequence;
import org.junit.Assert;
import org.junit.Test;


import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;

/**
 * Unit test for FastqTools.
 */
public final class FastqToolsTest {
	private final FastqBuilder builder = new FastqBuilder().withDescription("foo").withSequence("ACTG").withQuality("ZZZZ");

	@Test
	public void testCreateDNASequence() throws CompoundNotFoundException
	{
		DNASequence sequence = FastqTools.createDNASequence(builder.build());
		Assert.assertNotNull(sequence);
	}

	@Test
	public void testCreateDNASequenceNullFastq() throws CompoundNotFoundException
	{
		try
		{
			FastqTools.createDNASequence(null);
			Assert.fail("createDNASequence(null) expected IllegalArgumentException");
		}
		catch (IllegalArgumentException e)
		{
			// expected
		}
	}

	@Test
	public void testCreateDNASequenceWithQualityScores() throws CompoundNotFoundException
	{
		DNASequence sequence = FastqTools.createDNASequenceWithQualityScores(builder.build());
		Assert.assertNotNull(sequence);

		List<FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound>> features = sequence.getFeaturesByType("qualityScores");
		Assert.assertNotNull(features);
		Assert.assertEquals(1, features.size());
		QualityFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> qualityScores = (QualityFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound>) features.get(0);
		Assert.assertEquals(sequence.getLength(), qualityScores.getQualities().size());
		Assert.assertEquals(sequence.getLength(), qualityScores.getLocations().getLength());
	}

	@Test
	public void testCreateDNASequenceWithQualityScoresNullFastq() throws CompoundNotFoundException
	{
		try
		{
			FastqTools.createDNASequenceWithQualityScores(null);
			Assert.fail("createDNASequenceWithQualityScores(null) expected IllegalArgumentException");
		}
		catch (IllegalArgumentException e)
		{
			// expected
		}
	}

	@Test
	public void testCreateDNASequenceWithErrorProbabilies() throws CompoundNotFoundException
	{
		DNASequence sequence = FastqTools.createDNASequenceWithErrorProbabilities(builder.build());
		Assert.assertNotNull(sequence);

		List<FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound>> features = sequence.getFeaturesByType("errorProbabilities");
		Assert.assertNotNull(features);
		Assert.assertEquals(1, features.size());
		QuantityFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> errorProbabilities = (QuantityFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound>) features.get(0);
		Assert.assertEquals(sequence.getLength(), errorProbabilities.getQuantities().size());
		Assert.assertEquals(sequence.getLength(), errorProbabilities.getLocations().getLength());
	}

	@Test
	public void testCreateDNASequenceWithErrorProbabilitiesNullFastq() throws CompoundNotFoundException
	{
		try
		{
			FastqTools.createDNASequenceWithErrorProbabilities(null);
			Assert.fail("createDNASequenceWithErrorProbabilities(null) expected IllegalArgumentException");
		}
		catch (IllegalArgumentException e)
		{
			// expected
		}
	}

	@Test
	public void testCreateDNASequenceWithQualityScoresAndErrorProbabilities() throws CompoundNotFoundException
	{
		DNASequence sequence = FastqTools.createDNASequenceWithQualityScoresAndErrorProbabilities(builder.build());
		Assert.assertNotNull(sequence);

		List<FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound>> qualityScoresFeatures = sequence.getFeaturesByType("qualityScores");
		Assert.assertNotNull(qualityScoresFeatures);
		Assert.assertEquals(1, qualityScoresFeatures.size());
		QualityFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> qualityScores = (QualityFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound>) qualityScoresFeatures.get(0);
		Assert.assertEquals(sequence.getLength(), qualityScores.getQualities().size());
		Assert.assertEquals(sequence.getLength(), qualityScores.getLocations().getLength());

		List<FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound>> errorProbabilitiesFeatures = sequence.getFeaturesByType("errorProbabilities");
		Assert.assertNotNull(errorProbabilitiesFeatures);
		Assert.assertEquals(1, errorProbabilitiesFeatures.size());
		QuantityFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> errorProbabilities = (QuantityFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound>) errorProbabilitiesFeatures.get(0);
		Assert.assertEquals(sequence.getLength(), errorProbabilities.getQuantities().size());
		Assert.assertEquals(sequence.getLength(), errorProbabilities.getLocations().getLength());
	}

	@Test
	public void testCreateDNASequenceWithQualityScoresAndErrorProbabilitiesNullFastq() throws CompoundNotFoundException
	{
		try
		{
			FastqTools.createDNASequenceWithQualityScoresAndErrorProbabilities(null);
			Assert.fail("createDNASequenceWithQualityScoresAndErrorProbabilities(null) expected IllegalArgumentException");
		}
		catch (IllegalArgumentException e)
		{
			// expected
		}
	}

	@Test
	public void testCreateQualityScores()
	{
		Fastq fastq = builder.build();
		QualityFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> qualityScores = FastqTools.createQualityScores(fastq);
		Assert.assertNotNull(qualityScores);
		Assert.assertEquals(fastq.getSequence().length(), qualityScores.getQualities().size());
	}

	@Test
	public void testCreateQualityScoresNullFastq()
	{
		try
		{
			FastqTools.createQualityScores(null);
			Assert.fail("createQualityScores(null) expected IllegalArgumentException");
		}
		catch (IllegalArgumentException e)
		{
			// expected
		}
	}

	@Test
	public void testCreateErrorProbabilities()
	{
		Fastq fastq = builder.build();
		QuantityFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> errorProbabilities = FastqTools.createErrorProbabilities(fastq);
		Assert.assertNotNull(errorProbabilities);
		Assert.assertEquals(fastq.getSequence().length(), errorProbabilities.getQuantities().size());
	}

	@Test
	public void testCreateErrorProbabilitiesNullFastq()
	{
		try
		{
			FastqTools.createErrorProbabilities(null);
			Assert.fail("createErrorProbabilities(null) expected IllegalArgumentException");
		}
		catch (IllegalArgumentException e)
		{
			// expected
		}
	}

	@Test
	public void testQualityScores()
	{
		Iterable<Number> qualityScores = FastqTools.qualityScores(builder.build());
		Assert.assertNotNull(qualityScores);
		int count = 0;
		for (Number qualityScore : qualityScores)
		{
			Assert.assertNotNull(qualityScore);
			count++;
		}
		Assert.assertEquals(4, count);
	}

	@Test
	public void testQualityScoresNullFastq()
	{
		try
		{
			FastqTools.qualityScores(null);
			Assert.fail("qualityScores(null) expected IllegalArgumentException");
		}
		catch (IllegalArgumentException e)
		{
			// expected
		}
	}

	@Test
	public void testQualityScoresIntArray()
	{
		int[] qualityScores = new int[4];
		FastqTools.qualityScores(builder.build(), qualityScores);
		for (int i = 0; i < 4; i++)
		{
			Assert.assertTrue(qualityScores[i] != 0);
		}
	}

	@Test
	public void testQualityScoresIntArrayNullFastq()
	{
		try
		{
			FastqTools.qualityScores(null, new int[0]);
			Assert.fail("qualityScores(null, int[]) expected IllegalArgumentException");
		}
		catch (IllegalArgumentException e)
		{
			// expected
		}
	}

	@Test
	public void testQualityScoresNullIntArray()
	{
		try
		{
			FastqTools.qualityScores(builder.build(), null);
			Assert.fail("qualityScores(fastq, null) expected IllegalArgumentException");
		}
		catch (IllegalArgumentException e)
		{
			// expected
		}
	}

	@Test
	public void testQualityScoresQualityScoresTooSmall()
	{
		try
		{
			FastqTools.qualityScores(builder.build(), new int[3]);
			Assert.fail("expected IllegalArgumentException");
		}
		catch (IllegalArgumentException e)
		{
			// expected
		}
	}

	@Test
	public void testQualityScoresQualityScoresTooLarge()
	{
		try
		{
			FastqTools.qualityScores(builder.build(), new int[5]);
			Assert.fail("expected IllegalArgumentException");
		}
		catch (IllegalArgumentException e)
		{
			// expected
		}
	}

	@Test
	public void testErrorProbabilities()
	{
		Iterable<Number> errorProbabilities = FastqTools.errorProbabilities(builder.build());
		Assert.assertNotNull(errorProbabilities);
		int count = 0;
		for (Number errorProbability : errorProbabilities)
		{
			Assert.assertNotNull(errorProbability);
			count++;
		}
		Assert.assertEquals(4, count);
	}

	@Test
	public void testErrorProbabilitiesNullFastq()
	{
		try
		{
			FastqTools.errorProbabilities(null);
			Assert.fail("errorProbabilities(null) expected IllegalArgumentException");
		}
		catch (IllegalArgumentException e)
		{
			// expected
		}
	}

	@Test
	public void testErrorProbabilitiesDoubleArray()
	{
		double[] errorProbabilities = new double[4];
		FastqTools.errorProbabilities(builder.build(), errorProbabilities);
		for (int i = 0; i < 0; i++)
		{
			Assert.assertTrue(errorProbabilities[i] > 0.0d);
		}
	}

	@Test
	public void testErrorProbabilitiesDoubleArrayNullFastq()
	{
		try
		{
			FastqTools.errorProbabilities(null, new double[0]);
			Assert.fail("errorProbabilities(null, double[]) expected IllegalArgumentException");
		}
		catch (IllegalArgumentException e)
		{
			// expected
		}
	}

	@Test
	public void testErrorProbabilitiesNullErrorProbabilities()
	{
		try
		{
			FastqTools.errorProbabilities(builder.build(), null);
			Assert.fail("errorProbabilities(fastq, null) expected IllegalArgumentException");
		}
		catch (IllegalArgumentException e)
		{
			// expected
		}
	}

	@Test
	public void testErrorProbabilitiesErrorProbabilitiesTooSmall()
	{
		try
		{
			FastqTools.errorProbabilities(builder.build(), new double[3]);
			Assert.fail("expected IllegalArgumentException");
		}
		catch (IllegalArgumentException e)
		{
			// expected
		}
	}

	@Test
	public void testErrorProbabilitiesErrorProbabilitiesTooLarge()
	{
		try
		{
			FastqTools.errorProbabilities(builder.build(), new double[5]);
			Assert.fail("expected IllegalArgumentException");
		}
		catch (IllegalArgumentException e)
		{
			// expected
		}
	}

	@Test
	public void testConvertNullFastq()
	{
		try
		{
			FastqTools.convert(null, FastqVariant.FASTQ_SANGER);
			Assert.fail("expected IllegalArgumentException");
		}
		catch (IllegalArgumentException e)
		{
			// expected
		}
	}

	@Test
	public void testConvertNullVariant()
	{
		try
		{
			FastqTools.convert(builder.build(), null);
			Assert.fail("expected IllegalArgumentException");
		}
		catch (IllegalArgumentException e)
		{
			// expected
		}
	}

	@Test
	public void testConvertSameVariant()
	{
		Fastq fastq = builder.build();
		Assert.assertEquals(fastq, FastqTools.convert(fastq, fastq.getVariant()));
	}

	@Test
	public void testConvertQualitiesNullFastq()
	{
		try
		{
			FastqTools.convertQualities(null, FastqVariant.FASTQ_SANGER);
			Assert.fail("expected IllegalArgumentException");
		}
		catch (IllegalArgumentException e)
		{
			// expected
		}
	}

	@Test
	public void testConvertQualitiesNullVariant()
	{
		try
		{
			FastqTools.convertQualities(builder.build(), null);
			Assert.fail("expected IllegalArgumentException");
		}
		catch (IllegalArgumentException e)
		{
			// expected
		}
	}

	@Test
	public void testConvertQualitiesSameVariant()
	{
		Fastq fastq = builder.build();
		Assert.assertEquals(fastq.getQuality(), FastqTools.convertQualities(fastq, fastq.getVariant()));
	}

	@Test
	public void testConvertQualitiesSangerToSolexa()
	{
		Fastq fastq = builder.build();
		Assert.assertEquals("yyyy", FastqTools.convertQualities(fastq, FastqVariant.FASTQ_SOLEXA));
	}

	@Test
	public void testConvertQualitiesSangerToIllumina()
	{
		Fastq fastq = builder.build();
		Assert.assertEquals("yyyy", FastqTools.convertQualities(fastq, FastqVariant.FASTQ_ILLUMINA));
	}

	@Test
	public void testConvertQualitiesSolexaToSanger()
	{
		Fastq fastq = builder.withVariant(FastqVariant.FASTQ_SOLEXA).build();
		Assert.assertEquals(";;;;", FastqTools.convertQualities(fastq, FastqVariant.FASTQ_SANGER));
	}

	@Test
	public void testConvertQualitiesIlluminaToSanger()
	{
		Fastq fastq = builder.withVariant(FastqVariant.FASTQ_ILLUMINA).build();
		Assert.assertEquals(";;;;", FastqTools.convertQualities(fastq, FastqVariant.FASTQ_SANGER));
	}

	@Test
	public void testConvertQualitiesSolexaToIllumina()
	{
		Fastq fastq = builder.withVariant(FastqVariant.FASTQ_SOLEXA).build();
		Assert.assertEquals("ZZZZ", FastqTools.convertQualities(fastq, FastqVariant.FASTQ_ILLUMINA));
	}

	@Test
	public void testConvertQualitiesIlluminaToSolexa()
	{
		Fastq fastq = builder.withVariant(FastqVariant.FASTQ_ILLUMINA).build();
		Assert.assertEquals("ZZZZ", FastqTools.convertQualities(fastq, FastqVariant.FASTQ_SOLEXA));
	}

	@Test
	public void testToList()
	{
		List<String> list = new ArrayList<String>();
		Assert.assertSame(list, FastqTools.toList(list));
	}

	@Test
	public void testToListNotAList()
	{
		Collection<String> collection = new HashSet<String>();
		Assert.assertTrue(FastqTools.toList(collection) instanceof List);
		Assert.assertNotSame(collection, FastqTools.toList(collection));
	}
}

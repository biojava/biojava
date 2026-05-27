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
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

/**
 * Unit test for FastqTools.
 */
final class FastqToolsTest {
	private final FastqBuilder builder = new FastqBuilder().withDescription("foo").withSequence("ACTG").withQuality("ZZZZ");

	@Test
    void testCreateDNASequence() throws CompoundNotFoundException
	{
		DNASequence sequence = FastqTools.createDNASequence(builder.build());
		Assertions.assertNotNull(sequence);
	}

	@Test
    void testCreateDNASequenceNullFastq()
	{
		Assertions.assertThrows(IllegalArgumentException.class, () -> FastqTools.createDNASequence(null));
	}

	@Test
    void testCreateDNASequenceWithQualityScores() throws CompoundNotFoundException
	{
		DNASequence sequence = FastqTools.createDNASequenceWithQualityScores(builder.build());
		Assertions.assertNotNull(sequence);

		List<FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound>> features = sequence.getFeaturesByType("qualityScores");
		Assertions.assertNotNull(features);
		Assertions.assertEquals(1, features.size());
		QualityFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> qualityScores = (QualityFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound>) features.get(0);
		Assertions.assertEquals(sequence.getLength(), qualityScores.getQualities().size());
		Assertions.assertEquals(sequence.getLength(), qualityScores.getLocations().getLength());
	}

	@Test
    void testCreateDNASequenceWithQualityScoresNullFastq()
	{
		Assertions.assertThrows(IllegalArgumentException.class, () -> FastqTools.createDNASequenceWithQualityScores(null));
	}

	@Test
    void testCreateDNASequenceWithErrorProbabilies() throws CompoundNotFoundException
	{
		DNASequence sequence = FastqTools.createDNASequenceWithErrorProbabilities(builder.build());
		Assertions.assertNotNull(sequence);

		List<FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound>> features = sequence.getFeaturesByType("errorProbabilities");
		Assertions.assertNotNull(features);
		Assertions.assertEquals(1, features.size());
		QuantityFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> errorProbabilities = (QuantityFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound>) features.get(0);
		Assertions.assertEquals(sequence.getLength(), errorProbabilities.getQuantities().size());
		Assertions.assertEquals(sequence.getLength(), errorProbabilities.getLocations().getLength());
	}

	@Test
    void testCreateDNASequenceWithErrorProbabilitiesNullFastq()
	{
		Assertions.assertThrows(IllegalArgumentException.class, () -> FastqTools.createDNASequenceWithErrorProbabilities(null));
	}

	@Test
    void testCreateDNASequenceWithQualityScoresAndErrorProbabilities() throws CompoundNotFoundException
	{
		DNASequence sequence = FastqTools.createDNASequenceWithQualityScoresAndErrorProbabilities(builder.build());
		Assertions.assertNotNull(sequence);

		List<FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound>> qualityScoresFeatures = sequence.getFeaturesByType("qualityScores");
		Assertions.assertNotNull(qualityScoresFeatures);
		Assertions.assertEquals(1, qualityScoresFeatures.size());
		QualityFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> qualityScores = (QualityFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound>) qualityScoresFeatures.get(0);
		Assertions.assertEquals(sequence.getLength(), qualityScores.getQualities().size());
		Assertions.assertEquals(sequence.getLength(), qualityScores.getLocations().getLength());

		List<FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound>> errorProbabilitiesFeatures = sequence.getFeaturesByType("errorProbabilities");
		Assertions.assertNotNull(errorProbabilitiesFeatures);
		Assertions.assertEquals(1, errorProbabilitiesFeatures.size());
		QuantityFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> errorProbabilities = (QuantityFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound>) errorProbabilitiesFeatures.get(0);
		Assertions.assertEquals(sequence.getLength(), errorProbabilities.getQuantities().size());
		Assertions.assertEquals(sequence.getLength(), errorProbabilities.getLocations().getLength());
	}

	@Test
    void testCreateDNASequenceWithQualityScoresAndErrorProbabilitiesNullFastq()
	{
		Assertions.assertThrows(IllegalArgumentException.class, () -> FastqTools.createDNASequenceWithQualityScoresAndErrorProbabilities(null));
	}

	@Test
    void testCreateQualityScores()
	{
		Fastq fastq = builder.build();
		QualityFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> qualityScores = FastqTools.createQualityScores(fastq);
		Assertions.assertNotNull(qualityScores);
		Assertions.assertEquals(fastq.getSequence().length(), qualityScores.getQualities().size());
	}

	@Test
    void testCreateQualityScoresNullFastq()
	{
		Assertions.assertThrows(IllegalArgumentException.class, () -> FastqTools.createQualityScores(null));
	}

	@Test
    void testCreateErrorProbabilities()
	{
		Fastq fastq = builder.build();
		QuantityFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> errorProbabilities = FastqTools.createErrorProbabilities(fastq);
		Assertions.assertNotNull(errorProbabilities);
		Assertions.assertEquals(fastq.getSequence().length(), errorProbabilities.getQuantities().size());
	}

	@Test
    void testCreateErrorProbabilitiesNullFastq()
	{
		Assertions.assertThrows(IllegalArgumentException.class, () -> FastqTools.createErrorProbabilities(null));
	}

	@Test
    void testQualityScores()
	{
		Iterable<Number> qualityScores = FastqTools.qualityScores(builder.build());
		List<Number> scoresList = StreamSupport.stream(qualityScores.spliterator(), false)
				.collect(Collectors.toList());
		Assertions.assertAll(
				() -> Assertions.assertEquals(4, scoresList.size()),
				() -> Assertions.assertFalse(scoresList.contains(null))
		);
	}

	@Test
    void testQualityScoresNullFastq()
	{
		Assertions.assertThrows(IllegalArgumentException.class, () -> FastqTools.qualityScores(null));
	}

	@Test
    void testQualityScoresIntArray()
	{
		int[] qualityScores = new int[4];
		FastqTools.qualityScores(builder.build(), qualityScores);

		Assertions.assertTrue(Arrays.stream(qualityScores).allMatch(score -> score != 0), () ->
				"Array contains zero at some position: " + Arrays.toString(qualityScores));
	}

	@Test
    void testQualityScoresIntArrayNullFastq()
	{
		Assertions.assertThrows(IllegalArgumentException.class, () -> FastqTools.qualityScores(null, new int[0]));
	}

	@Test
    void testQualityScoresNullIntArray()
	{
		Fastq fastq = builder.build();
		Assertions.assertThrows(IllegalArgumentException.class, () -> FastqTools.qualityScores(fastq, null));
	}

	@Test
    void testQualityScoresQualityScoresTooSmall()
	{
		Fastq fastq = builder.build();
		Assertions.assertThrows(IllegalArgumentException.class, () -> FastqTools.qualityScores(fastq, new int[3]));
	}

	@Test
    void testQualityScoresQualityScoresTooLarge()
	{
		Fastq fastq = builder.build();
		Assertions.assertThrows(IllegalArgumentException.class, () -> FastqTools.qualityScores(fastq, new int[5]));
	}

	@Test
    void testErrorProbabilities()
	{
		Iterable<Number> errorProbabilities = FastqTools.errorProbabilities(builder.build());
		List<Number> scores = StreamSupport.stream(errorProbabilities.spliterator(), false)
				.collect(Collectors.toList());

		Assertions.assertNotNull(scores);
		Assertions.assertEquals(4, scores.size());
		Assertions.assertTrue(scores.stream().allMatch(Objects::nonNull));
	}

	@Test
    void testErrorProbabilitiesNullFastq()
	{
		Assertions.assertThrows(IllegalArgumentException.class, () -> FastqTools.errorProbabilities(null));
	}

	@Test
    void testErrorProbabilitiesDoubleArray()
	{
		double[] errorProbabilities = new double[4];
		FastqTools.errorProbabilities(builder.build(), errorProbabilities);
		Assertions.assertTrue(
				Arrays.stream(errorProbabilities).allMatch(p -> p > 0.0),
				() -> "Expected all probabilities to be > 0.0, but got: " + Arrays.toString(errorProbabilities)
		);
	}

	@Test
    void testErrorProbabilitiesDoubleArrayNullFastq()
	{
		Assertions.assertThrows(IllegalArgumentException.class, () -> FastqTools.errorProbabilities(null, new double[0]));
	}

	@Test
    void testErrorProbabilitiesNullErrorProbabilities()
	{
		Fastq fastq = builder.build();
		Assertions.assertThrows(IllegalArgumentException.class, () -> FastqTools.errorProbabilities(fastq, null));
	}

	@Test
    void testErrorProbabilitiesErrorProbabilitiesTooSmall()
	{
		Fastq fastq = builder.build();
		Assertions.assertThrows(IllegalArgumentException.class, () -> FastqTools.errorProbabilities(fastq, new double[3]));
	}

	@Test
    void testErrorProbabilitiesErrorProbabilitiesTooLarge()
	{
		Fastq fastq = builder.build();
		Assertions.assertThrows(IllegalArgumentException.class, () -> FastqTools.errorProbabilities(fastq, new double[5]));
	}

	@Test
    void testConvertNullFastq()
	{
		Assertions.assertThrows(IllegalArgumentException.class, () -> FastqTools.convert(null, FastqVariant.FASTQ_SANGER));
	}

	@Test
    void testConvertNullVariant()
	{
		Fastq fastq = builder.build();
		Assertions.assertThrows(IllegalArgumentException.class, () -> FastqTools.convert(fastq, null));
	}

	@Test
    void testConvertSameVariant()
	{
		Fastq fastq = builder.build();
		Assertions.assertEquals(fastq, FastqTools.convert(fastq, fastq.getVariant()));
	}

	@Test
    void testConvertQualitiesNullFastq()
	{
		Assertions.assertThrows(IllegalArgumentException.class, () -> FastqTools.convertQualities(null, FastqVariant.FASTQ_SANGER));
	}

	@Test
    void testConvertQualitiesNullVariant()
	{
		Fastq fastq = builder.build();
		Assertions.assertThrows(IllegalArgumentException.class, () -> FastqTools.convertQualities(fastq, null));
	}

	@Test
    void testConvertQualitiesSameVariant()
	{
		Fastq fastq = builder.build();
		Assertions.assertEquals(fastq.getQuality(), FastqTools.convertQualities(fastq, fastq.getVariant()));
	}

	@Test
    void testConvertQualitiesSangerToSolexa()
	{
		Fastq fastq = builder.build();
		Assertions.assertEquals("yyyy", FastqTools.convertQualities(fastq, FastqVariant.FASTQ_SOLEXA));
	}

	@Test
    void testConvertQualitiesSangerToIllumina()
	{
		Fastq fastq = builder.build();
		Assertions.assertEquals("yyyy", FastqTools.convertQualities(fastq, FastqVariant.FASTQ_ILLUMINA));
	}

	@Test
    void testConvertQualitiesSolexaToSanger()
	{
		Fastq fastq = builder.withVariant(FastqVariant.FASTQ_SOLEXA).build();
		Assertions.assertEquals(";;;;", FastqTools.convertQualities(fastq, FastqVariant.FASTQ_SANGER));
	}

	@Test
    void testConvertQualitiesIlluminaToSanger()
	{
		Fastq fastq = builder.withVariant(FastqVariant.FASTQ_ILLUMINA).build();
		Assertions.assertEquals(";;;;", FastqTools.convertQualities(fastq, FastqVariant.FASTQ_SANGER));
	}

	@Test
    void testConvertQualitiesSolexaToIllumina()
	{
		Fastq fastq = builder.withVariant(FastqVariant.FASTQ_SOLEXA).build();
		Assertions.assertEquals("ZZZZ", FastqTools.convertQualities(fastq, FastqVariant.FASTQ_ILLUMINA));
	}

	@Test
    void testConvertQualitiesIlluminaToSolexa()
	{
		Fastq fastq = builder.withVariant(FastqVariant.FASTQ_ILLUMINA).build();
		Assertions.assertEquals("ZZZZ", FastqTools.convertQualities(fastq, FastqVariant.FASTQ_SOLEXA));
	}

	@Test
    void testToList()
	{
		List<String> list = new ArrayList<>();
		Assertions.assertSame(list, FastqTools.toList(list));
	}

	@Test
    void testToListNotAList()
	{
		Collection<String> collection = new HashSet<>();
		Assertions.assertTrue(FastqTools.toList(collection) instanceof List);
		Assertions.assertNotSame(collection, FastqTools.toList(collection));
	}
}

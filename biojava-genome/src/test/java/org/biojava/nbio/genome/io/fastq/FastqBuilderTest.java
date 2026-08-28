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

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

/**
 * Unit test for FastqBuilder.
 */
final class FastqBuilderTest {

	@Test
    void testConstructor()
	{
		FastqBuilder fastqBuilder = new FastqBuilder();
		Assertions.assertNotNull(fastqBuilder);
	}

	@Test
    void testConstructorFastq()
	{
		FastqBuilder fastqBuilder = new FastqBuilder()
			.withDescription("description")
			.withSequence("sequence")
			.withQuality("quality_")
			.withVariant(FastqVariant.FASTQ_SOLEXA);

		Fastq fastq = fastqBuilder.build();

		FastqBuilder fastqBuilder2 = new FastqBuilder(fastq);
		Assertions.assertNotNull(fastqBuilder2);

		Fastq fastq2 = fastqBuilder2.build();
		Assertions.assertEquals("description", fastq2.getDescription());
		Assertions.assertEquals("sequence", fastq2.getSequence());
		Assertions.assertEquals("quality_", fastq2.getQuality());
		Assertions.assertEquals(FastqVariant.FASTQ_SOLEXA, fastq2.getVariant());
	}

	@Test
    void testConstructorNullFastq()
	{
		Assertions.assertThrows(IllegalArgumentException.class, () -> new FastqBuilder(null));
	}

	@Test
    void testBuildDefault()
	{
		try
		{
			FastqBuilder fastqBuilder = new FastqBuilder();
			fastqBuilder.build();
			Assertions.fail("build default expected IllegalStateException");
		}
		catch (IllegalStateException e)
		{
			// expected
		}
	}

	@Test
    void testBuildNullDescription()
	{
		try
		{
			FastqBuilder fastqBuilder = new FastqBuilder()
				.withDescription(null)
				.withSequence("sequence")
				.withQuality("quality_")
				.withVariant(FastqVariant.FASTQ_SOLEXA);

			fastqBuilder.build();
			Assertions.fail("build null description expected IllegalArgumentException");
		}
		catch (IllegalArgumentException e)
		{
			// expected
		}
	}

	@Test
    void testBuildNullSequence()
	{
		try
		{
			FastqBuilder fastqBuilder = new FastqBuilder()
				.withDescription("description")
				.withSequence(null)
				.withQuality("quality_")
				.withVariant(FastqVariant.FASTQ_SOLEXA);

			fastqBuilder.build();
			Assertions.fail("build null sequence expected IllegalArgumentException");
		}
		catch (IllegalArgumentException e)
		{
			// expected
		}
	}

	@Test
    void testBuildNullAppendSequence()
	{
		try
		{
			FastqBuilder fastqBuilder = new FastqBuilder()
				.withDescription("description")
				.appendSequence(null)
				.withQuality("quality_")
				.withVariant(FastqVariant.FASTQ_SOLEXA);

			fastqBuilder.build();
			Assertions.fail("build null append sequence expected IllegalArgumentException");
		}
		catch (IllegalArgumentException e)
		{
			// expected
		}
	}

	@Test
    void testBuildNullQuality()
	{
		try
		{
			FastqBuilder fastqBuilder = new FastqBuilder()
				.withDescription("description")
				.withSequence("sequence")
				.withQuality(null)
				.withVariant(FastqVariant.FASTQ_SOLEXA);

			fastqBuilder.build();
			Assertions.fail("build null quality expected IllegalArgumentException");
		}
		catch (IllegalArgumentException e)
		{
			// expected
		}
	}

	@Test
    void testBuildNullAppendQuality()
	{
		try
		{
			FastqBuilder fastqBuilder = new FastqBuilder()
				.withDescription("description")
				.withSequence("sequence")
				.appendQuality(null)
				.withVariant(FastqVariant.FASTQ_SOLEXA);

			fastqBuilder.build();
			Assertions.fail("build null append quality expected IllegalArgumentException");
		}
		catch (IllegalArgumentException e)
		{
			// expected
		}
	}

	@Test
    void testBuildNullVariant()
	{
		try
		{
			FastqBuilder fastqBuilder = new FastqBuilder()
				.withDescription("description")
				.withSequence("sequence")
				.withQuality("quality_")
				.withVariant(null);

			fastqBuilder.build();
			Assertions.fail("build null variant expected IllegalArgumentException");
		}
		catch (IllegalArgumentException e)
		{
			// expected
		}
	}

	@Test
    void testBuildMissingDescription()
	{
		try
		{
			FastqBuilder fastqBuilder = new FastqBuilder()
				.withSequence("sequence")
				.withQuality("quality_")
				.withVariant(FastqVariant.FASTQ_SOLEXA);

			fastqBuilder.build();
			Assertions.fail("build missing description expected IllegalStateException");
		}
		catch (IllegalStateException e)
		{
			// expected
		}
	}

	@Test
    void testBuildMissingSequence()
	{
		try
		{
			FastqBuilder fastqBuilder = new FastqBuilder()
				.withDescription("description")
				.withQuality("quality_")
				.withVariant(FastqVariant.FASTQ_SOLEXA);

			fastqBuilder.build();
			Assertions.fail("build missing sequence expected IllegalStateException");
		}
		catch (IllegalStateException e)
		{
			// expected
		}
	}

	@Test
    void testBuildMissingQuality()
	{
		try
		{
			FastqBuilder fastqBuilder = new FastqBuilder()
				.withDescription("description")
				.withSequence("sequence")
				.withVariant(FastqVariant.FASTQ_SOLEXA);

			fastqBuilder.build();
			Assertions.fail("build missing quality expected IllegalStateException");
		}
		catch (IllegalStateException e)
		{
			// expected
		}
	}

	@Test
    void testBuildDefaultVariant()
	{
		FastqBuilder fastqBuilder = new FastqBuilder()
			.withDescription("description")
			.withSequence("sequence")
			.withQuality("quality_");

		Fastq fastq = fastqBuilder.build();
		Assertions.assertEquals("description", fastqBuilder.getDescription());
		Assertions.assertTrue(fastqBuilder.sequenceAndQualityLengthsMatch());
		Assertions.assertEquals("description", fastq.getDescription());
		Assertions.assertEquals("sequence", fastq.getSequence());
		Assertions.assertEquals("quality_", fastq.getQuality());
		Assertions.assertEquals(FastqBuilder.DEFAULT_VARIANT, fastq.getVariant());
	}

	@Test
    void testBuild()
	{
		FastqBuilder fastqBuilder = new FastqBuilder()
			.withDescription("description")
			.withSequence("sequence")
			.withQuality("quality_")
			.withVariant(FastqVariant.FASTQ_SOLEXA);
		Fastq fastq = fastqBuilder.build();
		Assertions.assertEquals("description", fastqBuilder.getDescription());
		Assertions.assertTrue(fastqBuilder.sequenceAndQualityLengthsMatch());
		Assertions.assertEquals("description", fastq.getDescription());
		Assertions.assertEquals("sequence", fastq.getSequence());
		Assertions.assertEquals("quality_", fastq.getQuality());
		Assertions.assertEquals(FastqVariant.FASTQ_SOLEXA, fastq.getVariant());
	}

	@Test
    void testBuildAppendSequence()
	{
		FastqBuilder fastqBuilder = new FastqBuilder()
			.withDescription("description")
			.appendSequence("seq")
			.appendSequence("uence")
			.withQuality("quality_")
			.withVariant(FastqVariant.FASTQ_SOLEXA);
		Fastq fastq = fastqBuilder.build();
		Assertions.assertEquals("description", fastqBuilder.getDescription());
		Assertions.assertTrue(fastqBuilder.sequenceAndQualityLengthsMatch());
		Assertions.assertEquals("description", fastq.getDescription());
		Assertions.assertEquals("sequence", fastq.getSequence());
		Assertions.assertEquals("quality_", fastq.getQuality());
		Assertions.assertEquals(FastqVariant.FASTQ_SOLEXA, fastq.getVariant());
	}

	@Test
    void testBuildAppendQuality()
	{
		FastqBuilder fastqBuilder = new FastqBuilder()
			.withDescription("description")
			.withSequence("sequence")
			.appendQuality("qual")
			.appendQuality("ity_")
			.withVariant(FastqVariant.FASTQ_SOLEXA);
		Fastq fastq = fastqBuilder.build();
		Assertions.assertEquals("description", fastqBuilder.getDescription());
		Assertions.assertTrue(fastqBuilder.sequenceAndQualityLengthsMatch());
		Assertions.assertEquals("description", fastq.getDescription());
		Assertions.assertEquals("sequence", fastq.getSequence());
		Assertions.assertEquals("quality_", fastq.getQuality());
		Assertions.assertEquals(FastqVariant.FASTQ_SOLEXA, fastq.getVariant());
	}

	@Test
    void testBuildNonMatchingSequenceQualityScoreLengthsBothNull()
	{
		FastqBuilder fastqBuilder = new FastqBuilder()
			.withDescription("description")
			.withVariant(FastqVariant.FASTQ_SOLEXA);

		Assertions.assertTrue(fastqBuilder.sequenceAndQualityLengthsMatch());
	}

	@Test
    void testBuildNonMatchingSequenceQualityScoreLengthsSequenceNull()
	{
		FastqBuilder fastqBuilder = new FastqBuilder()
			.withDescription("description")
			.withQuality("0123")
			.withVariant(FastqVariant.FASTQ_SOLEXA);

		Assertions.assertEquals(false, fastqBuilder.sequenceAndQualityLengthsMatch());
	}

	@Test
    void testBuildNonMatchingSequenceQualityScoreLengthsQualityNull()
	{
		FastqBuilder fastqBuilder = new FastqBuilder()
			.withDescription("description")
			.withSequence("ACTG")
			.withVariant(FastqVariant.FASTQ_SOLEXA);

		Assertions.assertEquals(false, fastqBuilder.sequenceAndQualityLengthsMatch());
	}

	@Test
    void testBuildNonMatchingSequenceQualityScoreLengths0()
	{
		try
		{
			FastqBuilder fastqBuilder = new FastqBuilder()
				.withDescription("description")
				.withSequence("01234")
				.withQuality("0123")
				.withVariant(FastqVariant.FASTQ_SOLEXA);

			fastqBuilder.build();
			Assertions.fail("build sequence length > quality length expected IllegalStateException");
		}
		catch (IllegalStateException e)
		{
			// expected
		}
	}

	@Test
    void testBuildNonMatchingSequenceQualityScoreLengths1()
	{
		try
		{
			FastqBuilder fastqBuilder = new FastqBuilder()
				.withDescription("description")
				.withSequence("0123")
				.withQuality("01234")
				.withVariant(FastqVariant.FASTQ_SOLEXA);

			fastqBuilder.build();
			Assertions.fail("build sequence length < quality length expected IllegalStateException");
		}
		catch (IllegalStateException e)
		{
			// expected
		}
	}

	@Test
    void testBuildMultiple()
	{
		FastqBuilder fastqBuilder = new FastqBuilder()
			.withDescription("description")
			.withQuality("quality__")
			.withVariant(FastqVariant.FASTQ_SOLEXA);

		for (int i = 0; i < 10; i++)
		{
			Fastq fastq = fastqBuilder.withSequence("sequence" + i).build();
			Assertions.assertEquals("description", fastqBuilder.getDescription());
			Assertions.assertTrue(fastqBuilder.sequenceAndQualityLengthsMatch());
			Assertions.assertEquals("description", fastq.getDescription());
			Assertions.assertEquals("sequence" + i, fastq.getSequence());
			Assertions.assertEquals("quality__", fastq.getQuality());
			Assertions.assertEquals(FastqVariant.FASTQ_SOLEXA, fastq.getVariant());
		}
	}
}

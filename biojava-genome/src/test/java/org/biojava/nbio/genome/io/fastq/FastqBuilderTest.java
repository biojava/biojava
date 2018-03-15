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


import org.junit.Assert;
import org.junit.Test;

/**
 * Unit test for FastqBuilder.
 */
public final class FastqBuilderTest {

	@Test
	public void testConstructor()
	{
		FastqBuilder fastqBuilder = new FastqBuilder();
		Assert.assertNotNull(fastqBuilder);
	}

	@Test
	public void testBuildDefault()
	{
		try
		{
			FastqBuilder fastqBuilder = new FastqBuilder();
			fastqBuilder.build();
			Assert.fail("build default expected IllegalStateException");
		}
		catch (IllegalStateException e)
		{
			// expected
		}
	}

	@Test
	public void testBuildNullDescription()
	{
		try
		{
			FastqBuilder fastqBuilder = new FastqBuilder()
				.withDescription(null)
				.withSequence("sequence")
				.withQuality("quality_")
				.withVariant(FastqVariant.FASTQ_SOLEXA);

			fastqBuilder.build();
			Assert.fail("build null description expected IllegalArgumentException");
		}
		catch (IllegalArgumentException e)
		{
			// expected
		}
	}

	@Test
	public void testBuildNullSequence()
	{
		try
		{
			FastqBuilder fastqBuilder = new FastqBuilder()
				.withDescription("description")
				.withSequence(null)
				.withQuality("quality_")
				.withVariant(FastqVariant.FASTQ_SOLEXA);

			fastqBuilder.build();
			Assert.fail("build null sequence expected IllegalArgumentException");
		}
		catch (IllegalArgumentException e)
		{
			// expected
		}
	}

	@Test
	public void testBuildNullAppendSequence()
	{
		try
		{
			FastqBuilder fastqBuilder = new FastqBuilder()
				.withDescription("description")
				.appendSequence(null)
				.withQuality("quality_")
				.withVariant(FastqVariant.FASTQ_SOLEXA);

			fastqBuilder.build();
			Assert.fail("build null append sequence expected IllegalArgumentException");
		}
		catch (IllegalArgumentException e)
		{
			// expected
		}
	}

	@Test
	public void testBuildNullQuality()
	{
		try
		{
			FastqBuilder fastqBuilder = new FastqBuilder()
				.withDescription("description")
				.withSequence("sequence")
				.withQuality(null)
				.withVariant(FastqVariant.FASTQ_SOLEXA);

			fastqBuilder.build();
			Assert.fail("build null quality expected IllegalArgumentException");
		}
		catch (IllegalArgumentException e)
		{
			// expected
		}
	}

	@Test
	public void testBuildNullAppendQuality()
	{
		try
		{
			FastqBuilder fastqBuilder = new FastqBuilder()
				.withDescription("description")
				.withSequence("sequence")
				.appendQuality(null)
				.withVariant(FastqVariant.FASTQ_SOLEXA);

			fastqBuilder.build();
			Assert.fail("build null append quality expected IllegalArgumentException");
		}
		catch (IllegalArgumentException e)
		{
			// expected
		}
	}

	@Test
	public void testBuildNullVariant()
	{
		try
		{
			FastqBuilder fastqBuilder = new FastqBuilder()
				.withDescription("description")
				.withSequence("sequence")
				.withQuality("quality_")
				.withVariant(null);

			fastqBuilder.build();
			Assert.fail("build null variant expected IllegalArgumentException");
		}
		catch (IllegalArgumentException e)
		{
			// expected
		}
	}

	@Test
	public void testBuildMissingDescription()
	{
		try
		{
			FastqBuilder fastqBuilder = new FastqBuilder()
				.withSequence("sequence")
				.withQuality("quality_")
				.withVariant(FastqVariant.FASTQ_SOLEXA);

			fastqBuilder.build();
			Assert.fail("build missing description expected IllegalStateException");
		}
		catch (IllegalStateException e)
		{
			// expected
		}
	}

	@Test
	public void testBuildMissingSequence()
	{
		try
		{
			FastqBuilder fastqBuilder = new FastqBuilder()
				.withDescription("description")
				.withQuality("quality_")
				.withVariant(FastqVariant.FASTQ_SOLEXA);

			fastqBuilder.build();
			Assert.fail("build missing sequence expected IllegalStateException");
		}
		catch (IllegalStateException e)
		{
			// expected
		}
	}

	@Test
	public void testBuildMissingQuality()
	{
		try
		{
			FastqBuilder fastqBuilder = new FastqBuilder()
				.withDescription("description")
				.withSequence("sequence")
				.withVariant(FastqVariant.FASTQ_SOLEXA);

			fastqBuilder.build();
			Assert.fail("build missing quality expected IllegalStateException");
		}
		catch (IllegalStateException e)
		{
			// expected
		}
	}

	@Test
	public void testBuildDefaultVariant()
	{
		FastqBuilder fastqBuilder = new FastqBuilder()
			.withDescription("description")
			.withSequence("sequence")
			.withQuality("quality_");

		Fastq fastq = fastqBuilder.build();
		Assert.assertEquals("description", fastqBuilder.getDescription());
		Assert.assertTrue(fastqBuilder.sequenceAndQualityLengthsMatch());
		Assert.assertEquals("description", fastq.getDescription());
		Assert.assertEquals("sequence", fastq.getSequence());
		Assert.assertEquals("quality_", fastq.getQuality());
		Assert.assertEquals(FastqBuilder.DEFAULT_VARIANT, fastq.getVariant());
	}

	@Test
	public void testBuild()
	{
		FastqBuilder fastqBuilder = new FastqBuilder()
			.withDescription("description")
			.withSequence("sequence")
			.withQuality("quality_")
			.withVariant(FastqVariant.FASTQ_SOLEXA);
		Fastq fastq = fastqBuilder.build();
		Assert.assertEquals("description", fastqBuilder.getDescription());
		Assert.assertTrue(fastqBuilder.sequenceAndQualityLengthsMatch());
		Assert.assertEquals("description", fastq.getDescription());
		Assert.assertEquals("sequence", fastq.getSequence());
		Assert.assertEquals("quality_", fastq.getQuality());
		Assert.assertEquals(FastqVariant.FASTQ_SOLEXA, fastq.getVariant());
	}

	@Test
	public void testBuildAppendSequence()
	{
		FastqBuilder fastqBuilder = new FastqBuilder()
			.withDescription("description")
			.appendSequence("seq")
			.appendSequence("uence")
			.withQuality("quality_")
			.withVariant(FastqVariant.FASTQ_SOLEXA);
		Fastq fastq = fastqBuilder.build();
		Assert.assertEquals("description", fastqBuilder.getDescription());
		Assert.assertTrue(fastqBuilder.sequenceAndQualityLengthsMatch());
		Assert.assertEquals("description", fastq.getDescription());
		Assert.assertEquals("sequence", fastq.getSequence());
		Assert.assertEquals("quality_", fastq.getQuality());
		Assert.assertEquals(FastqVariant.FASTQ_SOLEXA, fastq.getVariant());
	}

	@Test
	public void testBuildAppendQuality()
	{
		FastqBuilder fastqBuilder = new FastqBuilder()
			.withDescription("description")
			.withSequence("sequence")
			.appendQuality("qual")
			.appendQuality("ity_")
			.withVariant(FastqVariant.FASTQ_SOLEXA);
		Fastq fastq = fastqBuilder.build();
		Assert.assertEquals("description", fastqBuilder.getDescription());
		Assert.assertTrue(fastqBuilder.sequenceAndQualityLengthsMatch());
		Assert.assertEquals("description", fastq.getDescription());
		Assert.assertEquals("sequence", fastq.getSequence());
		Assert.assertEquals("quality_", fastq.getQuality());
		Assert.assertEquals(FastqVariant.FASTQ_SOLEXA, fastq.getVariant());
	}

	@Test
	public void testBuildNonMatchingSequenceQualityScoreLengthsBothNull()
	{
		FastqBuilder fastqBuilder = new FastqBuilder()
			.withDescription("description")
			.withVariant(FastqVariant.FASTQ_SOLEXA);

		Assert.assertTrue(fastqBuilder.sequenceAndQualityLengthsMatch());
	}

	@Test
	public void testBuildNonMatchingSequenceQualityScoreLengthsSequenceNull()
	{
		FastqBuilder fastqBuilder = new FastqBuilder()
			.withDescription("description")
			.withQuality("0123")
			.withVariant(FastqVariant.FASTQ_SOLEXA);

		Assert.assertEquals(false, fastqBuilder.sequenceAndQualityLengthsMatch());
	}

	@Test
	public void testBuildNonMatchingSequenceQualityScoreLengthsQualityNull()
	{
		FastqBuilder fastqBuilder = new FastqBuilder()
			.withDescription("description")
			.withSequence("ACTG")
			.withVariant(FastqVariant.FASTQ_SOLEXA);

		Assert.assertEquals(false, fastqBuilder.sequenceAndQualityLengthsMatch());
	}

	@Test
	public void testBuildNonMatchingSequenceQualityScoreLengths0()
	{
		try
		{
			FastqBuilder fastqBuilder = new FastqBuilder()
				.withDescription("description")
				.withSequence("01234")
				.withQuality("0123")
				.withVariant(FastqVariant.FASTQ_SOLEXA);

			fastqBuilder.build();
			Assert.fail("build sequence length > quality length expected IllegalStateException");
		}
		catch (IllegalStateException e)
		{
			// expected
		}
	}

	@Test
	public void testBuildNonMatchingSequenceQualityScoreLengths1()
	{
		try
		{
			FastqBuilder fastqBuilder = new FastqBuilder()
				.withDescription("description")
				.withSequence("0123")
				.withQuality("01234")
				.withVariant(FastqVariant.FASTQ_SOLEXA);

			fastqBuilder.build();
			Assert.fail("build sequence length < quality length expected IllegalStateException");
		}
		catch (IllegalStateException e)
		{
			// expected
		}
	}

	@Test
	public void testBuildMultiple()
	{
		FastqBuilder fastqBuilder = new FastqBuilder()
			.withDescription("description")
			.withQuality("quality__")
			.withVariant(FastqVariant.FASTQ_SOLEXA);

		for (int i = 0; i < 10; i++)
		{
			Fastq fastq = fastqBuilder.withSequence("sequence" + i).build();
			Assert.assertEquals("description", fastqBuilder.getDescription());
			Assert.assertTrue(fastqBuilder.sequenceAndQualityLengthsMatch());
			Assert.assertEquals("description", fastq.getDescription());
			Assert.assertEquals("sequence" + i, fastq.getSequence());
			Assert.assertEquals("quality__", fastq.getQuality());
			Assert.assertEquals(FastqVariant.FASTQ_SOLEXA, fastq.getVariant());
		}
	}
}

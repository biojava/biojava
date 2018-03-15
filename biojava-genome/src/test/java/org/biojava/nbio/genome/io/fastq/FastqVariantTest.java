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

import static org.biojava.nbio.genome.io.fastq.FastqVariant.*;

import org.junit.Assert;
import org.junit.Test;

/**
 * Unit test for FastqVariant.
 */
public final class FastqVariantTest {

	@Test
	public void testDescription()
	{
		for (FastqVariant variant : values())
		{
			Assert.assertNotNull(variant.getDescription());
		}
	}

	@Test
	public void testIsSanger()
	{
		Assert.assertTrue(FASTQ_SANGER.isSanger());
		Assert.assertFalse(FASTQ_SOLEXA.isSanger());
		Assert.assertFalse(FASTQ_ILLUMINA.isSanger());
	}

	@Test
	public void testIsSolexa()
	{
		Assert.assertFalse(FASTQ_SANGER.isSolexa());
		Assert.assertTrue(FASTQ_SOLEXA.isSolexa());
		Assert.assertFalse(FASTQ_ILLUMINA.isSolexa());
	}

	@Test
	public void testIsIllumina()
	{
		Assert.assertFalse(FASTQ_SANGER.isIllumina());
		Assert.assertFalse(FASTQ_SOLEXA.isIllumina());
		Assert.assertTrue(FASTQ_ILLUMINA.isIllumina());
	}

	@Test
	public void testParseFastqVariant()
	{
		Assert.assertEquals(null, parseFastqVariant(null));
		Assert.assertEquals(null, parseFastqVariant(""));
		Assert.assertEquals(null, parseFastqVariant("not a valid FASTQ variant"));
		Assert.assertEquals(FASTQ_SANGER, parseFastqVariant("FASTQ_SANGER"));
		Assert.assertEquals(FASTQ_SANGER, parseFastqVariant("fastq-sanger"));
	}

	@Test
	public void testQualityLessThanMinimumQualityScore()
	{
		for (FastqVariant variant : values())
		{
			try
			{
				variant.quality(variant.minimumQualityScore() - 1);
				Assert.fail("expected IllegalArgumentException");
			}
			catch (IllegalArgumentException e)
			{
				// expected
			}
		}
	}

	@Test
	public void testQualityMoreThanMaximumQualityScore()
	{
		for (FastqVariant variant : values())
		{
			try
			{
				variant.quality(variant.maximumQualityScore() + 1);
				Assert.fail("expected IllegalArgumentException");
			}
			catch (IllegalArgumentException e)
			{
				// expected
			}
		}
	}

	@Test
	public void testQualityQualityScoreRoundTrip()
	{
		for (FastqVariant variant : values())
		{
			for (int i = variant.minimumQualityScore(); i < (variant.maximumQualityScore() + 1); i++)
			{
				Assert.assertEquals(i, variant.qualityScore(variant.quality(i)));
			}
		}
	}
}

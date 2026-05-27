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

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

/**
 * Unit test for FastqVariant.
 */
final class FastqVariantTest {

	@Test
    void testDescription()
	{
		for (FastqVariant variant : values())
		{
			Assertions.assertNotNull(variant.getDescription());
		}
	}

	@Test
    void testIsSanger()
	{
		Assertions.assertTrue(FASTQ_SANGER.isSanger());
		Assertions.assertFalse(FASTQ_SOLEXA.isSanger());
		Assertions.assertFalse(FASTQ_ILLUMINA.isSanger());
	}

	@Test
    void testIsSolexa()
	{
		Assertions.assertFalse(FASTQ_SANGER.isSolexa());
		Assertions.assertTrue(FASTQ_SOLEXA.isSolexa());
		Assertions.assertFalse(FASTQ_ILLUMINA.isSolexa());
	}

	@Test
    void testIsIllumina()
	{
		Assertions.assertFalse(FASTQ_SANGER.isIllumina());
		Assertions.assertFalse(FASTQ_SOLEXA.isIllumina());
		Assertions.assertTrue(FASTQ_ILLUMINA.isIllumina());
	}

	@Test
    void testParseFastqVariant()
	{
		Assertions.assertEquals(null, parseFastqVariant(null));
		Assertions.assertEquals(null, parseFastqVariant(""));
		Assertions.assertEquals(null, parseFastqVariant("not a valid FASTQ variant"));
		Assertions.assertEquals(FASTQ_SANGER, parseFastqVariant("FASTQ_SANGER"));
		Assertions.assertEquals(FASTQ_SANGER, parseFastqVariant("fastq-sanger"));
	}

	@Test
    void testQualityLessThanMinimumQualityScore()
	{
		for (FastqVariant variant : values())
		{
			try
			{
				variant.quality(variant.minimumQualityScore() - 1);
				Assertions.fail("expected IllegalArgumentException");
			}
			catch (IllegalArgumentException e)
			{
				// expected
			}
		}
	}

	@Test
    void testQualityMoreThanMaximumQualityScore()
	{
		for (FastqVariant variant : values())
		{
			try
			{
				variant.quality(variant.maximumQualityScore() + 1);
				Assertions.fail("expected IllegalArgumentException");
			}
			catch (IllegalArgumentException e)
			{
				// expected
			}
		}
	}

	@Test
    void testQualityQualityScoreRoundTrip()
	{
		for (FastqVariant variant : values())
		{
			for (int i = variant.minimumQualityScore(); i < (variant.maximumQualityScore() + 1); i++)
			{
				Assertions.assertEquals(i, variant.qualityScore(variant.quality(i)));
			}
		}
	}
}

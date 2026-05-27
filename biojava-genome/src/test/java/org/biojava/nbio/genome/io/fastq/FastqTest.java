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
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import org.junit.function.ThrowingRunnable;

/**
 * Unit test for Fastq.
 */
final class FastqTest {

	@Test
    void testConstructor()
	{
		Fastq fastq = new Fastq("description", "sequence", "quality_", FastqVariant.FASTQ_SANGER);
		Assertions.assertNotNull(fastq);

		try
		{
			new Fastq(null, "sequence", "quality_", FastqVariant.FASTQ_SANGER);
			Assertions.fail("ctr(null description) expected IllegalArgumentException");
		}
		catch (IllegalArgumentException e)
		{
			// expected
		}
		try
		{
			new Fastq("description", null, "quality_", FastqVariant.FASTQ_SANGER);
			Assertions.fail("ctr(null sequence) expected IllegalArgumentException");
		}
		catch (IllegalArgumentException e)
		{
			// expected
		}
		try
		{
			new Fastq("description", "sequence", null, FastqVariant.FASTQ_SANGER);
			Assertions.fail("ctr(null quality) expected IllegalArgumentException");
		}
		catch (IllegalArgumentException e)
		{
			// expected
		}
		try
		{
			new Fastq("description", "sequence", "quality_", null);
			Assertions.fail("ctr(null variant) expected IllegalArgumentException");
		}
		catch (IllegalArgumentException e)
		{
			// expected
		}
	}

	@Test
    void testDescription()
	{
		Fastq fastq = new Fastq("description", "sequence", "quality_", FastqVariant.FASTQ_SANGER);
		Assertions.assertTrue(fastq.getDescription() != null);
		Assertions.assertEquals("description", fastq.getDescription());
	}

	@Test
    void testSequence()
	{
		Fastq fastq = new Fastq("description", "sequence", "quality_", FastqVariant.FASTQ_SANGER);
		Assertions.assertTrue(fastq.getSequence() != null);
		Assertions.assertEquals("sequence", fastq.getSequence());
	}

	@Test
    void testQuality()
	{
		Fastq fastq = new Fastq("description", "sequence", "quality_", FastqVariant.FASTQ_SANGER);
		Assertions.assertTrue(fastq.getQuality() != null);
		Assertions.assertEquals("quality_", fastq.getQuality());
	}

	@Test
    void testVariant()
	{
		Fastq fastq = new Fastq("description", "sequence", "quality_", FastqVariant.FASTQ_SANGER);
		Assertions.assertTrue(fastq.getVariant() != null);
		Assertions.assertEquals(FastqVariant.FASTQ_SANGER, fastq.getVariant());
	}

	@Test
    void testBuilder()
	{
		Assertions.assertNotNull(Fastq.builder());
	}

	@Test
    void testBuilderNullFastq()
	{
		Assert.assertThrows(IllegalArgumentException.class, new ThrowingRunnable() {
			@Override
			public void run() {
				Fastq.builder(null);
			}
		});
	}

	@Test
    void testEquals()
	{
		Fastq fastq0 = new Fastq("description", "sequence", "quality_", FastqVariant.FASTQ_SANGER);
		Fastq fastq1 = new Fastq("description", "sequence", "quality_", FastqVariant.FASTQ_SANGER);

		Assertions.assertFalse(fastq0.equals(null));
		Assertions.assertFalse(fastq1.equals(null));
		Assertions.assertFalse(fastq0.equals(new Object()));
		Assertions.assertFalse(fastq1.equals(new Object()));
		Assertions.assertTrue(fastq0.equals(fastq0));
		Assertions.assertTrue(fastq1.equals(fastq1));
		Assertions.assertFalse(fastq0 == fastq1);
		Assertions.assertFalse(fastq0.equals(fastq1));
		Assertions.assertFalse(fastq1.equals(fastq0));
	}

	@Test
    void testHashCode()
	{
		Fastq fastq0 = new Fastq("description", "sequence", "quality_", FastqVariant.FASTQ_SANGER);
		Fastq fastq1 = new Fastq("description", "sequence", "quality_", FastqVariant.FASTQ_SANGER);

		Assertions.assertEquals(fastq0.hashCode(), fastq0.hashCode());
		Assertions.assertEquals(fastq1.hashCode(), fastq1.hashCode());
		if (fastq0.equals(fastq1))
		{
			Assertions.assertEquals(fastq0.hashCode(), fastq1.hashCode());
			Assertions.assertEquals(fastq1.hashCode(), fastq0.hashCode());
		}
		if (fastq1.equals(fastq0))
		{
			Assertions.assertEquals(fastq0.hashCode(), fastq1.hashCode());
			Assertions.assertEquals(fastq1.hashCode(), fastq0.hashCode());
		}
	}
}

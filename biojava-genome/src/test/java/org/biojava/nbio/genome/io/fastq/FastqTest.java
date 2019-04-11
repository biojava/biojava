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
 * Unit test for Fastq.
 */
public final class FastqTest {

	@Test
	public void testConstructor()
	{
		Fastq fastq = new Fastq("description", "sequence", "quality_", FastqVariant.FASTQ_SANGER);
		Assert.assertNotNull(fastq);

		try
		{
			new Fastq(null, "sequence", "quality_", FastqVariant.FASTQ_SANGER);
			Assert.fail("ctr(null description) expected IllegalArgumentException");
		}
		catch (IllegalArgumentException e)
		{
			// expected
		}
		try
		{
			new Fastq("description", null, "quality_", FastqVariant.FASTQ_SANGER);
			Assert.fail("ctr(null sequence) expected IllegalArgumentException");
		}
		catch (IllegalArgumentException e)
		{
			// expected
		}
		try
		{
			new Fastq("description", "sequence", null, FastqVariant.FASTQ_SANGER);
			Assert.fail("ctr(null quality) expected IllegalArgumentException");
		}
		catch (IllegalArgumentException e)
		{
			// expected
		}
		try
		{
			new Fastq("description", "sequence", "quality_", null);
			Assert.fail("ctr(null variant) expected IllegalArgumentException");
		}
		catch (IllegalArgumentException e)
		{
			// expected
		}
	}

	@Test
	public void testDescription()
	{
		Fastq fastq = new Fastq("description", "sequence", "quality_", FastqVariant.FASTQ_SANGER);
		Assert.assertTrue(fastq.getDescription() != null);
		Assert.assertEquals("description", fastq.getDescription());
	}

	@Test
	public void testSequence()
	{
		Fastq fastq = new Fastq("description", "sequence", "quality_", FastqVariant.FASTQ_SANGER);
		Assert.assertTrue(fastq.getSequence() != null);
		Assert.assertEquals("sequence", fastq.getSequence());
	}

	@Test
	public void testQuality()
	{
		Fastq fastq = new Fastq("description", "sequence", "quality_", FastqVariant.FASTQ_SANGER);
		Assert.assertTrue(fastq.getQuality() != null);
		Assert.assertEquals("quality_", fastq.getQuality());
	}

	@Test
	public void testVariant()
	{
		Fastq fastq = new Fastq("description", "sequence", "quality_", FastqVariant.FASTQ_SANGER);
		Assert.assertTrue(fastq.getVariant() != null);
		Assert.assertEquals(FastqVariant.FASTQ_SANGER, fastq.getVariant());
	}

	@Test
	public void testBuilder()
	{
		Assert.assertNotNull(Fastq.builder());
	}

	@Test
	public void testEquals()
	{
		Fastq fastq0 = new Fastq("description", "sequence", "quality_", FastqVariant.FASTQ_SANGER);
		Fastq fastq1 = new Fastq("description", "sequence", "quality_", FastqVariant.FASTQ_SANGER);

		Assert.assertFalse(fastq0.equals(null));
		Assert.assertFalse(fastq1.equals(null));
		Assert.assertFalse(fastq0.equals(new Object()));
		Assert.assertFalse(fastq1.equals(new Object()));
		Assert.assertTrue(fastq0.equals(fastq0));
		Assert.assertTrue(fastq1.equals(fastq1));
		Assert.assertFalse(fastq0 == fastq1);
		Assert.assertFalse(fastq0.equals(fastq1));
		Assert.assertFalse(fastq1.equals(fastq0));
	}

	@Test
	public void testHashCode()
	{
		Fastq fastq0 = new Fastq("description", "sequence", "quality_", FastqVariant.FASTQ_SANGER);
		Fastq fastq1 = new Fastq("description", "sequence", "quality_", FastqVariant.FASTQ_SANGER);

		Assert.assertEquals(fastq0.hashCode(), fastq0.hashCode());
		Assert.assertEquals(fastq1.hashCode(), fastq1.hashCode());
		if (fastq0.equals(fastq1))
		{
			Assert.assertEquals(fastq0.hashCode(), fastq1.hashCode());
			Assert.assertEquals(fastq1.hashCode(), fastq0.hashCode());
		}
		if (fastq1.equals(fastq0))
		{
			Assert.assertEquals(fastq0.hashCode(), fastq1.hashCode());
			Assert.assertEquals(fastq1.hashCode(), fastq0.hashCode());
		}
	}
}

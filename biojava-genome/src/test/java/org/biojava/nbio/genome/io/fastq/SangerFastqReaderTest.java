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

import org.junit.Test;

import java.io.IOException;
import java.io.InputStream;
import java.net.URL;

import static org.junit.Assert.*;

/**
 * Unit test for SangerFastqReader.
 */
public final class SangerFastqReaderTest
	extends AbstractFastqReaderTest
{

	@Override
	public Fastq createFastq()
	{
		return new FastqBuilder()
			.withDescription("description")
			.withSequence("sequence")
			.withQuality("quality_")
			.withVariant(FastqVariant.FASTQ_SANGER)
			.build();
	}

	@Override
	public FastqReader createFastqReader()
	{
		return new SangerFastqReader();
	}

	@Override
	public FastqWriter createFastqWriter()
	{
		return new SangerFastqWriter();
	}

	public void testValidateDescription() throws Exception
	{
		SangerFastqReader reader = new SangerFastqReader();
		URL invalidDescription = getClass().getResource("sanger-invalid-description.fastq");
		try
		{
			reader.read(invalidDescription);
			fail("read(invalidDescription) expected IOException");
		}
		catch (IOException e)
		{
			assertTrue(e.getMessage().contains("description must begin with a '@' character"));
		}
	}

	@Test
	public void testValidateRepeatDescription() throws Exception
	{
		SangerFastqReader reader = new SangerFastqReader();
		URL invalidRepeatDescription = getClass().getResource("sanger-invalid-repeat-description.fastq");
		try
		{
			reader.read(invalidRepeatDescription);
			fail("read(invalidRepeatDescription) expected IOException");
		}
		catch (IOException e)
		{
			assertTrue(e.getMessage().contains("repeat description must match description"));
		}
	}

	@Test
	public void testWrappingOriginal() throws Exception
	{
		FastqReader reader = createFastqReader();
		InputStream inputStream = getClass().getResourceAsStream("wrapping_original_sanger.fastq");
		Iterable<Fastq> iterable = reader.read(inputStream);
		assertNotNull(iterable);
		int count = 0;
		for (Fastq f : iterable)
		{
			assertNotNull(f);
			count++;
		}
		assertEquals(3, count);
		inputStream.close();
	}

	@Test
	public void testWrappingAsSanger() throws Exception
	{
		FastqReader reader = createFastqReader();
		InputStream inputStream = getClass().getResourceAsStream("wrapping_as_sanger.fastq");
		Iterable<Fastq> iterable = reader.read(inputStream);
		assertNotNull(iterable);
		int count = 0;
		for (Fastq f : iterable)
		{
			assertNotNull(f);
			count++;
		}
		assertEquals(3, count);
		inputStream.close();
	}

	@Test
	public void testFullRangeOriginal() throws Exception
	{
		FastqReader reader = createFastqReader();
		InputStream inputStream = getClass().getResourceAsStream("sanger_full_range_original_sanger.fastq");
		Iterable<Fastq> iterable = reader.read(inputStream);
		assertNotNull(iterable);
		int count = 0;
		for (Fastq f : iterable)
		{
			assertNotNull(f);
			count++;
		}
		assertEquals(2, count);
		inputStream.close();
	}

	@Test
	public void testFullRangeAsSanger() throws Exception
	{
		FastqReader reader = createFastqReader();
		InputStream inputStream = getClass().getResourceAsStream("sanger_full_range_as_sanger.fastq");
		Iterable<Fastq> iterable = reader.read(inputStream);
		assertNotNull(iterable);
		int count = 0;
		for (Fastq f : iterable)
		{
			assertNotNull(f);
			count++;
		}
		assertEquals(2, count);
		inputStream.close();
	}

	@Test
	public void testMiscDnaOriginal() throws Exception
	{
		FastqReader reader = createFastqReader();
		InputStream inputStream = getClass().getResourceAsStream("misc_dna_original_sanger.fastq");
		Iterable<Fastq> iterable = reader.read(inputStream);
		assertNotNull(iterable);
		int count = 0;
		for (Fastq f : iterable)
		{
			assertNotNull(f);
			count++;
		}
		assertEquals(4, count);
		inputStream.close();
	}

	@Test
	public void testMiscDnaAsSanger() throws Exception
	{
		FastqReader reader = createFastqReader();
		InputStream inputStream = getClass().getResourceAsStream("misc_dna_as_sanger.fastq");
		Iterable<Fastq> iterable = reader.read(inputStream);
		assertNotNull(iterable);
		int count = 0;
		for (Fastq f : iterable)
		{
			assertNotNull(f);
			count++;
		}
		assertEquals(4, count);
		inputStream.close();
	}

	@Test
	public void testMiscRnaOriginal() throws Exception
	{
		FastqReader reader = createFastqReader();
		InputStream inputStream = getClass().getResourceAsStream("misc_rna_original_sanger.fastq");
		Iterable<Fastq> iterable = reader.read(inputStream);
		assertNotNull(iterable);
		int count = 0;
		for (Fastq f : iterable)
		{
			assertNotNull(f);
			count++;
		}
		assertEquals(4, count);
		inputStream.close();
	}

	@Test
	public void testMiscRnaAsSanger() throws Exception
	{
		FastqReader reader = createFastqReader();
		InputStream inputStream = getClass().getResourceAsStream("misc_rna_as_sanger.fastq");
		Iterable<Fastq> iterable = reader.read(inputStream);
		assertNotNull(iterable);
		int count = 0;
		for (Fastq f : iterable)
		{
			assertNotNull(f);
			count++;
		}
		assertEquals(4, count);
		inputStream.close();
	}

	@Test
	public void testLongReadsOriginal() throws Exception
	{
		FastqReader reader = createFastqReader();
		InputStream inputStream = getClass().getResourceAsStream("longreads_original_sanger.fastq");
		Iterable<Fastq> iterable = reader.read(inputStream);
		assertNotNull(iterable);
		int count = 0;
		for (Fastq f : iterable)
		{
			assertNotNull(f);
			count++;
		}
		assertEquals(10, count);
		inputStream.close();
	}

	@Test
	public void testLongReadsAsSanger() throws Exception
	{
		FastqReader reader = createFastqReader();
		InputStream inputStream = getClass().getResourceAsStream("longreads_as_sanger.fastq");
		Iterable<Fastq> iterable = reader.read(inputStream);
		assertNotNull(iterable);
		int count = 0;
		for (Fastq f : iterable)
		{
			assertNotNull(f);
			count++;
		}
		assertEquals(10, count);
		inputStream.close();
	}
}

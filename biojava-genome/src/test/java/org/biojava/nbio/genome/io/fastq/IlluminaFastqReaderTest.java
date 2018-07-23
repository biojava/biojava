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
import static org.junit.Assert.*;

import java.io.IOException;
import java.io.InputStream;
import java.net.URL;


/**
 * Unit test for IlluminaFastqReader.
 */
public final class IlluminaFastqReaderTest
	extends AbstractFastqReaderTest
{

	@Override
	public Fastq createFastq()
	{
		return new FastqBuilder()
			.withDescription("description")
			.withSequence("sequence")
			.withQuality("quality_")
			.withVariant(FastqVariant.FASTQ_ILLUMINA)
			.build();
	}

	@Override
	public FastqReader createFastqReader()
	{
		return new IlluminaFastqReader();
	}

	@Override
	public FastqWriter createFastqWriter()
	{
		return new IlluminaFastqWriter();
	}

	@Test
	public void testValidateDescription() throws Exception
	{
		IlluminaFastqReader reader = new IlluminaFastqReader();
		URL invalidDescription = getClass().getResource("illumina-invalid-description.fastq");
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
		IlluminaFastqReader reader = new IlluminaFastqReader();
		URL invalidRepeatDescription = getClass().getResource("illumina-invalid-repeat-description.fastq");
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
	public void testWrappingAsIllumina() throws Exception
	{
		FastqReader reader = createFastqReader();
		InputStream inputStream = getClass().getResourceAsStream("wrapping_as_illumina.fastq");
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
	public void testFullRangeAsIllumina() throws Exception
	{
		FastqReader reader = createFastqReader();
		InputStream inputStream = getClass().getResourceAsStream("illumina_full_range_as_illumina.fastq");
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
	public void testMiscDnaAsIllumina() throws Exception
	{
		FastqReader reader = createFastqReader();
		InputStream inputStream = getClass().getResourceAsStream("misc_dna_as_illumina.fastq");
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
	public void testMiscRnaAsIllumina() throws Exception
	{
		FastqReader reader = createFastqReader();
		InputStream inputStream = getClass().getResourceAsStream("misc_rna_as_illumina.fastq");
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
	public void testLongReadsAsIllumina() throws Exception
	{
		FastqReader reader = createFastqReader();
		InputStream inputStream = getClass().getResourceAsStream("longreads_as_illumina.fastq");
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

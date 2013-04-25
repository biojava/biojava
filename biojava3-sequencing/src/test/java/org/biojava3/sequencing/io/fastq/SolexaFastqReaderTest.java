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
package org.biojava3.sequencing.io.fastq;

import java.io.InputStream;
import java.io.IOException;

import java.net.URL;

/**
 * Unit test for SolexaFastqReader.
 */
public final class SolexaFastqReaderTest
    extends AbstractFastqReaderTest
{

    @Override
    public Fastq createFastq()
    {
        return new FastqBuilder()
            .withDescription("description")
            .withSequence("sequence")
            .withQuality("quality_")
            .withVariant(FastqVariant.FASTQ_SOLEXA)
            .build();
    }

    @Override
    public FastqReader createFastqReader()
    {
        return new SolexaFastqReader();
    }

    @Override
    public FastqWriter createFastqWriter()
    {
        return new SolexaFastqWriter();
    }

    public void testValidateDescription() throws Exception
    {
        SolexaFastqReader reader = new SolexaFastqReader();
        URL invalidDescription = getClass().getResource("solexa-invalid-description.fastq");
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

    public void testValidateRepeatDescription() throws Exception
    {
        SolexaFastqReader reader = new SolexaFastqReader();
        URL invalidRepeatDescription = getClass().getResource("solexa-invalid-repeat-description.fastq");
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

    public void testWrappingAsSolexa() throws Exception
    {
        FastqReader reader = createFastqReader();
        InputStream inputStream = getClass().getResourceAsStream("wrapping_as_solexa.fastq");
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

    public void testFullRangeAsSolexa() throws Exception
    {
        FastqReader reader = createFastqReader();
        InputStream inputStream = getClass().getResourceAsStream("solexa_full_range_as_solexa.fastq");
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

    public void testMiscDnaAsSolexa() throws Exception
    {
        FastqReader reader = createFastqReader();
        InputStream inputStream = getClass().getResourceAsStream("misc_dna_as_solexa.fastq");
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

    public void testMiscRnaAsSolexa() throws Exception
    {
        FastqReader reader = createFastqReader();
        InputStream inputStream = getClass().getResourceAsStream("misc_rna_as_solexa.fastq");
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

    public void testLongReadsAsSolexa() throws Exception
    {
        FastqReader reader = createFastqReader();
        InputStream inputStream = getClass().getResourceAsStream("longreads_as_solexa.fastq");
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
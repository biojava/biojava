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
package org.biojava.bio.program.fastq;

import java.io.InputStream;
import java.io.IOException;
import java.io.StringReader;

import java.net.URL;

import com.google.common.io.InputSupplier;

/**
 * Unit test for SangerFastqReader.
 */
public final class SangerFastqReaderTest
    extends AbstractFastqReaderTest
{

    /** {@inheritDoc} */
    public Fastq createFastq()
    {
        return new FastqBuilder()
            .withDescription("description")
            .withSequence("sequence")
            .withQuality("quality_")
            .withVariant(FastqVariant.FASTQ_SANGER)
            .build();
    }

    /** {@inheritDoc} */
    public FastqReader createFastqReader()
    {
        return new SangerFastqReader();
    }

    /** {@inheritDoc} */
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

    public void testParse() throws Exception
    {
        SangerFastqReader sangerFastqReader = (SangerFastqReader) createFastqReader();
        final String input = "";
        sangerFastqReader.parse(new InputSupplier<StringReader>()
                                {
                                    /** {@inheritDoc} */
                                    public StringReader getInput() throws IOException {
                                        return new StringReader(input);
                                    }
                                },
                                new ParseListener() {
                                    /** {@inheritDoc} */
                                    public void description(final String description) throws IOException {
                                        // empty
                                    }
 
                                    /** {@inheritDoc} */
                                    public void sequence(final String sequence) throws IOException {
                                        // empty
                                    }
 
                                    /** {@inheritDoc} */
                                    public void appendSequence(final String sequence) throws IOException {
                                        // empty
                                    }
 
                                    /** {@inheritDoc} */
                                    public void repeatDescription(final String repeatDescription) throws IOException {
                                        // empty
                                    }
 
                                    /** {@inheritDoc} */
                                    public void quality(final String quality) throws IOException {
                                        // empty
                                    }
 
                                    /** {@inheritDoc} */
                                    public void appendQuality(final String quality) throws IOException {
                                        // empty
                                    }
 
                                    /** {@inheritDoc} */
                                    public void complete() throws IOException {
                                        // empty
                                    }
                                });
    }
 
    public void testParseNullInputSupplier() throws Exception
    {
        SangerFastqReader sangerFastqReader = (SangerFastqReader) createFastqReader();
        try
        {
            sangerFastqReader.parse((InputSupplier<StringReader>) null,
                                    new ParseListener() {
                                        /** {@inheritDoc} */
                                        public void description(final String description) throws IOException {
                                            // empty
                                        }
 
                                        /** {@inheritDoc} */
                                        public void sequence(final String sequence) throws IOException {
                                            // empty
                                        }
 
                                        /** {@inheritDoc} */
                                        public void appendSequence(final String sequence) throws IOException {
                                            // empty
                                        }
 
                                        /** {@inheritDoc} */
                                        public void repeatDescription(final String repeatDescription) throws IOException {
                                            // empty
                                        }
 
                                        /** {@inheritDoc} */
                                        public void quality(final String quality) throws IOException {
                                            // empty
                                        }
 
                                        /** {@inheritDoc} */
                                        public void appendQuality(final String quality) throws IOException {
                                            // empty
                                        }

                                        /** {@inheritDoc} */
                                        public void complete() throws IOException {
                                            // empty
                                        }
                                    });
            fail("parse(null, ) expected IllegalArgumentException");
        }
        catch (IllegalArgumentException e)
        {
            // expected
        }
    }
 
    public void testParseNullParseListener() throws Exception
    {
        SangerFastqReader sangerFastqReader = (SangerFastqReader) createFastqReader();
        final String input = "";
        try
        {
            sangerFastqReader.parse(new InputSupplier<StringReader>()
                                    {
                                        /** {@inheritDoc} */
                                        public StringReader getInput() throws IOException {
                                            return new StringReader(input);
                                        }
                                    }, null);
            fail("parse(, null) expected IllegalArgumentException");
        }
        catch (IllegalArgumentException e)
        {
            // expected
        }
    }
}
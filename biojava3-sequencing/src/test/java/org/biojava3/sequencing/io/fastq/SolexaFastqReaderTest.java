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
import java.io.StringReader;

import java.net.URL;

import com.google.common.io.InputSupplier;

/**
 * Unit test for SolexaFastqReader.
 */
public final class SolexaFastqReaderTest
    extends AbstractFastqReaderTest
{

    /** {@inheritDoc} */
    public Fastq createFastq()
    {
        return new FastqBuilder()
            .withDescription("description")
            .withSequence("sequence")
            .withQuality("quality_")
            .withVariant(FastqVariant.FASTQ_SOLEXA)
            .build();
    }

    /** {@inheritDoc} */
    public FastqReader createFastqReader()
    {
        return new SolexaFastqReader();
    }

    /** {@inheritDoc} */
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

    public void testParse() throws Exception
    {
        SolexaFastqReader solexaFastqReader = (SolexaFastqReader) createFastqReader();
        final String input = "";
        solexaFastqReader.parse(new InputSupplier<StringReader>()
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
        SolexaFastqReader solexaFastqReader = (SolexaFastqReader) createFastqReader();
        try
        {
            solexaFastqReader.parse((InputSupplier<StringReader>) null,
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
        SolexaFastqReader solexaFastqReader = (SolexaFastqReader) createFastqReader();
        final String input = "";
        try
        {
            solexaFastqReader.parse(new InputSupplier<StringReader>()
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
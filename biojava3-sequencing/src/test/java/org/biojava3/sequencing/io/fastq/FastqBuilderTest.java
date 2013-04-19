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

import junit.framework.TestCase;

/**
 * Unit test for FastqBuilder.
 */
public final class FastqBuilderTest
    extends TestCase
{

    public void testConstructor()
    {
        FastqBuilder fastqBuilder = new FastqBuilder();
        assertNotNull(fastqBuilder);
    }

    public void testBuildDefault()
    {
        try
        {
            FastqBuilder fastqBuilder = new FastqBuilder();
            fastqBuilder.build();
            fail("build default expected IllegalStateException");
        }
        catch (IllegalStateException e)
        {
            // expected
        }
    }

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
            fail("build null description expected IllegalArgumentException");
        }
        catch (IllegalArgumentException e)
        {
            // expected
        }
    }

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
            fail("build null sequence expected IllegalArgumentException");
        }
        catch (IllegalArgumentException e)
        {
            // expected
        }
    }

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
            fail("build null append sequence expected IllegalArgumentException");
        }
        catch (IllegalArgumentException e)
        {
            // expected
        }
    }

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
            fail("build null quality expected IllegalArgumentException");
        }
        catch (IllegalArgumentException e)
        {
            // expected
        }
    }

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
            fail("build null append quality expected IllegalArgumentException");
        }
        catch (IllegalArgumentException e)
        {
            // expected
        }
    }

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
            fail("build null variant expected IllegalArgumentException");
        }
        catch (IllegalArgumentException e)
        {
            // expected
        }
    }

    public void testBuildMissingDescription()
    {
        try
        {
            FastqBuilder fastqBuilder = new FastqBuilder()
                .withSequence("sequence")
                .withQuality("quality_")
                .withVariant(FastqVariant.FASTQ_SOLEXA);

            fastqBuilder.build();
            fail("build missing description expected IllegalStateException");
        }
        catch (IllegalStateException e)
        {
            // expected
        }
    }

    public void testBuildMissingSequence()
    {
        try
        {
            FastqBuilder fastqBuilder = new FastqBuilder()
                .withDescription("description")
                .withQuality("quality_")
                .withVariant(FastqVariant.FASTQ_SOLEXA);

            fastqBuilder.build();
            fail("build missing sequence expected IllegalStateException");
        }
        catch (IllegalStateException e)
        {
            // expected
        }
    }

    public void testBuildMissingQuality()
    {
        try
        {
            FastqBuilder fastqBuilder = new FastqBuilder()
                .withDescription("description")
                .withSequence("sequence")
                .withVariant(FastqVariant.FASTQ_SOLEXA);

            fastqBuilder.build();
            fail("build missing quality expected IllegalStateException");
        }
        catch (IllegalStateException e)
        {
            // expected
        }
    }

    public void testBuildDefaultVariant()
    {
        FastqBuilder fastqBuilder = new FastqBuilder()
            .withDescription("description")
            .withSequence("sequence")
            .withQuality("quality_");

        Fastq fastq = fastqBuilder.build();
        assertEquals("description", fastqBuilder.getDescription());
        assertTrue(fastqBuilder.sequenceAndQualityLengthsMatch());
        assertEquals("description", fastq.getDescription());
        assertEquals("sequence", fastq.getSequence());
        assertEquals("quality_", fastq.getQuality());
        assertEquals(FastqBuilder.DEFAULT_VARIANT, fastq.getVariant());
    }

    public void testBuild()
    {
        FastqBuilder fastqBuilder = new FastqBuilder()
            .withDescription("description")
            .withSequence("sequence")
            .withQuality("quality_")
            .withVariant(FastqVariant.FASTQ_SOLEXA);
        Fastq fastq = fastqBuilder.build();
        assertEquals("description", fastqBuilder.getDescription());
        assertTrue(fastqBuilder.sequenceAndQualityLengthsMatch());
        assertEquals("description", fastq.getDescription());
        assertEquals("sequence", fastq.getSequence());
        assertEquals("quality_", fastq.getQuality());
        assertEquals(FastqVariant.FASTQ_SOLEXA, fastq.getVariant());
    }

    public void testBuildAppendSequence()
    {
        FastqBuilder fastqBuilder = new FastqBuilder()
            .withDescription("description")
            .appendSequence("seq")
            .appendSequence("uence")
            .withQuality("quality_")
            .withVariant(FastqVariant.FASTQ_SOLEXA);
        Fastq fastq = fastqBuilder.build();
        assertEquals("description", fastqBuilder.getDescription());
        assertTrue(fastqBuilder.sequenceAndQualityLengthsMatch());
        assertEquals("description", fastq.getDescription());
        assertEquals("sequence", fastq.getSequence());
        assertEquals("quality_", fastq.getQuality());
        assertEquals(FastqVariant.FASTQ_SOLEXA, fastq.getVariant());
    }

    public void testBuildAppendQuality()
    {
        FastqBuilder fastqBuilder = new FastqBuilder()
            .withDescription("description")
            .withSequence("sequence")
            .appendQuality("qual")
            .appendQuality("ity_")
            .withVariant(FastqVariant.FASTQ_SOLEXA);
        Fastq fastq = fastqBuilder.build();
        assertEquals("description", fastqBuilder.getDescription());
        assertTrue(fastqBuilder.sequenceAndQualityLengthsMatch());
        assertEquals("description", fastq.getDescription());
        assertEquals("sequence", fastq.getSequence());
        assertEquals("quality_", fastq.getQuality());
        assertEquals(FastqVariant.FASTQ_SOLEXA, fastq.getVariant());
    }

    public void testBuildNonMatchingSequenceQualityScoreLengthsBothNull()
    {
        FastqBuilder fastqBuilder = new FastqBuilder()
            .withDescription("description")
            .withVariant(FastqVariant.FASTQ_SOLEXA);

        assertTrue(fastqBuilder.sequenceAndQualityLengthsMatch());
    }

    public void testBuildNonMatchingSequenceQualityScoreLengthsSequenceNull()
    {
        FastqBuilder fastqBuilder = new FastqBuilder()
            .withDescription("description")
            .withQuality("0123")
            .withVariant(FastqVariant.FASTQ_SOLEXA);

        assertEquals(false, fastqBuilder.sequenceAndQualityLengthsMatch());
    }

    public void testBuildNonMatchingSequenceQualityScoreLengthsQualityNull()
    {
        FastqBuilder fastqBuilder = new FastqBuilder()
            .withDescription("description")
            .withSequence("ACTG")
            .withVariant(FastqVariant.FASTQ_SOLEXA);

        assertEquals(false, fastqBuilder.sequenceAndQualityLengthsMatch());
    }

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
            fail("build sequence length > quality length expected IllegalStateException");
        }
        catch (IllegalStateException e)
        {
            // expected
        }
    }

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
            fail("build sequence length < quality length expected IllegalStateException");
        }
        catch (IllegalStateException e)
        {
            // expected
        }
    }

    public void testBuildMultiple()
    {
        FastqBuilder fastqBuilder = new FastqBuilder()
            .withDescription("description")
            .withQuality("quality__")
            .withVariant(FastqVariant.FASTQ_SOLEXA);

        for (int i = 0; i < 10; i++)
        {
            Fastq fastq = fastqBuilder.withSequence("sequence" + i).build();
            assertEquals("description", fastqBuilder.getDescription());
            assertTrue(fastqBuilder.sequenceAndQualityLengthsMatch());
            assertEquals("description", fastq.getDescription());
            assertEquals("sequence" + i, fastq.getSequence());
            assertEquals("quality__", fastq.getQuality());
            assertEquals(FastqVariant.FASTQ_SOLEXA, fastq.getVariant());
        }
    }
}
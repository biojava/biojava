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

import static org.biojava3.sequencing.io.fastq.FastqVariant.*;

/**
 * Unit test for FastqVariant.
 */
public final class FastqVariantTest
    extends TestCase
{

    public void testDescription()
    {
        for (FastqVariant variant : values())
        {
            assertNotNull(variant.getDescription());
        }
    }

    public void testIsSanger()
    {
        assertTrue(FASTQ_SANGER.isSanger());
        assertFalse(FASTQ_SOLEXA.isSanger());
        assertFalse(FASTQ_ILLUMINA.isSanger());
    }

    public void testIsSolexa()
    {
        assertFalse(FASTQ_SANGER.isSolexa());
        assertTrue(FASTQ_SOLEXA.isSolexa());
        assertFalse(FASTQ_ILLUMINA.isSolexa());
    }

    public void testIsIllumina()
    {
        assertFalse(FASTQ_SANGER.isIllumina());
        assertFalse(FASTQ_SOLEXA.isIllumina());
        assertTrue(FASTQ_ILLUMINA.isIllumina());
    }

    public void testParseFastqVariant()
    {
        assertEquals(null, parseFastqVariant(null));
        assertEquals(null, parseFastqVariant(""));
        assertEquals(null, parseFastqVariant("not a valid FASTQ variant"));
        assertEquals(FASTQ_SANGER, parseFastqVariant("FASTQ_SANGER"));
        assertEquals(FASTQ_SANGER, parseFastqVariant("fastq-sanger"));
    }

    public void testQualityLessThanMinimumQualityScore()
    {
        for (FastqVariant variant : values())
        {
            try
            {
                variant.quality(variant.minimumQualityScore() - 1);
                fail("expected IllegalArgumentException");
            }
            catch (IllegalArgumentException e)
            {
                // expected
            }
        }
    }

    public void testQualityMoreThanMaximumQualityScore()
    {
        for (FastqVariant variant : values())
        {
            try
            {
                variant.quality(variant.maximumQualityScore() + 1);
                fail("expected IllegalArgumentException");
            }
            catch (IllegalArgumentException e)
            {
                // expected
            }
        }
    }

    public void testQualityQualityScoreRoundTrip()
    {
        for (FastqVariant variant : values())
        {
            for (int i = variant.minimumQualityScore(); i < (variant.maximumQualityScore() + 1); i++)
            {
                assertEquals(i, variant.qualityScore(variant.quality(i)));
            }
        }
    }
}
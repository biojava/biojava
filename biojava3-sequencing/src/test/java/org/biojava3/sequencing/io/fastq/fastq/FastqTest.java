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

import java.lang.reflect.Field;
import java.lang.reflect.Modifier;

import junit.framework.TestCase;

/**
 * Unit test for Fastq.
 */
public final class FastqTest
    extends TestCase
{

    public void testImmutable()
    {
        Class<Fastq> cls = Fastq.class;
        assertTrue(Modifier.isPublic(cls.getModifiers()));
        assertTrue(Modifier.isFinal(cls.getModifiers()));
        Field[] fields = cls.getDeclaredFields();
        for (Field field : fields)
        {
            assertTrue(Modifier.isPrivate(field.getModifiers()));
            assertTrue(Modifier.isFinal(field.getModifiers()) ||
                    (Modifier.isVolatile(field.getModifiers()) && Modifier.isTransient(field.getModifiers())));
        }
    }

    public void testConstructor()
    {
        Fastq fastq = new Fastq("description", "sequence", "quality_", FastqVariant.FASTQ_SANGER);
        assertNotNull(fastq);

        try
        {
            new Fastq(null, "sequence", "quality_", FastqVariant.FASTQ_SANGER);
            fail("ctr(null description) expected IllegalArgumentException");
        }
        catch (IllegalArgumentException e)
        {
            // expected
        }
        try
        {
            new Fastq("description", null, "quality_", FastqVariant.FASTQ_SANGER);
            fail("ctr(null sequence) expected IllegalArgumentException");
        }
        catch (IllegalArgumentException e)
        {
            // expected
        }
        try
        {
            new Fastq("description", "sequence", null, FastqVariant.FASTQ_SANGER);
            fail("ctr(null quality) expected IllegalArgumentException");
        }
        catch (IllegalArgumentException e)
        {
            // expected
        }
        try
        {
            new Fastq("description", "sequence", "quality_", null);
            fail("ctr(null variant) expected IllegalArgumentException");
        }
        catch (IllegalArgumentException e)
        {
            // expected
        }
    }

    public void testDescription()
    {
        Fastq fastq = new Fastq("description", "sequence", "quality_", FastqVariant.FASTQ_SANGER);
        assertTrue(fastq.getDescription() != null);
        assertEquals("description", fastq.getDescription());
    }

    public void testSequence()
    {
        Fastq fastq = new Fastq("description", "sequence", "quality_", FastqVariant.FASTQ_SANGER);
        assertTrue(fastq.getSequence() != null);
        assertEquals("sequence", fastq.getSequence());
    }

    public void testQuality()
    {
        Fastq fastq = new Fastq("description", "sequence", "quality_", FastqVariant.FASTQ_SANGER);
        assertTrue(fastq.getQuality() != null);
        assertEquals("quality_", fastq.getQuality());
    }

    public void testVariant()
    {
        Fastq fastq = new Fastq("description", "sequence", "quality_", FastqVariant.FASTQ_SANGER);
        assertTrue(fastq.getVariant() != null);
        assertEquals(FastqVariant.FASTQ_SANGER, fastq.getVariant());
    }

    public void testEquals()
    {
        Fastq fastq0 = new Fastq("description", "sequence", "quality_", FastqVariant.FASTQ_SANGER);
        Fastq fastq1 = new Fastq("description", "sequence", "quality_", FastqVariant.FASTQ_SANGER);

        assertFalse(fastq0.equals(null));
        assertFalse(fastq1.equals(null));
        assertFalse(fastq0.equals(new Object()));
        assertFalse(fastq1.equals(new Object()));
        assertTrue(fastq0.equals(fastq0));
        assertTrue(fastq1.equals(fastq1));
        assertFalse(fastq0 == fastq1);
        assertFalse(fastq0.equals(fastq1));
        assertFalse(fastq1.equals(fastq0));
    }

    public void testHashCode()
    {
        Fastq fastq0 = new Fastq("description", "sequence", "quality_", FastqVariant.FASTQ_SANGER);
        Fastq fastq1 = new Fastq("description", "sequence", "quality_", FastqVariant.FASTQ_SANGER);

        assertEquals(fastq0.hashCode(), fastq0.hashCode());
        assertEquals(fastq1.hashCode(), fastq1.hashCode());
        if (fastq0.equals(fastq1))
        {
            assertEquals(fastq0.hashCode(), fastq1.hashCode());
            assertEquals(fastq1.hashCode(), fastq0.hashCode());
        }
        if (fastq1.equals(fastq0))
        {
            assertEquals(fastq0.hashCode(), fastq1.hashCode());
            assertEquals(fastq1.hashCode(), fastq0.hashCode());
        }
    }
}
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
package org.biojava3.genome.parsers.fastq;

import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.util.ArrayList;
import java.util.List;
import junit.framework.TestCase;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.features.FeatureInterface;
import org.biojava3.core.sequence.features.QualityFeature;
import org.biojava3.core.sequence.template.AbstractSequence;

/**
 * Unit test for Fastq.
 */
public final class FastqTest extends TestCase {

    public void testImmutable() {
        Class<Fastq> cls = Fastq.class;
        assertTrue(Modifier.isPublic(cls.getModifiers()));
        assertTrue(Modifier.isFinal(cls.getModifiers()));
        Field[] fields = cls.getDeclaredFields();
        for (Field field : fields) {
            assertTrue(Modifier.isPrivate(field.getModifiers()));
            assertTrue(Modifier.isFinal(field.getModifiers())
                    || (Modifier.isVolatile(field.getModifiers()) && Modifier.isTransient(field.getModifiers())));
        }
    }

    public void testConstructor() {
        Fastq fastq = new Fastq("description", "sequence", "quality_", FastqVariant.FASTQ_SANGER);
        assertNotNull(fastq);

        try {
            fastq = new Fastq(null, "sequence", "quality_", FastqVariant.FASTQ_SANGER);
            fail("ctr(null description) expected IllegalArgumentException");
        } catch (IllegalArgumentException e) {
            // expected
        }
        try {
            fastq = new Fastq("description", null, "quality_", FastqVariant.FASTQ_SANGER);
            fail("ctr(null sequence) expected IllegalArgumentException");
        } catch (IllegalArgumentException e) {
            // expected
        }
        try {
            fastq = new Fastq("description", "sequence", null, FastqVariant.FASTQ_SANGER);
            fail("ctr(null quality) expected IllegalArgumentException");
        } catch (IllegalArgumentException e) {
            // expected
        }
        try {
            fastq = new Fastq("description", "sequence", "quality_", null);
            fail("ctr(null variant) expected IllegalArgumentException");
        } catch (IllegalArgumentException e) {
            // expected
        }
        try {
            fastq = new Fastq(null, FastqVariant.FASTQ_SANGER);
            fail("ctr(null DNASequence) expected IllegalArgumentException");
        } catch (IllegalArgumentException e) {
            // expected
        }
        DNASequence validseq = new DNASequence("ACGT");
        try {
            fastq = new Fastq(validseq, null);
        } catch (IllegalArgumentException e) {
            // expected
        }
        QualityFeature qual = new QualityFeature<DNASequence, NucleotideCompound>("quality", "sequencing");
        List<Number> quals = new ArrayList<Number>();
        quals.add(40);
        quals.add(30);
        quals.add(20);
        quals.add(10);
        qual.setQualities(quals);
        validseq.addFeature(1, validseq.getLength(), qual);
        try {
            fastq = new Fastq(validseq, FastqVariant.FASTQ_SANGER);
        } catch (Exception e) {
            fail("unexpected exception" + e.getMessage());
        }
    }

    public void testDescription() {
        Fastq fastq = new Fastq("description", "sequence", "quality_", FastqVariant.FASTQ_SANGER);
        assertTrue(fastq.getDescription() != null);
        assertEquals("description", fastq.getDescription());
    }

    public void testSequence() {
        Fastq fastq = new Fastq("description", "sequence", "quality_", FastqVariant.FASTQ_SANGER);
        assertTrue(fastq.getSequence() != null);
        assertEquals("sequence", fastq.getSequence());
    }

    public void testQuality() {
        Fastq fastq = new Fastq("description", "sequence", "quality_", FastqVariant.FASTQ_SANGER);
        assertTrue(fastq.getQuality() != null);
        assertEquals("quality_", fastq.getQuality());
    }

    public void testVariant() {
        Fastq fastq = new Fastq("description", "sequence", "quality_", FastqVariant.FASTQ_SANGER);
        assertTrue(fastq.getVariant() != null);
        assertEquals(FastqVariant.FASTQ_SANGER, fastq.getVariant());
    }
    
    public void testQualityConversions() {
        DNASequence seq = new DNASequence("ACGTA");
        QualityFeature qual = new QualityFeature<DNASequence, NucleotideCompound>("quality", "sequencing");
        List<Number> quals = new ArrayList<Number>();
        quals.add(40);
        quals.add(30);
        quals.add(20);
        quals.add(10);
        quals.add(1);
        qual.setQualities(quals);
        seq.addFeature(1, seq.getLength(), qual);
        Fastq fastq = new Fastq(seq, FastqVariant.FASTQ_SANGER);
        assertEquals("I?5+\"", fastq.getQuality());
        Fastq sangerfastq = new Fastq("description", "ACGTA", "I?5+\"", FastqVariant.FASTQ_SANGER);
        DNASequence dnaSequence = sangerfastq.getDNASequence();
        checkDNASequence(dnaSequence);
        Fastq illuminafastq = new Fastq(dnaSequence, FastqVariant.FASTQ_ILLUMINA);
        assertEquals("h^TJA", illuminafastq.getQuality());
        dnaSequence = illuminafastq.getDNASequence();
        checkDNASequence(dnaSequence);
        Fastq newilluminafastq = new Fastq(dnaSequence, FastqVariant.FASTQ_NEW_ILLUMINA);
        assertEquals("I?5+\"", newilluminafastq.getQuality());
        dnaSequence = newilluminafastq.getDNASequence();
        checkDNASequence(dnaSequence);
        Fastq solexafastq = new Fastq(dnaSequence, FastqVariant.FASTQ_SOLEXA);
        assertEquals("h^TJ;", solexafastq.getQuality());
        dnaSequence = solexafastq.getDNASequence();
        checkDNASequence(dnaSequence);
        sangerfastq = new Fastq(dnaSequence, FastqVariant.FASTQ_SANGER);
        assertEquals("I?5+\"", sangerfastq.getQuality());
        dnaSequence = sangerfastq.getDNASequence();
        checkDNASequence(dnaSequence);
    }

    private void checkDNASequence(DNASequence dnaSequence) {
        assertEquals("ACGTA", dnaSequence.getSequenceAsString());
        List<FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound>> features = dnaSequence.getFeaturesByType("quality");
        assertEquals(1, features.size());
        FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound> feature = features.get(0);
        assertTrue(feature instanceof QualityFeature);
        QualityFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> qualfeature = (QualityFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound>) feature;
        List<Number> qualities = qualfeature.getQualities();
        assertEquals(5, qualities.size());
        assertEquals(40, qualities.get(0));
        assertEquals(30, qualities.get(1));
        assertEquals(20, qualities.get(2));
        assertEquals(10, qualities.get(3));
        assertEquals(01, qualities.get(4));
    }

    public void testEquals() {
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

    public void testHashCode() {
        Fastq fastq0 = new Fastq("description", "sequence", "quality_", FastqVariant.FASTQ_SANGER);
        Fastq fastq1 = new Fastq("description", "sequence", "quality_", FastqVariant.FASTQ_SANGER);

        assertEquals(fastq0.hashCode(), fastq0.hashCode());
        assertEquals(fastq1.hashCode(), fastq1.hashCode());
        if (fastq0.equals(fastq1)) {
            assertEquals(fastq0.hashCode(), fastq1.hashCode());
            assertEquals(fastq1.hashCode(), fastq0.hashCode());
        }
        if (fastq1.equals(fastq0)) {
            assertEquals(fastq0.hashCode(), fastq1.hashCode());
            assertEquals(fastq1.hashCode(), fastq0.hashCode());
        }
    }
}
package org.biojava3.core;

import junit.framework.TestCase;
import org.biojava3.core.sequence.template.Sequence;
import org.biojava3.core.sequence.compound.DNACompound;
import org.biojava3.core.sequence.DNASequence;

public class DNATests extends TestCase {

    public DNATests(String testName) {
        super(testName);
    }

    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }

    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }

    public void testgetReverseComplement() {
        String s = getSeq(null).getReverseComplement().toString();
        assertEquals(s, "GCAT");

    }

    public void testgetComplement() {
        String s = getSeq(null).getComplement().toString();
        assertEquals(s, "TACG");
    }

    public void testReverse() {
        String s = getSeq(null).getReverse().toString();
        assertEquals(s, "CGTA");
    }

    public void translateToRna() {
        String s = getSeq("ATGGCGGCGCTGAGCGGT").getRNASequence().toString();
        assertEquals(s, "AUGGCGGCGCUGAGCGGU");
    }

    public void testRespectCase() {
        String s = "ATgc";
        assertEquals(s.toString(), s);
    }

    //@Test(expected=BadSequence.class)
    public void testBogusSequence() {
        getSeq("ATGCx");
    }

    public DNASequence getSeq(final String seq) {
        String target = (seq == null) ? "ATGC" : seq;
        return new DNASequence(target);
    }
}


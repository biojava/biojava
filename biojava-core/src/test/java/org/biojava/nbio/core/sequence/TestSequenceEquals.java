package org.biojava.nbio.core.sequence;

import org.junit.Test;

import static org.junit.Assert.*;

/**
 * Created by yana on 3/27/17.
 */
public class TestSequenceEquals {

    @Test
    public void testSameCompounds() throws Exception{



        ProteinSequence seq1 = new ProteinSequence("ARNDCEQGHILKMFPSTWYVBZJX");

        ProteinSequence seq2 = new ProteinSequence("ARNDCEQGHILKMFPSTWYVBZJXARNDCEQGHILKMFPSTWYVBZJX");


        assertFalse(seq1.equals(seq2));

        assertTrue(seq1.equals(seq1));

        assertTrue(seq2.equals(seq2));


        ProteinSequence seq3 = new ProteinSequence("ARNDCEQGHILKMFPSTWYVBZJX");

        assertTrue(seq3.equals(seq1));

        assertFalse(seq2.equals(seq3));


        DNASequence dnaSeq = new DNASequence("ATGGCGGCGCTGAGCGGT");

        assertFalse(seq1.equals(dnaSeq));



    }
}

package org.biojava.nbio.core.sequence;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class GeneSequenceTest {

    @Nested
    class AfterValidConstruction {
        ChromosomeSequence chromosomeSequence;
        GeneSequence geneSequence;

        @BeforeEach
        void before() throws CompoundNotFoundException {
            chromosomeSequence = new ChromosomeSequence(ChromosomeSequenceTest.CHROMOSOME_SEQ);
            geneSequence = new GeneSequence(chromosomeSequence, 10,19, Strand.POSITIVE);
        }

        @Test
        void lengthIsSetByBeginAndEnd() {
            assertEquals(10, geneSequence.getLength());
        }

        @Test
        void noExonsIntronsOrTranscripts() {
            assertEquals(0, geneSequence.getExonSequences().size());
            assertEquals(0, geneSequence.getIntronSequences().size());
            assertEquals(0, geneSequence.getTranscripts().size());
        }
        @Test
        void geneSequenceIsChromosomeSequence() {
            assertEquals(chromosomeSequence.getSequenceAsString(), geneSequence.getSequenceAsString());
        }
    }


    @Test
    void addIntronsUsingExons() {
    }

    @Test
    void addTranscript() {
    }

    @Test
    void removeExon() {
    }

    @Test
    void addExon() {
    }

    @Test
    void getSequence5PrimeTo3Prime() {
    }
}
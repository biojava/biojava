package org.biojava.nbio.core.sequence;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class GeneSequenceTest {
    GeneSequence geneSequence;
    ChromosomeSequence chromosomeSequence;
    @BeforeEach
    void before() throws CompoundNotFoundException {
        chromosomeSequence = new ChromosomeSequence(ChromosomeSequenceTest.CHROMOSOME_SEQ);
        geneSequence = new GeneSequence(chromosomeSequence, 10,19, Strand.POSITIVE);
    }

    @Nested
    class AfterValidConstruction {


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

        @Test
        void geneSequenceIsChromos2omeSequence() {
            geneSequence.setStrand(Strand.NEGATIVE);
            geneSequence.getSequenceAsString();
            geneSequence.getBioBegin();
        }
    }


    @Test
    void addIntronsUsingExonsPositiveStrand() throws Exception {
        geneSequence = new GeneSequence(chromosomeSequence, 10,150, Strand.POSITIVE);
        geneSequence.addExon( new AccessionID("a"), 10,29);
        geneSequence.addExon( new AccessionID("b"), 33,80);
        geneSequence.addExon( new AccessionID("c"), 100,120);
        geneSequence.setAccession(new AccessionID("geneId"));
        geneSequence.addIntronsUsingExons();
        assertEquals(2, geneSequence.getIntronSequences().size());
        assertEquals(30, geneSequence.getIntronSequences().get(0).getBioBegin());
        assertEquals(32, geneSequence.getIntronSequences().get(0).getBioEnd());
        assertEquals(81, geneSequence.getIntronSequences().get(1).getBioBegin());
        assertEquals(99, geneSequence.getIntronSequences().get(1).getBioEnd());
    }

    @Test
    @Disabled("gives odd results for intron coords")
    void addIntronsUsingExonsNegativeStrand() throws Exception {
        geneSequence = new GeneSequence(chromosomeSequence, 150,10, Strand.NEGATIVE);
        ExonSequence e1 = geneSequence.addExon( new AccessionID("c"), 120,100);

        geneSequence.addExon( new AccessionID("b"), 80,33);
        geneSequence.addExon( new AccessionID("a"), 29,10);

        // this MUST be set in order to avoid NPE when adding introns
        geneSequence.setAccession(new AccessionID("geneId"));
        geneSequence.addIntronsUsingExons();
        // actual values generated are (9,81) for I1 and (32,121) for I2
        assertEquals(2, geneSequence.getIntronSequences().size());
        assertEquals(99, geneSequence.getIntronSequences().get(0).getBioBegin());
        assertEquals(81, geneSequence.getIntronSequences().get(0).getBioEnd());
        assertEquals(32, geneSequence.getIntronSequences().get(1).getBioBegin());
        assertEquals(30, geneSequence.getIntronSequences().get(1).getBioEnd());
    }

    @Test
    void getPositiveStrandSequence5To3Prime() {
        geneSequence = new GeneSequence(chromosomeSequence, 10,150, Strand.POSITIVE);
        // this must be set to avoid NPE
        geneSequence.setAccession(new AccessionID("geneId"));
        assertEquals(chromosomeSequence.getSequenceAsString().substring(9,150),
                geneSequence.getSequence5PrimeTo3Prime().getSequenceAsString());
    }

    @Test
    @Disabled("not complementing - seems to complement twice???")
    void getNegativeStrandSequence5To3Prime() {
        geneSequence = new GeneSequence(chromosomeSequence, 5,15, Strand.NEGATIVE);
        // this must be set to avoid NPE
        geneSequence.setAccession(new AccessionID("geneId"));
        DNASequence seq = geneSequence.getSequence5PrimeTo3Prime();
        System.err.println(  geneSequence.getSequence5PrimeTo3Prime().getSequenceAsString());
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
package org.biojava.nbio.core.sequence;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class TranscriptSequenceTest {
    GeneSequence anyGeneSequence;
    GeneSequence anyNegativeGeneSequence;
    TranscriptSequence transcriptSeq;
    TranscriptSequence transcriptNegativeSeq;


    @BeforeEach
    void setUp() throws CompoundNotFoundException {
        anyGeneSequence = SequenceTestUtils.anyGeneSequence();
        transcriptSeq = new TranscriptSequence(anyGeneSequence, new AccessionID("T5"), 5, 100);
        anyNegativeGeneSequence = SequenceTestUtils.any3GeneSequence();
        transcriptNegativeSeq = new TranscriptSequence(anyNegativeGeneSequence, new AccessionID("T3"), 5, 100);
    }

    @Nested
    class AfterValidConstruction {
        @Test
        void lengthIsTranscriptLength() {
            assertEquals(96, transcriptSeq.getLength());
        }

        @Test
        void strandIsSameAsGene() {
            assertEquals(anyGeneSequence.getStrand(), transcriptSeq.getStrand());
            assertEquals(anyNegativeGeneSequence.getStrand(), transcriptNegativeSeq.getStrand());
        }

        @Test
        void CDSListIsEmpty() {
            assertEquals(0, transcriptSeq.getCDSSequences().size());
        }

        @Test
        void equals() {
            assertTrue(transcriptSeq.equals(transcriptSeq));
        }

        // whether it's -ve or +ve doesn't affect equals?
        void equalsDoesntDependOnStrand() {
            assertTrue(transcriptSeq.equals(transcriptNegativeSeq));
        }

        @Test
        void hashcode() {
            assertTrue(transcriptSeq.hashCode() == (transcriptNegativeSeq.hashCode()));
        }
    }

    @Test
    void addCDS() throws Exception {
        transcriptSeq.addCDS(new AccessionID("b"), 40, 50, 1);
        assertEquals(1, transcriptSeq.getCDSSequences().size());
    }

    @Test
    void getCDNASeqPositiveStrand() throws Exception {
        String chrSeq = ChromosomeSequenceTest.CHROMOSOME_SEQ;
        // must set this to avoid NPE when generating sequence


        // make 2 CDS that are contiguous. These can be added in any order and are sorted OK
        CDSSequence s1 = transcriptSeq.addCDS(new AccessionID("a"), 11, 20, 1);
        assertEquals(chrSeq, s1.getSequenceAsString());

        CDSSequence s2 = transcriptSeq.addCDS(new AccessionID("b"), 1, 10, 1);
        assertEquals(chrSeq, s2.getSequenceAsString());

        DNASequence cDNA = transcriptSeq.getDNACodingSequence();
        assertEquals(chrSeq.substring(0, 20), cDNA.getSequenceAsString());
        assertEquals(20, cDNA.getLength());
    }

    @Test
    @Disabled("is reversed, not complemented?")
    void getCDNASeqNegativeStrand() throws Exception {
        TranscriptSequence ts = SequenceTestUtils.transcriptFromSequence("AAAAACCCCCTTTTGGGGGG", 3, 10, Strand.NEGATIVE);
        CDSSequence s2 = ts.addCDS(new AccessionID("b"), 1, 10, 0);
        // this should be GGGGGTTTTT( ie the reverse complement of the chromosome sequence,
        // but is just reversed and generates CCCCCAAAAA
        //assertEquals("GGGGGTTTTT", ts.getDNACodingSequence());
    }

    @Test
    void removeCDS() throws Exception {
        transcriptSeq.addCDS(new AccessionID("a"), 50, 60, 1);
        assertEquals(1, transcriptSeq.getCDSSequences().size());
        // throws NPE
        transcriptSeq.removeCDS("a");
        assertEquals(0, transcriptSeq.getCDSSequences().size());
    }

    @Test
    void addGetStartCodonSequence () {
        assertNull(transcriptSeq.getStartCodonSequence());
        transcriptSeq.addStartCodonSequence(new AccessionID("cds"), 40,42);
        StartCodonSequence scs = transcriptSeq.getStartCodonSequence();
        assertEquals(3, scs.getLength());
    }

    @Test
    void addGetStopCodonSequence () {
        assertNull(transcriptSeq.getStopCodonSequence());
        transcriptSeq.addStopCodonSequence(new AccessionID("cds"), 40,42);
        StopCodonSequence scs = transcriptSeq.getStopCodonSequence();
        assertEquals(3, scs.getLength());
    }
}
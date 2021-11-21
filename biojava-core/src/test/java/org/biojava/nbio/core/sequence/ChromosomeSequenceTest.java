package org.biojava.nbio.core.sequence;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.EnumSource;
import org.junit.jupiter.params.provider.ValueSource;

import static org.junit.jupiter.api.Assertions.*;

class ChromosomeSequenceTest {

    static final String CHROMOSOME_SEQ = "ATATCGACTTATATATATATATATATATATATATACGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCATATATATATATATATATATATACGCGCGCGCGCGCGCGCATATATATATATATATATATATATATATATATACGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCATATATATATATATATATATATACGCGCGCGCGCGCGCGC";

    @Nested
    class AfterValidConstruction {
        ChromosomeSequence seq = null;

        @BeforeEach
        void before() throws CompoundNotFoundException {
            seq = new ChromosomeSequence(CHROMOSOME_SEQ);
        }

        @Test
        void beginAndEndAreLengthOfSequence() {
            assertEquals(1, seq.getBioBegin());
            assertEquals(210, seq.getBioEnd());
            assertEquals(210, seq.getLength());
        }

        @Test
        void noGenesAreDefined() {
            assertEquals(0, seq.getGeneSequences().size());
        }

        @Test
        void chromosomeNumberIsZero() {
            assertEquals(0, seq.getChromosomeNumber());
        }

        @Test
        void sequenceTypeIsUnknown() {
            assertEquals(DNASequence.DNAType.UNKNOWN, seq.getDNAType());
        }
    }

    @Nested
    class AfterConstructionWithEmptyString {
        ChromosomeSequence seq = null;

        @BeforeEach
        void before() throws CompoundNotFoundException {
            seq = new ChromosomeSequence("");
        }

        @Test
        void lengthIsZero() {
            assertEquals(0, seq.getLength());
        }

        @Test
        void endIsBeforeBeginning() {
            assertEquals(0, seq.getBioEnd());
            assertEquals(1, seq.getBioBegin());
        }

    }

    @Test
    void nullSequenceNotAllowed() throws CompoundNotFoundException {
        assertThrows(NullPointerException.class, () -> new ChromosomeSequence((String) null));
    }

    @ParameterizedTest
    @ValueSource(ints = {Integer.MAX_VALUE, Integer.MIN_VALUE, 100, 0, -1, -100})
    void anyIntegerIsValidChromosomeNumber(int value) throws CompoundNotFoundException {
        ChromosomeSequence seq = new ChromosomeSequence(CHROMOSOME_SEQ);
        seq.setChromosomeNumber(value);
        assertEquals(value, seq.getChromosomeNumber());
    }

    @ParameterizedTest
    @EnumSource(DNASequence.DNAType.class)
    void anyDNATypeIsValid(DNASequence.DNAType dnaType) throws CompoundNotFoundException {
        ChromosomeSequence seq = new ChromosomeSequence(CHROMOSOME_SEQ);
        seq.setDNAType(dnaType);
        assertEquals(dnaType, seq.getDNAType());
    }

    @Nested
    class AddingAndRemovingGeneSequences {
        ChromosomeSequence seq = null;
        @BeforeEach
        void before() throws CompoundNotFoundException {
            seq = new ChromosomeSequence(CHROMOSOME_SEQ);
        }
        @Test
        void canAddSameGeneTwice(){
            seq.addGene(new AccessionID("ABCDE1"), 1, 20, Strand.POSITIVE);
            assertEquals(1, seq.getGeneSequences().size());
            seq.addGene(new AccessionID("ABCDE1"), 1, 20, Strand.POSITIVE);
            assertEquals(1, seq.getGeneSequences().size());
        }

        @Test
        void isOKToRemoveNonExistentSequence(){
            seq.removeGeneSequence("XXX");
        }

        @Test
        void addAndRemove(){
            final String accessionId = "ABCDE1";
            GeneSequence geneSequence = seq.addGene(new AccessionID(accessionId), 1, 20, Strand.POSITIVE);
            assertEquals(geneSequence.getAccession(), seq.getGene(accessionId).getAccession());
            assertEquals(1, seq.getGeneSequences().size());
            seq.removeGeneSequence(accessionId);
            assertEquals(0, seq.getGeneSequences().size());
        }

        @Test
        void geneSequenceHasCorrectLength(){
            final String accessionId = "ABCDE1";
            GeneSequence geneSequence = seq.addGene(new AccessionID(accessionId), 1, 20, Strand.POSITIVE);
            assertEquals(20, geneSequence.getLength());
        }

        @Test
        void geneSequenceCanHaveBeginAndEndOutsideOfChromosomeSeq(){
            final String accessionId = "ABCDE1";
            GeneSequence geneSequence = seq.addGene(new AccessionID(accessionId), Integer.MAX_VALUE-10, Integer.MAX_VALUE, Strand.POSITIVE);
            assertEquals(11, geneSequence.getLength());

        }
    }

    @Test
    void addAndRemoveGeneSequence() throws CompoundNotFoundException {
        ChromosomeSequence seq = new ChromosomeSequence(CHROMOSOME_SEQ);
        seq.addGene(new AccessionID("ABCDE1"), 1, 20, Strand.POSITIVE);
        assertEquals(1, seq.getGeneSequences().size());


        // still present
        assertEquals(1, seq.getGeneSequences().size());
        // can be added again with sam
        seq.addGene(new AccessionID("ABCDE1"), 1, 20, Strand.POSITIVE);
        assertEquals(1, seq.getGeneSequences().size());

    }



    @Test
    void addGene() {
    }

    @Test
    void getGene() {
    }
}
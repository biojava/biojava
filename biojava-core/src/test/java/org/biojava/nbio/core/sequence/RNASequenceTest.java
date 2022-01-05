package org.biojava.nbio.core.sequence;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.SequenceView;
import org.biojava.nbio.core.sequence.transcription.TranscriptionEngine;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import static java.util.stream.Collectors.joining;
import static org.junit.jupiter.api.Assertions.*;

class RNASequenceTest {

    // AUG start, then 3 AA, then stop codon
    final String rnaSeq    = "AUGGUCGAACUCUGA";
    final String rnaSeqCompl = "UACCAGCUUGAGACU";
    final String rnaSeqReversed = "AGUCUCAAGCUGGUA";
    final String rnaSeqReversedComplement = "UCAGAGUUCGACCAU";
    RNASequence rna;
    @BeforeEach
    void before() throws CompoundNotFoundException {
        rna = new RNASequence(rnaSeq);
    }

    @Test
    void translateToProteinSequence()  {
        ProteinSequence protein =  rna.getProteinSequence(TranscriptionEngine.getDefault());
        assertEquals(4, protein.getLength());
        assertEquals("MVEL", protein.getSequenceAsString());
    }

    @Test
    void complement()  {
        SequenceView<NucleotideCompound> complement = rna.getComplement();
        assertEquals(rnaSeqCompl, complement.getSequenceAsString());
        assertEquals(rnaSeq, complement.getViewedSequence().getSequenceAsString());
        assertEquals(rna.getLength(), complement.getLength());
    }

    @Test
    void reverse()  {
        SequenceView<NucleotideCompound> reversed = rna.getInverse();
        assertEquals(rnaSeqReversed, reversed.getSequenceAsString());
        assertEquals(rnaSeq, reversed.getViewedSequence().getSequenceAsString());
        assertEquals(rna.getLength(), reversed.getLength());
    }

    @Test
    void reverseComplement()  {
        SequenceView<NucleotideCompound> reverseComplement = rna.getReverseComplement();
        assertEquals(rnaSeqReversedComplement, reverseComplement.getSequenceAsString());
        assertEquals(rna.getLength(), reverseComplement.getLength());
        StringBuilder sb = new StringBuilder();
        for (int i = 1; i <= rna.getLength(); i++) {
            sb.append(reverseComplement.getCompoundAt(i).toString());
        }
        assertEquals(rnaSeqReversedComplement, sb.toString());

        sb = new StringBuilder();
        for (Compound c: reverseComplement) {
            sb.append(c.toString());
        }
        assertEquals(rnaSeqReversedComplement, sb.toString());
        assertEquals(rnaSeqReversedComplement, reverseComplement.getAsList().stream().map(Compound::toString).collect(joining("")));
    }

    @Test
    void rejectThymineInSequence()  {
        String dna = rnaSeq.replaceAll("U", "T");
        assertThrows(CompoundNotFoundException.class, ()->new RNASequence(dna));
    }

}
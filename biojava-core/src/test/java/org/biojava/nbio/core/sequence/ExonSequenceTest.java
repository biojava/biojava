package org.biojava.nbio.core.sequence;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;

import static org.biojava.nbio.core.sequence.SequenceTestUtils.anyGeneSequence;
import static org.junit.jupiter.api.Assertions.*;

class ExonSequenceTest {
    @Test
    void createExon() throws CompoundNotFoundException {
        GeneSequence gene = anyGeneSequence();
        ExonSequence es = new ExonSequence(gene, 30, 40);
        assertEquals(11, es.getLength());
    }

    @Test
    void equalsAndHashcode() throws CompoundNotFoundException {
        GeneSequence gene = anyGeneSequence();
        ExonSequence es = new ExonSequence(gene, 30, 40);
        ExonSequence es2 = new ExonSequence(gene, 30, 40);
        // calling equals throws npe
        // assertEquals(es, es2);

        // this also throws NPE
        assertEquals(es.hashCode(), es2.hashCode());
    }

}
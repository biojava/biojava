package org.biojava.nbio.core.sequence;

import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import static org.biojava.nbio.core.sequence.SequenceTestUtils.*;
import static org.junit.jupiter.api.Assertions.*;

class ExonComparatorTest {

    @Test
    void sortPositiveStrandExons() throws Exception {
        GeneSequence geneSequence = anyGeneSequence();
        // added in order 2,3,1
        ExonSequence e2 = geneSequence.addExon(new AccessionID("b"), 40, 60);
        ExonSequence e3 = geneSequence.addExon(new AccessionID("c"), 80, 100);
        ExonSequence e1 = geneSequence.addExon(new AccessionID("a"), 10, 30);
        List<ExonSequence> exonsToSort = new ArrayList<>();
        exonsToSort.add(e2);
        exonsToSort.add(e3);
        exonsToSort.add(e1);
        Collections.sort(exonsToSort, new ExonComparator());
        // sorted by starting position, in 5' to 3' order
        assertEquals("a", exonsToSort.get(0).getAccession().getID());
        assertEquals("b", exonsToSort.get(1).getAccession().getID());
        assertEquals("c", exonsToSort.get(2).getAccession().getID());
    }

    @Test
    void sortNegativeStrandExons() throws Exception {
        GeneSequence geneSequence = anyGeneSequence();
        geneSequence.setStrand(Strand.NEGATIVE);
        // added in order 2,3,1
        ExonSequence e2 = geneSequence.addExon(new AccessionID("b"), 60, 40);
        ExonSequence e3 = geneSequence.addExon(new AccessionID("c"), 100, 80);
        ExonSequence e1 = geneSequence.addExon(new AccessionID("a"), 30, 10);
        List<ExonSequence> exonsToSort = new ArrayList<>();
        exonsToSort.add(e2);
        exonsToSort.add(e3);
        exonsToSort.add(e1);
        Collections.sort(exonsToSort, new ExonComparator());
        // sorted by starting position - this is 3' - 5' order
        assertEquals("a", exonsToSort.get(0).getAccession().getID());
        assertEquals("b", exonsToSort.get(1).getAccession().getID());
        assertEquals("c", exonsToSort.get(2).getAccession().getID());
    }

}
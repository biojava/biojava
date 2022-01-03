package org.biojava.nbio.core.sequence;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;

/**
 * Test utility methods for classes in this package
 */
public class SequenceTestUtils {

    /**
     * A gene sequence of 190 bp length on + strand
     *
     * @return
     * @throws CompoundNotFoundException
     */
    static GeneSequence anyGeneSequence() throws CompoundNotFoundException {
        ChromosomeSequence chr = new ChromosomeSequence(ChromosomeSequenceTest.CHROMOSOME_SEQ);
        return new GeneSequence(chr, new AccessionID("someGeneId"), 10, 200, Strand.POSITIVE);
    }

    /**
     * A gene sequence of 190 bp length on MINUS strand
     *
     * @return
     * @throws CompoundNotFoundException
     */
    static GeneSequence any3GeneSequence() throws CompoundNotFoundException {
        ChromosomeSequence chr = new ChromosomeSequence(ChromosomeSequenceTest.CHROMOSOME_SEQ);
        GeneSequence gene = new GeneSequence(chr, new AccessionID("some3PrimeGeneId"),10, 200, Strand.NEGATIVE);
        return gene;
    }

    /**
     * Generate a GeneSequence as a subsequence of defined chromosome sequence.
     *
     * @param chromosomeSequence
     * @param bioStart
     * @param bioEnd
     * @param strand
     * @return
     * @throws CompoundNotFoundException
     */
    static GeneSequence fromSequence(String chromosomeSequence, int bioStart, int bioEnd, Strand strand) throws CompoundNotFoundException {
        ChromosomeSequence chr = new ChromosomeSequence(chromosomeSequence);
        GeneSequence gene = new GeneSequence(chr, new AccessionID("Gene"), bioStart, bioEnd, strand);
        gene.setAccession(new AccessionID("Gene1"));
        return gene;
    }

    /**
     * Creates a transcript from coordinates on the supplied chromosome sequence.
     * The GeneSequence is set to same length as Chromosomal sequence for simplicity.
     *
     * @param chromosomeSequence
     * @param bioStart
     * @param bioEnd
     * @param strand
     * @return
     * @throws CompoundNotFoundException
     */
    static TranscriptSequence transcriptFromSequence(String chromosomeSequence, int bioStart, int bioEnd, Strand strand) throws CompoundNotFoundException {
        GeneSequence gene = fromSequence(chromosomeSequence, 1, chromosomeSequence.length(), strand);
        TranscriptSequence ts = new TranscriptSequence(gene, new AccessionID("Transcript"), bioStart, bioEnd);
        ts.setAccession(new AccessionID("Transcript1"));
        return ts;
    }
}

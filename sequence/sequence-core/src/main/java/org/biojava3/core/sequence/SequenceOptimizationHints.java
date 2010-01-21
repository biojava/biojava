/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.core.sequence;

/**
 * A static class that provides optimization hints for memory or performance handling of sequence data.
 * If you are working with genomic sequence data and the use case is sub sequence then the sequence proxy loader
 * implementation may do file reads for each sub-sequence request instead of loading everything into memory.
 * If a large collection of protein sequences is being loaded from a fasta file but only a few sequences will be selected
 * then the file loader could create a sequence proxy loader that retains the offset in the file for each sequence
 * and it is loaded when the sequence data is requested. This way you could load a pointer to all sequences but delay
 * loading until the user actually needs the data. The sequence loader could also have an internal sequence management
 * algorithm that goes through and returns sequence data freeing up memory where a future request for sequence data will
 * then be reloaded.
 *
 * @author Scooter
 */
public class SequenceOptimizationHints {

    /**
     * @return the sequenceUsage
     */
    public static SequenceUsage getSequenceUsage() {
        return sequenceUsage;
    }

    /**
     * @param aSequenceUsage the sequenceUsage to set
     */
    public static void setSequenceUsage(SequenceUsage aSequenceUsage) {
        sequenceUsage = aSequenceUsage;
    }

    /**
     * @return the sequenceColection
     */
    public static SequenceCollection getSequenceCollection() {
        return sequenceCollection;
    }

    /**
     * @param aSequenceColection the sequenceColection to set
     */
    public static void setSequenceCollection(SequenceCollection aSequenceColection) {
        sequenceCollection = aSequenceColection;
    }

    public enum SequenceUsage {

        FULL_SEQUENCE_DATA, SUB_SEQUENCE_DATA, MINIMAL_SEQUENCE_DATA;
    }

    public enum SequenceCollection {

        ALL_SEQUENCES, VARIABLE_SEQUENCES, MINIMINAL_SEQUENCES;
    }

    static private SequenceUsage sequenceUsage = SequenceUsage.FULL_SEQUENCE_DATA;
    static private SequenceCollection sequenceCollection = SequenceCollection.ALL_SEQUENCES;



    

}

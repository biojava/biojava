/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.core.sequence;

/**
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

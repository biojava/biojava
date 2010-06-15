/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 * Created on June 7, 2010
 * Author: Mark Chapman
 */

package org.biojava3.alignment;

import java.util.ArrayList;
import java.util.List;

import org.biojava3.alignment.template.*;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.Sequence;

public class Alignments {

    // prevents instantiation
    private Alignments() { }

    public static <S extends Sequence<C>, C extends Compound>
            List<SequencePair<S, C>> getAllPairsAlignments(List<S> sequences) {
        List<SequencePair<S, C>> allPairs = new ArrayList<SequencePair<S, C>>();
        for (int i = 0; i < sequences.size(); i++) {
            for (int j = i+1; j < sequences.size(); j++) {
                allPairs.add(getPairwiseAlignment(sequences.get(i), sequences.get(j)));
            }
        }
        return allPairs;
    }

    public static <S extends Sequence<C>, C extends Compound>
            int[] getAllPairsScores(List<S> sequences) {
        int[] allPairs = new int[sequences.size()*(sequences.size()-1)/2];
        for (int i = 0, p = 0; i < sequences.size(); i++) {
            for (int j = i+1; p < allPairs.length && j < sequences.size(); j++, p++) {
                allPairs[p] = getPairwiseScore(sequences.get(i), sequences.get(j));
            }
        }
        return allPairs;
    }

    public static <S extends Sequence<C>, C extends Compound>
            SequencePair<S, C> getPairwiseAlignment(S sequence1, S sequence2) {
        // TODO
        return null;
    }

    public static <S extends Sequence<C>, C extends Compound>
            int getPairwiseScore(S sequence1, S sequence2) {
        // TODO
        return 0;
    }

    public static <S extends Sequence<C>, C extends Compound>
            Profile<S, C> getMultipleSequenceAlignment(List<S> sequences) {
        // TODO
        return null;
    }

    public static enum MSAEmulation { CLUSTALW, MUSCLE, KALIGN, CUSTOM }

    /**
     * Stores the default values for the alignment algorithms and data structures
     */
    public static class Default {

        private static MSAEmulation emulation;
        private static Class<? extends GapPenalty> gapPenalty;
        private static Class<? extends HierarchicalClusterer> clusterer;
        private static Class<? extends PairwiseSequenceAligner<?, ?>> pwAligner;
        private static Class<? extends PairwiseSequenceScorer<?, ?>> pwScorer;
        private static Class<? extends PartitionRefiner<?, ?>> pRefiner;
        private static Class<? extends ProfileProfileAligner<?, ?>> ppAligner;
        private static Class<? extends RescoreRefiner<?, ?>> rRefiner;
        private static Class<? extends SubstitutionMatrix<?>> subMatrix;

        public static MSAEmulation getEmulation() {
            return Default.emulation;
        }

        public static Class<? extends GapPenalty> getGapPenalty() {
            return Default.gapPenalty;
        }

        public static Class<? extends HierarchicalClusterer> getHierarchicalClusterer() {
            return Default.clusterer;
        }

        public static Class<? extends PairwiseSequenceAligner<?, ?>> getPairwiseSequenceAligner() {
            return Default.pwAligner;
        }

        public static Class<? extends PairwiseSequenceScorer<?, ?>> getPairwiseSequenceScorer() {
            return Default.pwScorer;
        }

        public static Class<? extends PartitionRefiner<?, ?>> getPartitionRefiner() {
            return Default.pRefiner;
        }

        public static Class<? extends ProfileProfileAligner<?, ?>> getProfileProfileAligner() {
            return Default.ppAligner;
        }

        public static Class<? extends RescoreRefiner<?, ?>> getRescoreRefiner() {
            return Default.rRefiner;
        }

        public static Class<? extends SubstitutionMatrix<?>> getSubstitutionMatrix() {
            return Default.subMatrix;
        }

        public static void setEmulation(MSAEmulation emulation) {
            // TODO set defaults according to emulation value
            Default.emulation = MSAEmulation.CUSTOM;
        }

        public static void setGapPenalty(Class<? extends GapPenalty> gapPenalty) {
            Default.gapPenalty = gapPenalty;
        }

        public static void setHierarchicalClusterer(Class<? extends HierarchicalClusterer> clusterer) {
            Default.clusterer = clusterer;
        }

        public static void setPairwiseSequenceAligner(Class<? extends PairwiseSequenceAligner<?, ?>> pwAligner) {
            Default.pwAligner = pwAligner;
        }

        public static void setPairwiseSequenceScorer(Class<? extends PairwiseSequenceScorer<?, ?>> pwScorer) {
            Default.pwScorer = pwScorer;
        }

        public static void setPartitionRefiner(Class<? extends PartitionRefiner<?, ?>> pRefiner) {
            Default.pRefiner = pRefiner;
        }

        public static void setProfileProfileAligner(Class<? extends ProfileProfileAligner<?, ?>> ppAligner) {
            Default.ppAligner = ppAligner;
        }

        public static void setRescoreRefiner(Class<? extends RescoreRefiner<?, ?>> rRefiner) {
            Default.rRefiner = rRefiner;
        }

        public static void setSubstitutionMatrix(Class<? extends SubstitutionMatrix<?>> subMatrix) {
            Default.subMatrix = subMatrix;
        }

    }

}

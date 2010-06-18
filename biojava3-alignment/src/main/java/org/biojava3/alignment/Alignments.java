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
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

import org.biojava3.alignment.template.*;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.Sequence;
import org.biojava3.core.util.ConcurrencyTools;

/**
 * Static utility to easily run alignment routines.
 *
 * @author Mark Chapman
 */
public class Alignments {

    // prevents instantiation
    private Alignments() { }

    /**
     * Factory method to run a list of alignments concurrently.
     *
     * @param <A> each {@link Aligner} is of type A
     * @param <S> each {@link Sequence} of an alignment pair is of type S
     * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
     * @param alignments
     * @return
     */
    public static <A extends PairwiseSequenceAligner<S, C>, S extends Sequence<C>, C extends Compound>
            List<SequencePair<S, C>> getAllPairsAlignments(List<A> alignments) {
        // submit alignment tasks to the shared thread pool
        int n = 1, all = alignments.size();
        List<Future<SequencePair<S, C>>> allPairsFutures = new ArrayList<Future<SequencePair<S, C>>>();
        for (PairwiseSequenceAligner<S, C> alignment : alignments) {
            allPairsFutures.add(ConcurrencyTools.submit(new CallablePairwiseSequenceAligner<S, C>(alignment),
                    String.format("Aligning pair %d of %d", n++, all)));
        }
        // store completed alignments
        List<SequencePair<S, C>> allPairs = new ArrayList<SequencePair<S, C>>();
        for (Future<SequencePair<S, C>> f : allPairsFutures) {
            // TODO when added to ConcurrencyTools, log completions and exceptions instead of printing stack traces
            try {
                allPairs.add(f.get());
            } catch (InterruptedException e) {
                e.printStackTrace();
            } catch (ExecutionException e) {
                e.printStackTrace();
            }
        }
        return allPairs;
    }

    /**
     * Factory method which computes a global sequence alignment for all {@link Sequence} pairs in the given
     * {@link List}.  This method runs the alignments in parallel by submitting all of the alignments to the shared
     * thread pool of the {@link ConcurrencyTools} utility.
     *
     * @param <S> each {@link Sequence} of the alignment pair is of type S
     * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
     * @param sequences the {@link List} of {@link Sequence}s to align
     * @param gapPenalty the gap penalties used during alignment
     * @param subMatrix the set of substitution scores used during alignment
     * @return list of sequence alignment pairs
     */
    public static <S extends Sequence<C>, C extends Compound> List<SequencePair<S, C>> getAllPairsGlobalAlignments(
            List<S> sequences, GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix) {
        List<NeedlemanWunsch<S, C>> allPairs = new ArrayList<NeedlemanWunsch<S, C>>();
        for (int i = 0; i < sequences.size(); i++) {
            for (int j = i+1; j < sequences.size(); j++) {
                allPairs.add(new NeedlemanWunsch<S, C>(sequences.get(i), sequences.get(j), gapPenalty, subMatrix));
            }
        }
        return getAllPairsAlignments(allPairs);
    }

    public static <S extends Sequence<C>, C extends Compound>
            int[] getAllPairsScores(List<S> sequences) {
        // TODO similar to getAllPairsAlignments
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

    /* TODO abandoned idea?
    public static enum MSAEmulation { CLUSTALW, MUSCLE, KALIGN, CUSTOM }

    /**
     * Stores the default values for the alignment algorithms and data structures
     *
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

    } */

}

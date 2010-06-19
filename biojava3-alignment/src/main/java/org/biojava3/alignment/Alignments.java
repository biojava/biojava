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

    /**
     * List of multiple sequence alignment routines that can (mostly) be emulated.
     */
    public static enum MSAEmulation {
        CLUSTALW,
        MUSCLE,
        KALIGN
    }

    /**
     * List of implemented pairwise sequence alignment routines.
     */
    public static enum PairwiseAligner {
        GLOBAL,
        LOCAL
    }

    /**
     * List of implemented pairwise sequence scoring routines.
     */
    public static enum PairwiseScorer {
        GLOBAL,
        GLOBAL_IDENTITIES,
        GLOBAL_SIMILARITIES,
        LOCAL,
        LOCAL_IDENTITIES,
        LOCAL_SIMILARITIES,
        KMERS
    }

    // prevents instantiation
    private Alignments() { }

    /**
     * Factory method which computes a sequence alignment for all {@link Sequence} pairs in the given {@link List}.
     * This method runs the alignments in parallel by submitting all of the alignments to the shared thread pool of the
     * {@link ConcurrencyTools} utility.
     *
     * @param <S> each {@link Sequence} of an alignment pair is of type S
     * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
     * @param sequences the {@link List} of {@link Sequence}s to align
     * @param type chosen type from list of pairwise sequence alignment routines
     * @param gapPenalty the gap penalties used during alignment
     * @param subMatrix the set of substitution scores used during alignment
     * @return list of sequence alignment pairs
     */
    public static <S extends Sequence<C>, C extends Compound> List<SequencePair<S, C>> getAllPairsAlignments(
            List<S> sequences, PairwiseAligner type, GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix) {
        return runAligners(getAlignerList(sequences, type, gapPenalty, subMatrix));
    }

    /**
     * Factory method which computes a sequence pair score for all {@link Sequence} pairs in the given {@link List}.
     * This method runs the scorings in parallel by submitting all of the scorings to the shared thread pool of the
     * {@link ConcurrencyTools} utility.
     *
     * @param <S> each {@link Sequence} of a pair is of type S
     * @param <C> each element of a {@link Sequence} is a {@link Compound} of type C
     * @param sequences the {@link List} of {@link Sequence}s to align
     * @param type chosen type from list of pairwise sequence scoring routines
     * @param gapPenalty the gap penalties used during alignment
     * @param subMatrix the set of substitution scores used during alignment
     * @return list of sequence pair scores
     */
    public static <S extends Sequence<C>, C extends Compound> int[] getAllPairsScores(
            List<S> sequences, PairwiseScorer type, GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix) {
        return runScorers(getScorerList(sequences, type, gapPenalty, subMatrix));
    }

    /**
     * Factory method which computes a sequence alignment for the given {@link Sequence} pair.
     *
     * @param <S> each {@link Sequence} of the pair is of type S
     * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
     * @param query the first {@link Sequence}s to align
     * @param target the second {@link Sequence}s to align
     * @param type chosen type from list of pairwise sequence alignment routines
     * @param gapPenalty the gap penalties used during alignment
     * @param subMatrix the set of substitution scores used during alignment
     * @return sequence alignment pair
     */
    public static <S extends Sequence<C>, C extends Compound> SequencePair<S, C> getPairwiseAlignment(
            S query, S target, PairwiseAligner type, GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix) {
        return getAligner(query, target, type, gapPenalty, subMatrix).getPair();
    }

    /**
     * Factory method which computes a similarity score for the given {@link Sequence} pair.
     *
     * @param <S> each {@link Sequence} of the pair is of type S
     * @param <C> each element of a {@link Sequence} is a {@link Compound} of type C
     * @param query the first {@link Sequence} to score
     * @param target the second {@link Sequence} to score
     * @param type chosen type from list of pairwise sequence scoring routines
     * @param gapPenalty the gap penalties used during alignment
     * @param subMatrix the set of substitution scores used during alignment
     * @return sequence pair score
     */
    public static <S extends Sequence<C>, C extends Compound> int getPairwiseScore(
            S query, S target, PairwiseScorer type, GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix) {
        return getScorer(query, target, type, gapPenalty, subMatrix).getScore();
    }

    public static <S extends Sequence<C>, C extends Compound> Profile<S, C> getMultipleSequenceAlignment(
            List<S> sequences, MSAEmulation type, Object... settings) {
        // TODO multiple sequence alignments, convert other factories to this parameter style?
        return null;
    }

    /**
     * Factory method to run a list of alignments concurrently.  This method runs the alignments in parallel by
     * submitting all of the alignment tasks to the shared thread pool of the {@link ConcurrencyTools} utility.
     *
     * @param <A> each {@link Aligner} is of type A
     * @param <S> each {@link Sequence} of an alignment pair is of type S
     * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
     * @param alignments list of alignments to run
     * @return list of {@link SequencePair} results from running alignments
     */
    public static <A extends PairwiseSequenceAligner<S, C>, S extends Sequence<C>, C extends Compound>
            List<SequencePair<S, C>> runAligners(List<A> aligners) {
        int n = 1, all = aligners.size();
        List<Future<SequencePair<S, C>>> futures = new ArrayList<Future<SequencePair<S, C>>>();
        for (PairwiseSequenceAligner<S, C> aligner : aligners) {
            futures.add(ConcurrencyTools.submit(new CallablePairwiseSequenceAligner<S, C>(aligner),
                    String.format("Aligning pair %d of %d", n++, all)));
        }
        return getListFromFutures(futures);
    }

    /**
     * Factory method to run a list of scorers concurrently.  This method runs the scorers in parallel by submitting
     * all of the scoring tasks to the shared thread pool of the {@link ConcurrencyTools} utility.
     *
     * @param <P> each {@link Scorer} is of type P
     * @param <S> each {@link Sequence} of an alignment pair is of type S
     * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
     * @param scorers list of scorers to run
     * @return list of score results from running scorers
     */
    public static <P extends PairwiseSequenceScorer<S, C>, S extends Sequence<C>, C extends Compound>
            int[] runScorers(List<P> scorers) {
        int n = 1, all = scorers.size();
        List<Future<Integer>> futures = new ArrayList<Future<Integer>>();
        for (PairwiseSequenceScorer<S, C> scorer : scorers) {
            futures.add(ConcurrencyTools.submit(new CallablePairwiseSequenceScorer<S, C>(scorer),
                    String.format("Scoring pair %d of %d", n++, all)));
        }
        List<Integer> results = getListFromFutures(futures);
        int[] scores = new int[results.size()];
        for (int i = 0; i < scores.length; i++) {
            scores[i] = results.get(i);
        }
        return scores;
    }

    // helper methods

    // constructs a pairwise sequence alignment
    private static <S extends Sequence<C>, C extends Compound> PairwiseSequenceAligner<S, C> getAligner(
            S query, S target, PairwiseAligner type, GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix) {
        switch (type) {
        default:
        case GLOBAL:
            return new NeedlemanWunsch<S, C>(query, target, gapPenalty, subMatrix);
        case LOCAL:
            // TODO local alignment option
            return null;
        }
    }

    // constructs a list of all pairwise sequence alignments from a list of sequences
    private static <S extends Sequence<C>, C extends Compound> List<PairwiseSequenceAligner<S, C>> getAlignerList(
            List<S> sequences, PairwiseAligner type, GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix) {
        List<PairwiseSequenceAligner<S, C>> allPairs = new ArrayList<PairwiseSequenceAligner<S, C>>();
        for (int i = 0; i < sequences.size(); i++) {
            for (int j = i+1; j < sequences.size(); j++) {
                allPairs.add(getAligner(sequences.get(i), sequences.get(j), type, gapPenalty, subMatrix));
            }
        }
        return allPairs;
    }

    // retrieves calculated elements from a list on the concurrent execution queue
    private static <E> List<E> getListFromFutures(List<Future<E>> futures) {
        List<E> list = new ArrayList<E>();
        for (Future<E> f : futures) {
            // TODO when added to ConcurrencyTools, log completions and exceptions instead of printing stack traces
            try {
                list.add(f.get());
            } catch (InterruptedException e) {
                e.printStackTrace();
            } catch (ExecutionException e) {
                e.printStackTrace();
            }
        }
        return list;
    }

    // constructs a pairwise sequence scorer
    private static <S extends Sequence<C>, C extends Compound> PairwiseSequenceScorer<S, C> getScorer(
            S query, S target, PairwiseScorer type, GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix) {
        switch(type) {
        default:
        case GLOBAL:
            return getAligner(query, target, PairwiseAligner.GLOBAL, gapPenalty, subMatrix);
        case GLOBAL_IDENTITIES:
        case GLOBAL_SIMILARITIES:
        case LOCAL:
        case LOCAL_IDENTITIES:
        case LOCAL_SIMILARITIES:
        case KMERS:
            // TODO other scoring options
            return null;
        }
    }

    // constructs a list of all pairwise sequence scorers from a list of sequences
    private static <S extends Sequence<C>, C extends Compound> List<PairwiseSequenceScorer<S, C>> getScorerList(
            List<S> sequences, PairwiseScorer type, GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix) {
        List<PairwiseSequenceScorer<S, C>> allPairs = new ArrayList<PairwiseSequenceScorer<S, C>>();
        for (int i = 0; i < sequences.size(); i++) {
            for (int j = i+1; j < sequences.size(); j++) {
                allPairs.add(getScorer(sequences.get(i), sequences.get(j), type, gapPenalty, subMatrix));
            }
        }
        return allPairs;
    }

}

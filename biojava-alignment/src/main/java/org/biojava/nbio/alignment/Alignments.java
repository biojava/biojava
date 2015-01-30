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

package org.biojava.nbio.alignment;

import org.biojava.nbio.alignment.template.*;
import org.biojava.nbio.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.compound.DNACompoundSet;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.CompoundSet;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.biojava.nbio.core.util.ConcurrencyTools;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

/**
 * Static utility to easily run alignment routines.  To exit cleanly after running any parallel method that mentions
 * use of the {@link ConcurrencyTools} utility, {@link ConcurrencyTools#shutdown()} or
 * {@link ConcurrencyTools#shutdownAndAwaitTermination()} must be called.
 *
 * @author Mark Chapman
 */
public class Alignments {

	private final static Logger logger = LoggerFactory.getLogger(Alignments.class);

    /**
     * List of implemented sequence pair in a profile scoring routines.
     */
    public static enum PairInProfileScorerType {
        IDENTITIES,  // similar to MUSCLE
        SIMILARITIES
    }

    /**
     * List of implemented pairwise sequence alignment routines.
     */
    public static enum PairwiseSequenceAlignerType {
        GLOBAL,              // Needleman-Wunsch/Gotoh
        GLOBAL_LINEAR_SPACE, // Guan-Uberbacher
        LOCAL,               // Smith-Waterman/Gotoh
        LOCAL_LINEAR_SPACE   // Smith-Waterman/Gotoh with smart traceback at each maximum
    }

    /**
     * List of implemented pairwise sequence scoring routines.
     */
    public static enum PairwiseSequenceScorerType {
        GLOBAL,
        GLOBAL_IDENTITIES,   // similar to CLUSTALW and CLUSTALW2
        GLOBAL_SIMILARITIES,
        LOCAL,
        LOCAL_IDENTITIES,
        LOCAL_SIMILARITIES,
        KMERS,               // similar to CLUSTAL and MUSCLE
        WU_MANBER            // similar to KALIGN
    }

    /**
     * List of implemented profile-profile alignment routines.
     */
    public static enum ProfileProfileAlignerType {
        GLOBAL,              // similar to MUSCLE and KALIGN
        GLOBAL_LINEAR_SPACE, // similar to CLUSTALW and CLUSTALW2
        GLOBAL_CONSENSUS,    // similar to CLUSTAL
        LOCAL,
        LOCAL_LINEAR_SPACE,
        LOCAL_CONSENSUS
    }

    /**
     * List of implemented profile refinement routines.
     */
    public static enum RefinerType {
        PARTITION_SINGLE,     // similar to CLUSTALW2
        PARTITION_SINGLE_ALL, // similar to CLUSTALW2
        PARTITION_TREE,       // similar to MUSCLE
        PARTITION_TREE_ALL,
        RESCORE_IDENTITIES,   // similar to MUSCLE
        RESCORE_SIMILARITIES
    }

    // prevents instantiation
    private Alignments() { }

    // public factory methods

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
            List<S> sequences, PairwiseSequenceAlignerType type, GapPenalty gapPenalty,
            SubstitutionMatrix<C> subMatrix) {
        return runPairwiseAligners(getAllPairsAligners(sequences, type, gapPenalty, subMatrix));
    }

    /**
     * Factory method which computes a multiple sequence alignment for the given {@link List} of {@link Sequence}s.
     *
     * @param <S> each {@link Sequence} of the {@link List} is of type S
     * @param <C> each element of a {@link Sequence} is a {@link Compound} of type C
     * @param sequences the {@link List} of {@link Sequence}s to align
     * @param settings optional settings that adjust the alignment
     * @return multiple sequence alignment {@link Profile}
     */
    public static <S extends Sequence<C>, C extends Compound> Profile<S, C> getMultipleSequenceAlignment(
            List<S> sequences, Object... settings) { // TODO convert other factories to this parameter style?
        CompoundSet<C> cs = sequences.get(0).getCompoundSet();
        PairwiseSequenceScorerType ps = PairwiseSequenceScorerType.GLOBAL_IDENTITIES;
        GapPenalty gapPenalty = new SimpleGapPenalty();
        SubstitutionMatrix<C> subMatrix = null;
        if (cs == AminoAcidCompoundSet.getAminoAcidCompoundSet()) {
            @SuppressWarnings("unchecked") // compound types must be equal since compound sets are equal
            SubstitutionMatrix<C> temp = (SubstitutionMatrix<C>) SubstitutionMatrixHelper.getBlosum62();
            subMatrix = temp;
        } else if (cs == DNACompoundSet.getDNACompoundSet()) {
            @SuppressWarnings("unchecked") // compound types must be equal since compound sets are equal
            SubstitutionMatrix<C> temp = (SubstitutionMatrix<C>) SubstitutionMatrixHelper.getNuc4_4();
            subMatrix = temp;
            
        } else if (cs == AmbiguityDNACompoundSet.getDNACompoundSet()) {
            @SuppressWarnings("unchecked") // compound types must be equal since compound sets are equal
            SubstitutionMatrix<C> temp = (SubstitutionMatrix<C>) SubstitutionMatrixHelper.getNuc4_4();
            subMatrix = temp;
            
        }
        ProfileProfileAlignerType pa = ProfileProfileAlignerType.GLOBAL;
        for (Object o : settings) {
            if (o instanceof PairwiseSequenceScorerType) {
                ps = (PairwiseSequenceScorerType) o;
            } else if (o instanceof GapPenalty) {
                gapPenalty = (GapPenalty) o;
            } else if (o instanceof SubstitutionMatrix<?>) {
                if (cs != ((SubstitutionMatrix<?>) o).getCompoundSet()) {
                    throw new IllegalArgumentException(
                            "Compound sets of the sequences and substitution matrix must match.");
                }
                @SuppressWarnings("unchecked") // compound types must be equal since compound sets are equal
                SubstitutionMatrix<C> temp = (SubstitutionMatrix<C>) o;
                subMatrix = temp;
            } else if (o instanceof ProfileProfileAlignerType) {
                pa = (ProfileProfileAlignerType) o;
            }
        }

        // stage 1: pairwise similarity calculation
        List<PairwiseSequenceScorer<S, C>> scorers = getAllPairsScorers(sequences, ps, gapPenalty, subMatrix);
        runPairwiseScorers(scorers);

        // stage 2: hierarchical clustering into a guide tree
        GuideTree<S, C> tree = new GuideTree<S, C>(sequences, scorers);
        scorers = null;

        // stage 3: progressive alignment
        Profile<S, C> msa = getProgressiveAlignment(tree, pa, gapPenalty, subMatrix);

        // TODO stage 4: refinement
        return msa;
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
            S query, S target, PairwiseSequenceAlignerType type, GapPenalty gapPenalty,
            SubstitutionMatrix<C> subMatrix) {
        return getPairwiseAligner(query, target, type, gapPenalty, subMatrix).getPair();
    }

    // default access (package private) factory methods

    /**
     * Factory method which sets up a sequence alignment for all {@link Sequence} pairs in the given {@link List}.
     *
     * @param <S> each {@link Sequence} of an alignment pair is of type S
     * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
     * @param sequences the {@link List} of {@link Sequence}s to align
     * @param type chosen type from list of pairwise sequence alignment routines
     * @param gapPenalty the gap penalties used during alignment
     * @param subMatrix the set of substitution scores used during alignment
     * @return list of pairwise sequence aligners
     */
    static <S extends Sequence<C>, C extends Compound> List<PairwiseSequenceAligner<S, C>> getAllPairsAligners(
            List<S> sequences, PairwiseSequenceAlignerType type, GapPenalty gapPenalty,
            SubstitutionMatrix<C> subMatrix) {
        List<PairwiseSequenceAligner<S, C>> allPairs = new ArrayList<PairwiseSequenceAligner<S, C>>();
        for (int i = 0; i < sequences.size(); i++) {
            for (int j = i+1; j < sequences.size(); j++) {
                allPairs.add(getPairwiseAligner(sequences.get(i), sequences.get(j), type, gapPenalty, subMatrix));
            }
        }
        return allPairs;
    }

    /**
     * Factory method which sets up a sequence pair scorer for all {@link Sequence} pairs in the given {@link List}.
     *
     * @param <S> each {@link Sequence} of a pair is of type S
     * @param <C> each element of a {@link Sequence} is a {@link Compound} of type C
     * @param sequences the {@link List} of {@link Sequence}s to align
     * @param type chosen type from list of pairwise sequence scoring routines
     * @param gapPenalty the gap penalties used during alignment
     * @param subMatrix the set of substitution scores used during alignment
     * @return list of sequence pair scorers
     */
    public static <S extends Sequence<C>, C extends Compound> List<PairwiseSequenceScorer<S, C>> getAllPairsScorers(
            List<S> sequences, PairwiseSequenceScorerType type, GapPenalty gapPenalty,
            SubstitutionMatrix<C> subMatrix) {
        List<PairwiseSequenceScorer<S, C>> allPairs = new ArrayList<PairwiseSequenceScorer<S, C>>();
        for (int i = 0; i < sequences.size(); i++) {
            for (int j = i+1; j < sequences.size(); j++) {
                allPairs.add(getPairwiseScorer(sequences.get(i), sequences.get(j), type, gapPenalty, subMatrix));
            }
        }
        return allPairs;
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
    public static <S extends Sequence<C>, C extends Compound> double[] getAllPairsScores( List<S> sequences,
            PairwiseSequenceScorerType type, GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix) {
        return runPairwiseScorers(getAllPairsScorers(sequences, type, gapPenalty, subMatrix));
    }

    /**
     * Factory method which retrieves calculated elements from a list of tasks on the concurrent execution queue.
     *
     * @param <E> each task calculates a value of type E
     * @param futures list of tasks
     * @return calculated elements
     */
    static <E> List<E> getListFromFutures(List<Future<E>> futures) {
        List<E> list = new ArrayList<E>();
        for (Future<E> f : futures) {
            // TODO when added to ConcurrencyTools, log completions and exceptions instead of printing stack traces
            try {
                list.add(f.get());
            } catch (InterruptedException e) {
                logger.error("Interrupted Exception: ", e);
            } catch (ExecutionException e) {
                logger.error("Execution Exception: ", e);
            }
        }
        return list;
    }

    /**
     * Factory method which constructs a pairwise sequence aligner.
     *
     * @param <S> each {@link Sequence} of an alignment pair is of type S
     * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
     * @param query the first {@link Sequence} to align
     * @param target the second {@link Sequence} to align
     * @param type chosen type from list of pairwise sequence alignment routines
     * @param gapPenalty the gap penalties used during alignment
     * @param subMatrix the set of substitution scores used during alignment
     * @return pairwise sequence aligner
     */
    public static <S extends Sequence<C>, C extends Compound> PairwiseSequenceAligner<S, C> getPairwiseAligner(
            S query, S target, PairwiseSequenceAlignerType type, GapPenalty gapPenalty,
            SubstitutionMatrix<C> subMatrix) {
    	if (!query.getCompoundSet().equals(target.getCompoundSet())) {
    		System.err.println(query.getCompoundSet().getClass().getName() + " != " + target.getCompoundSet().getClass().getName());
    		throw new IllegalArgumentException("Sequence compound sets must be the same");
    	}
        switch (type) {
        default:
        case GLOBAL:
            return new NeedlemanWunsch<S, C>(query, target, gapPenalty, subMatrix);
        case LOCAL:
            return new SmithWaterman<S, C>(query, target, gapPenalty, subMatrix);
        case GLOBAL_LINEAR_SPACE:
        case LOCAL_LINEAR_SPACE:
            // TODO other alignment options (Myers-Miller, Thompson)
            throw new UnsupportedOperationException(Alignments.class.getSimpleName() + " does not yet support " +
                    type + " alignment");
        }
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
    static <S extends Sequence<C>, C extends Compound> double getPairwiseScore(S query, S target,
            PairwiseSequenceScorerType type, GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix) {
        return getPairwiseScorer(query, target, type, gapPenalty, subMatrix).getScore();
    }

    /**
     * Factory method which constructs a pairwise sequence scorer.
     *
     * @param <S> each {@link Sequence} of a pair is of type S
     * @param <C> each element of a {@link Sequence} is a {@link Compound} of type C
     * @param query the first {@link Sequence} to score
     * @param target the second {@link Sequence} to score
     * @param type chosen type from list of pairwise sequence scoring routines
     * @param gapPenalty the gap penalties used during alignment
     * @param subMatrix the set of substitution scores used during alignment
     * @return sequence pair scorer
     */
    static <S extends Sequence<C>, C extends Compound> PairwiseSequenceScorer<S, C> getPairwiseScorer(
            S query, S target, PairwiseSequenceScorerType type, GapPenalty gapPenalty,
            SubstitutionMatrix<C> subMatrix) {
        switch (type) {
        default:
        case GLOBAL:
            return getPairwiseAligner(query, target, PairwiseSequenceAlignerType.GLOBAL, gapPenalty, subMatrix);
        case GLOBAL_IDENTITIES:
            return new FractionalIdentityScorer<S, C>(getPairwiseAligner(query, target,
                    PairwiseSequenceAlignerType.GLOBAL, gapPenalty, subMatrix));
        case GLOBAL_SIMILARITIES:
            return new FractionalSimilarityScorer<S, C>(getPairwiseAligner(query, target,
                    PairwiseSequenceAlignerType.GLOBAL, gapPenalty, subMatrix));
        case LOCAL:
            return getPairwiseAligner(query, target, PairwiseSequenceAlignerType.LOCAL, gapPenalty, subMatrix);
        case LOCAL_IDENTITIES:
            return new FractionalIdentityScorer<S, C>(getPairwiseAligner(query, target,
                    PairwiseSequenceAlignerType.LOCAL, gapPenalty, subMatrix));
        case LOCAL_SIMILARITIES:
            return new FractionalSimilarityScorer<S, C>(getPairwiseAligner(query, target,
                    PairwiseSequenceAlignerType.LOCAL, gapPenalty, subMatrix));
        case KMERS:
        case WU_MANBER:
            // TODO other scoring options
            throw new UnsupportedOperationException(Alignments.class.getSimpleName() + " does not yet support " +
                    type + " scoring");
        }
    }

    /**
     * Factory method which constructs a profile-profile aligner.
     *
     * @param <S> each {@link Sequence} of an alignment profile is of type S
     * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
     * @param profile1 the first {@link Profile} to align
     * @param profile2 the second {@link Profile} to align
     * @param type chosen type from list of profile-profile alignment routines
     * @param gapPenalty the gap penalties used during alignment
     * @param subMatrix the set of substitution scores used during alignment
     * @return profile-profile aligner
     */
    static <S extends Sequence<C>, C extends Compound> ProfileProfileAligner<S, C> getProfileProfileAligner(
            Profile<S, C> profile1, Profile<S, C> profile2, ProfileProfileAlignerType type, GapPenalty gapPenalty,
            SubstitutionMatrix<C> subMatrix) {
        switch (type) {
        default:
        case GLOBAL:
            return new SimpleProfileProfileAligner<S, C>(profile1, profile2, gapPenalty, subMatrix);
        case GLOBAL_LINEAR_SPACE:
        case GLOBAL_CONSENSUS:
        case LOCAL:
        case LOCAL_LINEAR_SPACE:
        case LOCAL_CONSENSUS:
            // TODO other alignment options (Myers-Miller, consensus, local)
            throw new UnsupportedOperationException(Alignments.class.getSimpleName() + " does not yet support " +
                    type + " alignment");
        }
    }

    /**
     * Factory method which constructs a profile-profile aligner.
     *
     * @param <S> each {@link Sequence} of an alignment profile is of type S
     * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
     * @param profile1 the first {@link Profile} to align
     * @param profile2 the second {@link Profile} to align
     * @param type chosen type from list of profile-profile alignment routines
     * @param gapPenalty the gap penalties used during alignment
     * @param subMatrix the set of substitution scores used during alignment
     * @return profile-profile aligner
     */
    static <S extends Sequence<C>, C extends Compound> ProfileProfileAligner<S, C> getProfileProfileAligner(
            Future<ProfilePair<S, C>> profile1, Future<ProfilePair<S, C>> profile2, ProfileProfileAlignerType type,
            GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix) {
        switch (type) {
        default:
        case GLOBAL:
            return new SimpleProfileProfileAligner<S, C>(profile1, profile2, gapPenalty, subMatrix);
        case GLOBAL_LINEAR_SPACE:
        case GLOBAL_CONSENSUS:
        case LOCAL:
        case LOCAL_LINEAR_SPACE:
        case LOCAL_CONSENSUS:
            // TODO other alignment options (Myers-Miller, consensus, local)
            throw new UnsupportedOperationException(Alignments.class.getSimpleName() + " does not yet support " +
                    type + " alignment");
        }
    }

    /**
     * Factory method which constructs a profile-profile aligner.
     *
     * @param <S> each {@link Sequence} of an alignment profile is of type S
     * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
     * @param profile1 the first {@link Profile} to align
     * @param profile2 the second {@link Profile} to align
     * @param type chosen type from list of profile-profile alignment routines
     * @param gapPenalty the gap penalties used during alignment
     * @param subMatrix the set of substitution scores used during alignment
     * @return profile-profile aligner
     */
    static <S extends Sequence<C>, C extends Compound> ProfileProfileAligner<S, C> getProfileProfileAligner(
            Profile<S, C> profile1, Future<ProfilePair<S, C>> profile2, ProfileProfileAlignerType type,
            GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix) {
        switch (type) {
        default:
        case GLOBAL:
            return new SimpleProfileProfileAligner<S, C>(profile1, profile2, gapPenalty, subMatrix);
        case GLOBAL_LINEAR_SPACE:
        case GLOBAL_CONSENSUS:
        case LOCAL:
        case LOCAL_LINEAR_SPACE:
        case LOCAL_CONSENSUS:
            // TODO other alignment options (Myers-Miller, consensus, local)
            throw new UnsupportedOperationException(Alignments.class.getSimpleName() + " does not yet support " +
                    type + " alignment");
        }
    }

    /**
     * Factory method which constructs a profile-profile aligner.
     *
     * @param <S> each {@link Sequence} of an alignment profile is of type S
     * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
     * @param profile1 the first {@link Profile} to align
     * @param profile2 the second {@link Profile} to align
     * @param type chosen type from list of profile-profile alignment routines
     * @param gapPenalty the gap penalties used during alignment
     * @param subMatrix the set of substitution scores used during alignment
     * @return profile-profile aligner
     */
    static <S extends Sequence<C>, C extends Compound> ProfileProfileAligner<S, C> getProfileProfileAligner(
            Future<ProfilePair<S, C>> profile1, Profile<S, C> profile2, ProfileProfileAlignerType type,
            GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix) {
        switch (type) {
        default:
        case GLOBAL:
            return new SimpleProfileProfileAligner<S, C>(profile1, profile2, gapPenalty, subMatrix);
        case GLOBAL_LINEAR_SPACE:
        case GLOBAL_CONSENSUS:
        case LOCAL:
        case LOCAL_LINEAR_SPACE:
        case LOCAL_CONSENSUS:
            // TODO other alignment options (Myers-Miller, consensus, local)
            throw new UnsupportedOperationException(Alignments.class.getSimpleName() + " does not yet support " +
                    type + " alignment");
        }
    }

    /**
     * Factory method which computes a profile alignment for the given {@link Profile} pair.
     *
     * @param <S> each {@link Sequence} of the {@link Profile} pair is of type S
     * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
     * @param profile1 the first {@link Profile} to align
     * @param profile2 the second {@link Profile} to align
     * @param type chosen type from list of profile-profile alignment routines
     * @param gapPenalty the gap penalties used during alignment
     * @param subMatrix the set of substitution scores used during alignment
     * @return alignment profile
     */
    static <S extends Sequence<C>, C extends Compound> ProfilePair<S, C> getProfileProfileAlignment(
            Profile<S, C> profile1, Profile<S, C> profile2, ProfileProfileAlignerType type, GapPenalty gapPenalty,
            SubstitutionMatrix<C> subMatrix) {
        return getProfileProfileAligner(profile1, profile2, type, gapPenalty, subMatrix).getPair();
    }

    /**
     * Factory method to run the profile-profile alignments of a progressive multiple sequence alignment concurrently.
     * This method runs the alignments in parallel by submitting all of the alignment tasks to the shared thread pool
     * of the {@link ConcurrencyTools} utility.
     *
     * @param <S> each {@link Sequence} of the {@link Profile} pair is of type S
     * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
     * @param tree guide tree to follow aligning profiles from leaves to root
     * @param type chosen type from list of profile-profile alignment routines
     * @param gapPenalty the gap penalties used during alignment
     * @param subMatrix the set of substitution scores used during alignment
     * @return multiple sequence alignment
     */
    public static <S extends Sequence<C>, C extends Compound> Profile<S, C> getProgressiveAlignment(GuideTree<S, C> tree,
            ProfileProfileAlignerType type, GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix) {

        // find inner nodes in post-order traversal of tree (each leaf node has a single sequence profile)
        List<GuideTreeNode<S, C>> innerNodes = new ArrayList<GuideTreeNode<S, C>>();
        for (GuideTreeNode<S, C> n : tree) {
            if (n.getProfile() == null) {
                innerNodes.add(n);
            }
        }

        // submit alignment tasks to the shared thread pool
        int i = 1, all = innerNodes.size();
        for (GuideTreeNode<S, C> n : innerNodes) {
            Profile<S, C> p1 = n.getChild1().getProfile(), p2 = n.getChild2().getProfile();
            Future<ProfilePair<S, C>> pf1 = n.getChild1().getProfileFuture(), pf2 = n.getChild2().getProfileFuture();
            ProfileProfileAligner<S, C> aligner =
                    (p1 != null) ? ((p2 != null) ? getProfileProfileAligner(p1, p2, type, gapPenalty, subMatrix) :
                            getProfileProfileAligner(p1, pf2, type, gapPenalty, subMatrix)) :
                    ((p2 != null) ? getProfileProfileAligner(pf1, p2, type, gapPenalty, subMatrix) :
                            getProfileProfileAligner(pf1, pf2, type, gapPenalty, subMatrix));
            n.setProfileFuture(ConcurrencyTools.submit(new CallableProfileProfileAligner<S, C>(aligner), String.format(
                    "Aligning pair %d of %d", i++, all)));
        }

        // retrieve the alignment results
        for (GuideTreeNode<S, C> n : innerNodes) {
            // TODO when added to ConcurrencyTools, log completions and exceptions instead of printing stack traces
            try {
                n.setProfile(n.getProfileFuture().get());
            } catch (InterruptedException e) {
                logger.error("Interrupted Exception: ", e);
            } catch (ExecutionException e) {
                logger.error("Execution Exception: ", e);
            }
        }

        // the alignment profile at the root of the tree is the full multiple sequence alignment
        return tree.getRoot().getProfile();
    }

    /**
     * Factory method to run a list of alignments concurrently.  This method runs the alignments in parallel by
     * submitting all of the alignment tasks to the shared thread pool of the {@link ConcurrencyTools} utility.
     *
     * @param <S> each {@link Sequence} of an alignment pair is of type S
     * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
     * @param aligners list of alignments to run
     * @return list of {@link SequencePair} results from running alignments
     */
    static <S extends Sequence<C>, C extends Compound> List<SequencePair<S, C>>
            runPairwiseAligners(List<PairwiseSequenceAligner<S, C>> aligners) {
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
     * @param <S> each {@link Sequence} of an alignment pair is of type S
     * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
     * @param scorers list of scorers to run
     * @return list of score results from running scorers
     */
    public static <S extends Sequence<C>, C extends Compound> double[] runPairwiseScorers(
            List<PairwiseSequenceScorer<S, C>> scorers) {
        int n = 1, all = scorers.size();
        List<Future<Double>> futures = new ArrayList<Future<Double>>();
        for (PairwiseSequenceScorer<S, C> scorer : scorers) {
            futures.add(ConcurrencyTools.submit(new CallablePairwiseSequenceScorer<S, C>(scorer),
                    String.format("Scoring pair %d of %d", n++, all)));
        }
        List<Double> results = getListFromFutures(futures);
        double[] scores = new double[results.size()];
        for (int i = 0; i < scores.length; i++) {
            scores[i] = results.get(i);
        }
        return scores;
    }

    /**
     * Factory method to run a list of alignments concurrently.  This method runs the alignments in parallel by
     * submitting all of the alignment tasks to the shared thread pool of the {@link ConcurrencyTools} utility.
     *
     * @param <S> each {@link Sequence} of the {@link Profile} pair is of type S
     * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
     * @param aligners list of alignments to run
     * @return list of {@link ProfilePair} results from running alignments
     */
    static <S extends Sequence<C>, C extends Compound> List<ProfilePair<S, C>>
            runProfileAligners(List<ProfileProfileAligner<S, C>> aligners) {
        int n = 1, all = aligners.size();
        List<Future<ProfilePair<S, C>>> futures = new ArrayList<Future<ProfilePair<S, C>>>();
        for (ProfileProfileAligner<S, C> aligner : aligners) {
            futures.add(ConcurrencyTools.submit(new CallableProfileProfileAligner<S, C>(aligner),
                    String.format("Aligning pair %d of %d", n++, all)));
        }
        return getListFromFutures(futures);
    }

}

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

import java.util.List;

import org.biojava3.alignment.template.GapPenalty;
import org.biojava3.alignment.template.HierarchicalClusterer;
import org.biojava3.alignment.template.PairwiseSequenceAligner;
import org.biojava3.alignment.template.PairwiseSequenceScorer;
import org.biojava3.alignment.template.PartitionRefiner;
import org.biojava3.alignment.template.Profile;
import org.biojava3.alignment.template.ProfileProfileAligner;
import org.biojava3.alignment.template.RescoreRefiner;
import org.biojava3.alignment.template.SequencePair;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.Sequence;

public class Alignments {

    private Alignments() { } // prevents instantiation

    public static <S extends Sequence<C>, C extends Compound> List<SequencePair<S, C>> getAllPairsAlignments(List<S> sequences) {
        // TODO
        return null;
    }

    public static <S extends Sequence<C>, C extends Compound> int[] getAllPairsScores(List<S> sequences) {
        // TODO
        return null;
    }

    public static <S extends Sequence<C>, C extends Compound> SequencePair<S, C> getPairwiseAlignment(S sequence1, S sequence2) {
        // TODO
        return null;
    }

    public static <S extends Sequence<C>, C extends Compound> int getPairwiseScore(S sequence1, S sequence2) {
        // TODO
        return 0;
    }

    public static <S extends Sequence<C>, C extends Compound> Profile<S, C> getMultipleSequenceAlignment(List<S> sequences) {
        // TODO
        return null;
    }

    public static enum MSAEmulation { CLUSTALW, MUSCLE, KALIGN, CUSTOM }

    public static class Defaults { // static inner class

        public static MSAEmulation getEmulation() {
            // TODO
            return null;
        }

        public static GapPenalty getGapPenalty() {
            // TODO
            return null;
        }

        public static Class<? extends HierarchicalClusterer> getHierarchicalClusterer() {
            // TODO
            return null;
        }

        public static <S extends Sequence<C>, C extends Compound> Class<? extends PairwiseSequenceAligner<S, C>> getPairwiseSequenceAligner() {
            // TODO
            return null;
        }

        public static <S extends Sequence<C>, C extends Compound> Class<? extends PairwiseSequenceScorer<S, C>> getPairwiseSequenceScorer() {
            // TODO
            return null;
        }

        public static <S extends Sequence<C>, C extends Compound> Class<? extends PartitionRefiner<S, C>> getPartitionRefiner() {
            // TODO
            return null;
        }

        public static <S extends Sequence<C>, C extends Compound>Class<? extends ProfileProfileAligner<S, C>> getProfileProfileAligner() {
            // TODO
            return null;
        }

        public static <S extends Sequence<C>, C extends Compound> Class<? extends RescoreRefiner<S, C>> getRescoreRefiner() {
            // TODO
            return null;
        }

        public static <S extends CompoundSet<C>, C extends Compound> SubstitutionMatrix<S, C> getSubstitutionMatrix() {
            // TODO
            return null;
        }

        public static void setEmulation(MSAEmulation emulation) {
            // TODO
        }

        public static void setGapPenalty(GapPenalty gapPenalty) {
            // TODO
        }

        public static void setHierarchicalClusterer(Class<? extends HierarchicalClusterer> clusterer) {
            // TODO
        }

        public static <S extends Sequence<C>, C extends Compound> void setPairwiseSequenceAligner(Class<? extends PairwiseSequenceAligner<S, C>> aligner) {
            // TODO
        }

        public static <S extends Sequence<C>, C extends Compound> void setPairwiseSequenceScorer(Class<? extends PairwiseSequenceScorer<S, C>> scorer) {
            // TODO
        }

        public static <S extends Sequence<C>, C extends Compound> void setPartitionRefiner(Class<? extends PartitionRefiner<S, C>> refiner) {
            // TODO
        }

        public static <S extends Sequence<C>, C extends Compound> void setProfileProfileAligner(Class<? extends ProfileProfileAligner<S, C>> aligner) {
            // TODO
        }

        public static <S extends Sequence<C>, C extends Compound> void setRescoreRefiner(Class<? extends RescoreRefiner<S, C>> refiner) {
            // TODO
        }

        public static <S extends CompoundSet<C>, C extends Compound> void setSubstitutionMatrix(SubstitutionMatrix<S, C> matrix) {
            // TODO
        }

    }

}

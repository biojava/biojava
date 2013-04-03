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
 * Created on August 11, 2010
 * Author: Mark Chapman
 */

package org.biojava3.alignment.routines;

import java.util.Arrays;
import java.util.List;

import org.biojava3.alignment.SimpleSequencePair;
import org.biojava3.alignment.template.AbstractPairwiseSequenceAligner;
import org.biojava3.alignment.template.AlignedSequence;
import org.biojava3.alignment.template.AlignedSequence.Step;
import org.biojava3.alignment.template.GapPenalty;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.Sequence;

/**
 * This algorithm uses a divide-and-conquer approach to find optimal pairwise global sequence alignments (from the
 * first until the last {@link Compound} of each {@link Sequence}) with the restriction that any alignment produced
 * will connect the query sequence to the target sequence at the <em>anchors</em>.  This class performs such global
 * sequence comparisons efficiently by dynamic programming with a space requirement reduced from quadratic (a multiple
 * of query sequence length times target sequence length) to only linear (a multiple of query sequence length).  The
 * counterpoint to this reduction in space complexity is a modest (a multiple < 2) increase in time.
 *
 * @author Mark Chapman
 * @param <S> each {@link Sequence} of the alignment pair is of type S
 * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
 */
public class AnchoredPairwiseSequenceAligner<S extends Sequence<C>, C extends Compound> extends
        AbstractPairwiseSequenceAligner<S, C> {

    private static int defaultCutsPerSection = 10;

    /**
     * Sets the default number of cuts added to each section during each pass.
     * @param defaultCutsPerSection the default number of cuts added to each section during each pass
     */
    public static void setDefaultCutsPerSection(int defaultCutsPerSection) {
        AnchoredPairwiseSequenceAligner.defaultCutsPerSection = Math.max(1, defaultCutsPerSection);
    }

    /**
     * Before running a pairwise global sequence alignment, data must be sent in via calls to
     * {@link #setQuery(Sequence)}, {@link #setTarget(Sequence)}, {@link #setGapPenalty(GapPenalty)}, and
     * {@link #setSubstitutionMatrix(SubstitutionMatrix)}.
     */
    public AnchoredPairwiseSequenceAligner() {
        setCutsPerSection(defaultCutsPerSection);
    }

    /**
     * Prepares for a pairwise global sequence alignment.
     *
     * @param query the first {@link Sequence} of the pair to align
     * @param target the second {@link Sequence} of the pair to align
     * @param gapPenalty the gap penalties used during alignment
     * @param subMatrix the set of substitution scores used during alignment
     */
    public AnchoredPairwiseSequenceAligner(S query, S target, GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix) {
        this(query, target, gapPenalty, subMatrix, defaultCutsPerSection, null);
    }

    /**
     * Prepares for a pairwise global sequence alignment.
     *
     * @param query the first {@link Sequence} of the pair to align
     * @param target the second {@link Sequence} of the pair to align
     * @param gapPenalty the gap penalties used during alignment
     * @param subMatrix the set of substitution scores used during alignment
     * @param cutsPerSection the number of cuts added to each section during each pass
     */
    public AnchoredPairwiseSequenceAligner(S query, S target, GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix,
            int cutsPerSection) {
        this(query, target, gapPenalty, subMatrix, cutsPerSection, null);
    }

    /**
     * Prepares for a pairwise global sequence alignment.
     *
     * @param query the first {@link Sequence} of the pair to align
     * @param target the second {@link Sequence} of the pair to align
     * @param gapPenalty the gap penalties used during alignment
     * @param subMatrix the set of substitution scores used during alignment
     * @param cutsPerSection the number of cuts added to each section during each pass
     * @param anchors the initial list of anchors
     */
    public AnchoredPairwiseSequenceAligner(S query, S target, GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix,
            int cutsPerSection, int[] anchors) {
        super(query, target, gapPenalty, subMatrix);
        setCutsPerSection(cutsPerSection);
        setAnchors(anchors);
    }

    /**
     * Returns the list of anchors.  The populated elements correspond to query compounds with a connection established
     * to a target compound.
     *
     * @return the list of anchors
     */
    public int[] getAnchors() {
        return Arrays.copyOf(anchors, anchors.length);
    }

    /**
     * Returns the number of cuts added to each section during each pass.
     *
     * @return the number of cuts added to each section during each pass
     */
    public int getCutsPerSection() {
        return cutsPerSection;
    }

    /**
     * Sets the starting list of anchors before running the alignment routine.
     *
     * @param anchors list of points that are tied to the given indices in the target
     */
    public void setAnchors(int[] anchors) {
        this.anchors = (anchors == null) ? null : Arrays.copyOf(anchors, anchors.length);
        reset();
    }

    /**
     * Sets the number of cuts added to each section during each pass.
     *
     * @param cutsPerSection the number of cuts added to each section during each pass
     */
    public void setCutsPerSection(int cutsPerSection) {
        this.cutsPerSection = Math.max(1, cutsPerSection);
    }

    // method for AbstractMatrixAligner

    @Override
    protected void reset() {
        super.reset();
        if (getQuery() != null && getTarget() != null) {
            resetAnchors();
        }
    }

    @Override
    protected void setProfile(List<Step> sx, List<Step> sy) {
        profile = pair = new SimpleSequencePair<S, C>(getQuery(), getTarget(), sx, sy);
    }

}

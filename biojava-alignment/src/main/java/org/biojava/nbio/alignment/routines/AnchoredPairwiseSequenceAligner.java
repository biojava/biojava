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

package org.biojava.nbio.alignment.routines;

import org.biojava.nbio.core.alignment.SimpleSequencePair;
import org.biojava.nbio.alignment.routines.AlignerHelper.Anchor;
import org.biojava.nbio.alignment.template.AbstractPairwiseSequenceAligner;
import org.biojava.nbio.core.alignment.template.AlignedSequence;
import org.biojava.nbio.core.alignment.template.AlignedSequence.Step;
import org.biojava.nbio.alignment.template.GapPenalty;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.Sequence;

import java.util.ArrayList;
import java.util.List;

/**
 * This algorithm uses a divide-and-conquer approach to find optimal pairwise global sequence alignments (from the
 * first until the last {@link Compound} of each {@link Sequence}) with the restriction that any alignment produced
 * will connect the query sequence to the target sequence at the <em>anchors</em>.  This class performs such global
 * sequence comparisons efficiently by dynamic programming with a space requirement reduced from quadratic (a multiple
 * of query sequence length times target sequence length) to only linear (a multiple of query sequence length).  The
 * counterpoint to this reduction in space complexity is a modest (a multiple < 2) increase in time.
 *
 * @author Mark Chapman
 * @author Daniel Cameron
 * @param <S> each {@link Sequence} of the alignment pair is of type S
 * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
 */
public class AnchoredPairwiseSequenceAligner<S extends Sequence<C>, C extends Compound> extends
        AbstractPairwiseSequenceAligner<S, C> {

    /**
     * Before running a pairwise global sequence alignment, data must be sent in via calls to
     * {@link #setQuery(Sequence)}, {@link #setTarget(Sequence)}, {@link #setGapPenalty(GapPenalty)}, and
     * {@link #setSubstitutionMatrix(SubstitutionMatrix)}.
     */
    public AnchoredPairwiseSequenceAligner() {
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
    public AnchoredPairwiseSequenceAligner(S query, S target, GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix) {
        this(query, target, gapPenalty, subMatrix, null);
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
    public AnchoredPairwiseSequenceAligner(S query, S target, GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix, int[] anchors) {
        super(query, target, gapPenalty, subMatrix);
        setAnchors(anchors);
    }

    /**
     * Returns the list of anchors.  The populated elements correspond to query compounds with a connection established
     * to a target compound.
     *
     * @return the list of anchors
     */
    public int[] getAnchors() {
    	int[] anchor = new int[getScoreMatrixDimensions()[0] - 1];
    	for (int i = 0; i < anchor.length; i++) {
    		anchor[i] = -1;
    	}
    	for (int i = 0; i < anchors.size(); i++) {
    		anchor[anchors.get(i).getQueryIndex()] = anchors.get(i).getTargetIndex();
    	}
    	return anchor;
    }

    /**
     * Sets the starting list of anchors before running the alignment routine.
     *
     * @param anchors list of points that are tied to the given indices in the target
     */
    public void setAnchors(int[] anchors) {
    	super.anchors = new ArrayList<Anchor>();
    	if (anchors != null) {
	    	for (int i = 0; i < anchors.length; i++) {
	    		if (anchors[i] >= 0) {
	    			addAnchor(i, anchors[i]);
	    		}
	    	}
    	}
    }
    /**
     * Adds an additional anchor to the set of anchored compounds
     * @param queryIndex 0-based index of query sequence compound
     * @param targetIndex 0-base index of target sequence compound to anchor to
     */
    public void addAnchor(int queryIndex, int targetIndex) {
    	anchors.add(new Anchor(queryIndex, targetIndex));
    }

    // method for AbstractMatrixAligner

    @Override
    protected void setProfile(List<Step> sx, List<Step> sy) {
        profile = pair = new SimpleSequencePair<S, C>(getQuery(), getTarget(), sx, sy);
    }

}

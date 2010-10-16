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

import org.biojava3.alignment.template.AlignedSequence;
import org.biojava3.alignment.template.GapPenalty;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.Sequence;

/**
 * Guan and Uberbacher defined an algorithm for pairwise global sequence alignments (from the first until the last
 * {@link Compound} of each {@link Sequence}).  This class performs such global sequence comparisons efficiently by
 * dynamic programming with a space requirement reduced from quadratic (a multiple of query sequence length times
 * target sequence length) to only linear (a multiple of query sequence length).  The counterpoint to this reduction in
 * space complexity is a modest (a multiple < 2) increase in time.
 *
 * @author Mark Chapman
 * @param <S> each {@link Sequence} of the alignment pair is of type S
 * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
 */
public class GuanUberbacher<S extends Sequence<C>, C extends Compound> extends AnchoredPairwiseSequenceAligner<S, C> {

    /**
     * Before running a pairwise global sequence alignment, data must be sent in via calls to
     * {@link #setQuery(Sequence)}, {@link #setTarget(Sequence)}, {@link #setGapPenalty(GapPenalty)}, and
     * {@link #setSubstitutionMatrix(SubstitutionMatrix)}.
     */
    public GuanUberbacher() {
    }

    /**
     * Prepares for a pairwise global sequence alignment.
     *
     * @param query the first {@link Sequence} of the pair to align
     * @param target the second {@link Sequence} of the pair to align
     * @param gapPenalty the gap penalties used during alignment
     * @param subMatrix the set of substitution scores used during alignment
     */
    public GuanUberbacher(S query, S target, GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix) {
        super(query, target, gapPenalty, subMatrix);
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
    public GuanUberbacher(S query, S target, GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix,
            int cutsPerSection) {
        super(query, target, gapPenalty, subMatrix, cutsPerSection);
    }

}

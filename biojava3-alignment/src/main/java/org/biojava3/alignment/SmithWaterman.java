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
 * Created on June 24, 2010
 * Author: Mark Chapman
 */

package org.biojava3.alignment;

import java.util.List;

import org.biojava3.alignment.template.AbstractPairwiseSequenceAligner;
import org.biojava3.alignment.template.AlignedSequence;
import org.biojava3.alignment.template.AlignedSequence.Step;
import org.biojava3.alignment.template.GapPenalty;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.Sequence;

/**
 * Smith and Waterman defined an algorithm for pairwise local sequence alignments (best match of sections from each
 * {@link Sequence}).  This class performs such local sequence comparisons efficiently by dynamic programming.
 *
 * @author Mark Chapman
 * @param <S> each {@link Sequence} of the alignment pair is of type S
 * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
 */
public class SmithWaterman<S extends Sequence<C>, C extends Compound> extends AbstractPairwiseSequenceAligner<S, C> {

    private int xMax, yMax, xStart, yStart;

    /**
     * Before running a pairwise local sequence alignment, data must be sent in via calls to
     * {@link #setQuery(Sequence)}, {@link #setTarget(Sequence)}, {@link #setGapPenalty(GapPenalty)}, and
     * {@link #setSubstitutionMatrix(SubstitutionMatrix)}.
     */
    public SmithWaterman() {
    }

    /**
     * Prepares for a pairwise local sequence alignment.
     *
     * @param query the first {@link Sequence} of the pair to align
     * @param target the second {@link Sequence} of the pair to align
     * @param gapPenalty the gap penalties used during alignment
     * @param subMatrix the set of substitution scores used during alignment
     */
    public SmithWaterman(S query, S target, GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix) {
        super(query, target, gapPenalty, subMatrix);
    }

    // helper methods

    // scores with linear gap penalty; saves memory by skipping allocation of separate matching and gap arrays
    @Override
    protected void alignScoreLinear() {
        for (int x = 1; x < scores.length; x++) {
            scores[x][0] = 0;
        }
        for (int y = 1; y < scores[0].length; y++) {
            scores[0][y] = 0;
        }
        for (int x = 1; x < scores.length; x++) {
            for (int y = 1; y < scores[0].length; y++) {
                scores[x][y] = (short) Math.max(Math.max(scores[x - 1][y] + getGapPenalty().getExtensionPenalty(),
                        scores[x][y - 1] + getGapPenalty().getExtensionPenalty()), Math.max(scores[x - 1][y - 1] +
                        alignScoreColumns(x, y), 0));
                if (scores[x][y] > scores[xMax][yMax]) {
                    xMax = x;
                    yMax = y;
                }
            }
        }
    }

    // traces back through score matrix; chooses highroad alignment
    @Override
    protected void alignTracebackLinear(List<Step> sx, List<Step> sy) {
        int x = xMax, y = yMax;
        while (scores[x][y] > 0) {
            if (scores[x][y] == scores[x - 1][y] + getGapPenalty().getExtensionPenalty()) {
                sx.add(0, Step.COMPOUND);
                sy.add(0, Step.GAP);
                x--;
            } else if (scores[x][y] == scores[x - 1][y - 1] + alignScoreColumns(x, y)) {
                sx.add(0, Step.COMPOUND);
                sy.add(0, Step.COMPOUND);
                x--;
                y--;
            } else {
                sx.add(0, Step.GAP);
                sy.add(0, Step.COMPOUND);
                y--;
            }
        }
        xStart = x;
        yStart = y;
    }

    // scores with affine gap penalty
    @Override
    protected void alignScoreAffine(short[][] ix, short[][] iy) {
        GapPenalty gapPenalty = getGapPenalty();
        short min = (short) (Short.MIN_VALUE - gapPenalty.getOpenPenalty() - gapPenalty.getExtensionPenalty());
        ix[0][0] = iy[0][0] = gapPenalty.getOpenPenalty();
        for (int x = 1; x < scores.length; x++) {
            scores[x][0] = 0;
            ix[x][0] = iy[x][0] = min;
        }
        for (int y = 1; y < scores[0].length; y++) {
            scores[0][y] = 0;
            iy[0][y] = ix[0][y] = min;
        }
        for (int x = 1; x < scores.length; x++) {
            for (int y = 1; y < scores[0].length; y++) {
                scores[x][y] = (short) Math.max(Math.max(Math.max(scores[x - 1][y - 1], ix[x - 1][y - 1]),
                        iy[x - 1][y - 1]) + alignScoreColumns(x, y), 0);
                ix[x][y] = (short) (Math.max(scores[x - 1][y] + gapPenalty.getOpenPenalty(), ix[x - 1][y]) +
                        gapPenalty.getExtensionPenalty());
                iy[x][y] = (short) (Math.max(scores[x][y - 1] + gapPenalty.getOpenPenalty(), iy[x][y - 1]) +
                        gapPenalty.getExtensionPenalty());
                if (scores[x][y] > scores[xMax][yMax]) {
                    xMax = x;
                    yMax = y;
                }
            }
        }
    }

    // traces back through score matrices; chooses highroad alignment
    @Override
    protected void alignTracebackAffine(List<Step> sx, List<Step> sy, short[][] ix, short[][] iy) {
        int x = xMax, y = yMax;
        Last last = Last.M;
        while (scores[x][y] > 0) {
            switch (last) {
            case IX:
                sx.add(0, Step.COMPOUND);
                sy.add(0, Step.GAP);
                x--;
                last = (scores[x][y] + getGapPenalty().getOpenPenalty() > ix[x][y]) ? Last.M : Last.IX;
                break;
            case M:
                sx.add(0, Step.COMPOUND);
                sy.add(0, Step.COMPOUND);
                x--;
                y--;
                int max = Math.max(Math.max(scores[x][y], ix[x][y]), iy[x][y]);
                last = (max == ix[x][y]) ? Last.IX : ((max == scores[x][y]) ? Last.M : Last.IY);
                break;
            case IY:
                sx.add(0, Step.GAP);
                sy.add(0, Step.COMPOUND);
                y--;
                last = (scores[x][y] + getGapPenalty().getOpenPenalty() >= iy[x][y]) ? Last.M : Last.IY;
            }
        }
        xStart = x;
        yStart = y;
    }

    // sets output fields
    @Override
    protected void alignSetOutputs(List<Step> sx, List<Step> sy) {
        score = scores[xMax][yMax];
        profile = pair = new SimpleSequencePair<S, C>(getQuery(), getTarget(), sx, xStart, scores.length - 1 - xMax,
                sy, yStart, scores[0].length - 1 - yMax);
    }

    // resets output fields
    @Override
    protected void reset() {
        super.reset();
        xMax = yMax = xStart = yStart = score = min = 0;
    }

}

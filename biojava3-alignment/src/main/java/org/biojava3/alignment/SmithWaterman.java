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

import java.util.ArrayList;
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

    // helper enum for alignment with affine gap penalties
    private enum Last { M, IX, IY }

    // helper method that performs alignment
    @Override
    protected void align() {
        reset();
        S query = getQuery(), target = getTarget();
        GapPenalty gapPenalty = getGapPenalty();
        SubstitutionMatrix<C> subMatrix = getSubstitutionMatrix();

        if (query == null || target == null || gapPenalty == null || subMatrix == null
                || !query.getCompoundSet().equals(target.getCompoundSet())) {
            return;
        }

        long timeStart = System.nanoTime();
        scores = new short[query.getLength() + 1][target.getLength() + 1];
        scores[0][0] = 0;
        int x, y, xMax = 0, yMax = 0;
        short[][] ix, iy;
        List<Step> sx = new ArrayList<Step>(), sy = new ArrayList<Step>();

        if (gapPenalty.getType() == GapPenalty.Type.LINEAR) {
            // scoring: saves memory by skipping allocation of separate matching and gap arrays
            for (x = 1; x < scores.length; x++) {
                scores[x][0] = 0;
            }
            for (y = 1; y < scores[0].length; y++) {
                scores[0][y] = 0;
            }
            for (x = 1; x < scores.length; x++) {
                for (y = 1; y < scores[0].length; y++) {
                    scores[x][y] = (short) Math.max(Math.max(scores[x - 1][y] + gapPenalty.getExtensionPenalty(),
                            scores[x][y - 1] + gapPenalty.getExtensionPenalty()), Math.max(scores[x - 1][y - 1] +
                            subMatrix.getValue(query.getCompoundAt(x), target.getCompoundAt(y)), 0));
                    if (scores[x][y] > scores[xMax][yMax]) {
                        xMax = x;
                        yMax = y;
                    }
                }
            }
            // traceback: chooses highroad alignment
            x = xMax;
            y = yMax;
            while (scores[x][y] > 0) {
                if (scores[x][y] == scores[x - 1][y] + gapPenalty.getExtensionPenalty()) {
                    sx.add(0, Step.COMPOUND);
                    sy.add(0, Step.GAP);
                    x--;
                } else if (scores[x][y] == scores[x - 1][y - 1] + subMatrix.getValue(query.getCompoundAt(x),
                        target.getCompoundAt(y))) {
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
        } else {
            // scoring
            ix = new short[scores.length][scores[0].length];
            iy = new short[scores.length][scores[0].length];
            short min = (short) (Short.MIN_VALUE - gapPenalty.getOpenPenalty() - gapPenalty.getExtensionPenalty());
            ix[0][0] = iy[0][0] = gapPenalty.getOpenPenalty();
            for (x = 1; x < scores.length; x++) {
                scores[x][0] = 0;
                ix[x][0] = iy[x][0] = min;
            }
            for (y = 1; y < scores[0].length; y++) {
                scores[0][y] = 0;
                iy[0][y] = ix[0][y] = min;
            }
            for (x = 1; x < scores.length; x++) {
                for (y = 1; y < scores[0].length; y++) {
                    scores[x][y] = (short) Math.max(Math.max(Math.max(scores[x - 1][y - 1], ix[x - 1][y - 1]),
                            iy[x - 1][y - 1]) + subMatrix.getValue(query.getCompoundAt(x), target.getCompoundAt(y)),
                            0);
                    ix[x][y] = (short) (Math.max(scores[x - 1][y] + gapPenalty.getOpenPenalty(), ix[x - 1][y])
                            + gapPenalty.getExtensionPenalty());
                    iy[x][y] = (short) (Math.max(scores[x][y - 1] + gapPenalty.getOpenPenalty(), iy[x][y - 1])
                            + gapPenalty.getExtensionPenalty());
                    if (scores[x][y] > scores[xMax][yMax]) {
                        xMax = x;
                        yMax = y;
                    }
                }
            }
            // traceback: chooses highroad alignment
            x = xMax;
            y = yMax;
            Last last = Last.M;
            while (scores[x][y] > 0) {
                switch (last) {
                case IX:
                    sx.add(0, Step.COMPOUND);
                    sy.add(0, Step.GAP);
                    x--;
                    last = (scores[x][y] + gapPenalty.getOpenPenalty() > ix[x][y]) ? Last.M : Last.IX;
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
                    last = (scores[x][y] + gapPenalty.getOpenPenalty() >= iy[x][y]) ? Last.M : Last.IY;
                }
            }
            // save maximum of three score matrices in scores
            for (int i = 0; i < scores.length; i++) {
                for (int j = 0; j < scores[0].length; j++) {
                    scores[i][j] = (short) Math.max(Math.max(scores[i][j], ix[i][j]), iy[i][j]);
                }
            }
        }

        // set output fields 
        score = scores[xMax][yMax];
        pair = new SimpleSequencePair<S, C>(query, target, sx, x, scores.length - 1 - xMax, sy, y,
                scores[0].length - 1 - yMax);
        time = System.nanoTime() - timeStart;

        // save memory by deleting score matrix
        if (!isStoringScoreMatrix()) {
            scores = null;
        }
    }

    // helper method that resets output fields
    @Override
    protected void reset() {
        super.reset();
        score = min = 0;
    }

}

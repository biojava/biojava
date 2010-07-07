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
 * Created on July 2, 2010
 * Author: Mark Chapman
 */

package org.biojava3.alignment.template;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.biojava3.alignment.template.AlignedSequence.Step;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.Sequence;

/**
 * Implements common code for an {@link Aligner} which builds a score matrix during computation.
 *
 * @author Mark Chapman
 * @param <S> each element of the alignment {@link Profile} is of type S
 * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
 */
public abstract class AbstractMatrixAligner<S extends Sequence<C>, C extends Compound> implements MatrixAligner<S, C> {

    // input fields
    private GapPenalty gapPenalty;
    private SubstitutionMatrix<C> subMatrix;
    private boolean storingScoreMatrix;

    // output fields
    protected short max, min, score;
    protected short[][] scores;
    protected long time = -1;
    protected Profile<S, C> profile;

    /**
     * Before running an alignment, data must be sent in via calls to {@link #setGapPenalty(GapPenalty)} and
     * {@link #setSubstitutionMatrix(SubstitutionMatrix)}.
     */
    protected AbstractMatrixAligner() {
    }

    /**
     * Prepares for an alignment.
     *
     * @param gapPenalty the gap penalties used during alignment
     * @param subMatrix the set of substitution scores used during alignment
     */
    protected AbstractMatrixAligner(GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix) {
        this.gapPenalty = gapPenalty;
        this.subMatrix = subMatrix;
        reset();
    }

    /**
     * Returns the gap penalties.
     *
     * @return the gap penalties used during alignment
     */
    public GapPenalty getGapPenalty() {
        return gapPenalty;
    }

    /**
     * Returns the substitution matrix.
     *
     * @return the set of substitution scores used during alignment
     */
    public SubstitutionMatrix<C> getSubstitutionMatrix() {
        return subMatrix;
    }

    /**
     * Returns choice to cache the score matrix or to save memory by deleting score matrix after alignment.
     *
     * @return choice to cache the score matrix
     */
    public boolean isStoringScoreMatrix() {
        return storingScoreMatrix;
    }

    /**
     * Sets the gap penalties.
     *
     * @param gapPenalty the gap penalties used during alignment
     */
    public void setGapPenalty(GapPenalty gapPenalty) {
        this.gapPenalty = gapPenalty;
        reset();
    }

    /**
     * Sets the substitution matrix.
     *
     * @param subMatrix the set of substitution scores used during alignment
     */
    public void setSubstitutionMatrix(SubstitutionMatrix<C> subMatrix) {
        this.subMatrix = subMatrix;
        reset();
    }

    /**
     * Sets choice to cache the score matrix or to save memory by deleting score matrix after alignment.
     *
     * @param storingScoreMatrix choice to cache the score matrix
     */
    public void setStoringScoreMatrix(boolean storingScoreMatrix) {
        this.storingScoreMatrix = storingScoreMatrix;
        if (!storingScoreMatrix) {
            scores = null;
        }
    }

    // methods for MatrixAligner

    @Override
    public short[][] getScoreMatrix() {
        boolean tempStoringScoreMatrix = storingScoreMatrix;
        if (scores == null) {
            storingScoreMatrix = true;
            align();
            if (scores == null) {
                return null;
            }
        }
        short[][] copy = scores;
        if (tempStoringScoreMatrix) {
            copy = new short[scores.length][scores[0].length];
            for (int i = 0; i < copy.length; i++) {
                copy[i] = Arrays.copyOf(scores[i], scores[i].length);
            }
        }
        setStoringScoreMatrix(tempStoringScoreMatrix);
        return copy;
    }

    @Override
    public short getScoreMatrixAt(int queryIndex, int targetIndex) {
        boolean tempStoringScoreMatrix = storingScoreMatrix;
        if (scores == null) {
            storingScoreMatrix = true;
            align();
            if (scores == null) {
                return Short.MIN_VALUE;
            }
        }
        short score = scores[queryIndex][targetIndex];
        setStoringScoreMatrix(tempStoringScoreMatrix);
        return score;
    }

    // methods for Aligner

    @Override
    public long getComputationTime() {
        if (profile == null) {
            align();
        }
        return time;
    }

    @Override
    public Profile<S, C> getProfile() {
        if (profile == null) {
            align();
        }
        return profile;
    }

    // methods for Scorer

    @Override
    public int getMaxScore() {
        if (profile == null) {
            align();
        }
        return max;
    }

    @Override
    public int getMinScore() {
        if (profile == null) {
            align();
        }
        return min;
    }

    @Override
    public int getScore() {
        if (profile == null) {
            align();
        }
        return score;
    }

    // performs alignment
    protected void align() {
        if (!alignReady()) {
            return;
        }

        long timeStart = System.nanoTime();

        List<Step> sx = new ArrayList<Step>(), sy = new ArrayList<Step>();
        if (getGapPenalty().getType() == GapPenalty.Type.LINEAR) {
            alignScoreLinear();
            alignTracebackLinear(sx, sy);
        } else {
            short[][] ix = new short[scores.length][scores[0].length], iy = new short[scores.length][scores[0].length];
            alignScoreAffine(ix, iy);
            alignTracebackAffine(sx, sy, ix, iy);
            // save maximum of three score matrices in scores
            if (isStoringScoreMatrix()) {
                for (int x = 0; x < scores.length; x++) {
                    for (int y = 0; y < scores[0].length; y++) {
                        scores[x][y] = (short) Math.max(Math.max(scores[x][y], ix[x][y]), iy[x][y]);
                    }
                }
            }
        }
        alignSetOutputs(sx, sy);

        time = System.nanoTime() - timeStart;

        // save memory by deleting score matrix
        if (!isStoringScoreMatrix()) {
            scores = null;
        }
    }

    // prepares for alignment; returns true if everything is set to run the alignment
    protected abstract boolean alignReady();

    // scores alignment of two columns
    protected abstract short alignScoreColumns(int queryColumn, int targetColumn);

    // scores with linear gap penalty; saves memory by skipping allocation of separate matching and gap arrays
    protected abstract void alignScoreLinear();

    // traces back through score matrix; chooses highroad alignment
    protected abstract void alignTracebackLinear(List<Step> sx, List<Step> sy);

    // scores with affine gap penalty
    protected abstract void alignScoreAffine(short[][] ix, short[][] iy);

    // traces back through score matrices; chooses highroad alignment
    protected abstract void alignTracebackAffine(List<Step> sx, List<Step> sy, short[][] ix, short[][] iy);

    // sets output fields
    protected abstract void alignSetOutputs(List<Step> sx, List<Step> sy);

    // resets output fields
    protected abstract void reset();

    // helper enum for alignment with affine gap penalties; TODO? move to subclasses
    public enum Last { M, IX, IY }

}

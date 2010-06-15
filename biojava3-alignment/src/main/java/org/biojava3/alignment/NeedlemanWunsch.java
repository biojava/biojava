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
 * Created on June 11, 2010
 * Author: Mark Chapman
 */

package org.biojava3.alignment;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.biojava3.alignment.template.AlignedSequence;
import org.biojava3.alignment.template.AlignedSequence.Step;
import org.biojava3.alignment.template.Aligner;
import org.biojava3.alignment.template.GapPenalty;
import org.biojava3.alignment.template.PairwiseSequenceAligner;
import org.biojava3.alignment.template.Profile;
import org.biojava3.alignment.template.SequencePair;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.Sequence;

/**
 * Needleman and Wunsch defined an {@link Aligner} for the problem of pairwise global sequence alignments, from the
 * first until the last {@link Compound} of each {@link Sequence}.  This class performs such global sequence
 * comparisons efficiently by dynamic programming.
 *
 * @author Mark Chapman
 * @param <S> each {@link Sequence} of the alignment pair is of type S
 * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
 */
public class NeedlemanWunsch<S extends Sequence<C>, C extends Compound> implements PairwiseSequenceAligner<S, C> {

    private static final String newLine = System.getProperty("line.separator");

    // input fields
    private S query, target;
    private GapPenalty gapPenalty;
    private SubstitutionMatrix<C> subMatrix;
    private boolean storingScoreMatrix;

    // output fields
    private SequencePair<S, C> pair;
    private short[][] scores;
    private short max, min, score;
    private long time;

    public NeedlemanWunsch() {
        this(null, null, null, null);
    }

    public NeedlemanWunsch(S query, S target) {
        this(query, target, null, null);
    }

    public NeedlemanWunsch(GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix) {
        this(null, null, gapPenalty, subMatrix);
    }

    public NeedlemanWunsch(S query, S target, GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix) {
        // TODO proper type checking
        NeedlemanWunsch<S, C> def = (NeedlemanWunsch<S, C>) Default.getInstance();
        this.query = query;
        this.target = target;
        this.gapPenalty = (gapPenalty != null) ? gapPenalty : ((def != null) ? def.getGapPenalty() : null);
        this.subMatrix = (subMatrix != null) ? subMatrix : ((def != null) ? def.getSubstitutionMatrix() : null);
    }

    public S getQuerySequence() {
        return query;
    }

    public S getTargetSequence() {
        return target;
    }

    public GapPenalty getGapPenalty() {
        return gapPenalty;
    }

    public SubstitutionMatrix<C> getSubstitutionMatrix() {
        return subMatrix;
    }

    public boolean isStoringScoreMatrix() {
        return storingScoreMatrix;
    }

    public void setQuerySequence(S query) {
        this.query = query;
        reset();
    }

    public void setTargetSequence(S target) {
        this.target = target;
        reset();
    }

    public void setGapPenalty(GapPenalty gapPenalty) {
        this.gapPenalty = gapPenalty;
        reset();
    }

    public void setSubstitutionMatrix(SubstitutionMatrix<C> subMatrix) {
        this.subMatrix = subMatrix;
        reset();
    }

    public void setStoringScoreMatrix(boolean storingScoreMatrix) {
        this.storingScoreMatrix = storingScoreMatrix;
        if (!storingScoreMatrix) {
            scores = null;
        }
    }

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
    public String getScoreMatrixAsString() {
        boolean tempStoringScoreMatrix = storingScoreMatrix;
        if (scores == null) {
            storingScoreMatrix = true;
            align();
            if (scores == null) {
                return null;
            }
        }
        StringBuilder s = new StringBuilder();
        CompoundSet<C> compoundSet = query.getCompoundSet();
        int lengthCompound = compoundSet.getMaxSingleCompoundStringLength(), lengthRest =
                Math.max(Math.max(Short.toString(min).length(), Short.toString(max).length()), lengthCompound) + 1;
        String padCompound = "%" + Integer.toString(lengthCompound) + "s",
                padRest = "%" + Integer.toString(lengthRest);
        s.append(String.format(padCompound, ""));
        s.append(String.format(padRest + "s", ""));
        for (C col : target.getAsList()) {
            s.append(String.format(padRest + "s", compoundSet.getStringForCompound(col)));
        }
        s.append(newLine);
        for (int row = 0; row <= query.getLength(); row++) {
            s.append(String.format(padCompound, (row == 0) ? "" :
                    compoundSet.getStringForCompound(query.getCompoundAt(row))));
            for (int col = 0; col <= target.getLength(); col++) {
                s.append(String.format(padRest + "d", scores[row][col]));
            }
            s.append(newLine);
        }
        setStoringScoreMatrix(tempStoringScoreMatrix);
        return s.toString();
    }

    @Override
    public short getScoreMatrixAt(int queryIndex, int targetIndex) {
        boolean tempStoringScoreMatrix = storingScoreMatrix;
        if (scores == null) {
            storingScoreMatrix = true;
            align();
        }
        short score = scores[queryIndex][targetIndex];
        setStoringScoreMatrix(tempStoringScoreMatrix);
        return score;
    }

    @Override
    public long getComputationTime() {
        if (pair == null) {
            align();
        }
        return time;
    }

    @Override
    public Profile<S, C> getProfile() {
        if (pair == null) {
            align();
        }
        return pair;
    }

    @Override
    public int getMaxScore() {
        if (pair == null) {
            align();
        }
        return max;
    }

    @Override
    public int getMinScore() {
        if (pair == null) {
            align();
        }
        return min;
    }

    @Override
    public int getScore() {
        if (pair == null) {
            align();
        }
        return score;
    }

    @Override
    public SequencePair<S, C> getPair() {
        if (pair == null) {
            align();
        }
        return pair;
    }

    // helper class for alignment with affine gap penalties
    private enum Last { M, IX, IY }

    // helper method that performs alignment
    private void align() {
        reset();
        if (query == null || target == null || gapPenalty == null || subMatrix == null
                || !query.getCompoundSet().equals(target.getCompoundSet())) {
            return;
        }

        long timeStart = System.nanoTime();
        scores = new short[query.getLength() + 1][target.getLength() + 1];
        scores[0][0] = 0;
        int x, y;
        short[][] ix, iy;
        List<Step> sx = new ArrayList<Step>(), sy = new ArrayList<Step>();

        if (gapPenalty.getType() == GapPenalty.Type.LINEAR) {
            // scoring: saves memory by skipping allocation of separate matching and gap arrays
            for (x = 1; x < scores.length; x++) {
                scores[x][0] = (short) (scores[x - 1][0] + gapPenalty.getExtensionPenalty());
            }
            for (y = 1; y < scores[0].length; y++) {
                scores[0][y] = (short) (scores[0][y - 1] + gapPenalty.getExtensionPenalty());
            }
            for (x = 1; x < scores.length; x++) {
                for (y = 1; y < scores[0].length; y++) {
                    scores[x][y] = (short) Math.max(Math.max(scores[x - 1][y] + gapPenalty.getExtensionPenalty(), scores[x][y - 1] + gapPenalty.getExtensionPenalty()),
                            scores[x - 1][y - 1] + subMatrix.getValue(query.getCompoundAt(x),
                                    target.getCompoundAt(y)));
                }
            }
            // traceback: chooses highroad alignment
            x = scores.length - 1;
            y = scores[0].length - 1;
            while (x > 0 || y > 0) {
                if (x == 0) {
                    sx.add(0, Step.GAP);
                    sy.add(0, Step.COMPOUND);
                    y--;
                } else if (y == 0 || scores[x][y] == scores[x - 1][y] + gapPenalty.getExtensionPenalty()) {
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
            ix[0][0] = iy[0][0] = gapPenalty.getOpenPenalty();
            for (x = 1; x < scores.length; x++) {
                scores[x][0] = iy[x][0] = Short.MIN_VALUE;
                ix[x][0] = (short) (ix[x - 1][0] + gapPenalty.getExtensionPenalty());
            }
            for (y = 1; y < scores[0].length; y++) {
                scores[0][y] = ix[0][y] = Short.MIN_VALUE;
                iy[0][y] = (short) (iy[0][y - 1] + gapPenalty.getExtensionPenalty());
            }
            for (x = 1; x < scores.length; x++) {
                for (y = 1; y < scores[0].length; y++) {
                    scores[x][y] = (short) (Math.max(Math.max(scores[x - 1][y - 1], ix[x - 1][y - 1]),
                            iy[x - 1][y - 1]) + subMatrix.getValue(query.getCompoundAt(x), target.getCompoundAt(y)));
                    ix[x][y] = (short) (Math.max(scores[x - 1][y] + gapPenalty.getOpenPenalty(), ix[x - 1][y])
                            + gapPenalty.getExtensionPenalty());
                    iy[x][y] = (short) (Math.max(scores[x][y - 1] + gapPenalty.getOpenPenalty(), iy[x][y - 1])
                            + gapPenalty.getExtensionPenalty());
                }
            }
            // traceback: chooses highroad alignment
            x = scores.length - 1;
            y = scores[0].length - 1;
            int max = Math.max(Math.max(scores[x][y], ix[x][y]), iy[x][y]);
            Last last = (max == ix[x][y]) ? Last.IX : ((max == scores[x][y]) ? Last.M : Last.IY);
            while (x > 0 || y > 0) {
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
                    max = Math.max(Math.max(scores[x][y], ix[x][y]), iy[x][y]);
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
            for (x = 0; x < scores.length; x++) {
                for (y = 0; y < scores[0].length; y++) {
                    scores[x][y] = (short) Math.max(Math.max(scores[x][y], ix[x][y]), iy[x][y]);
                }
            }
        }

        // set output fields 
        pair = new SimpleSequencePair<S, C>(query, target, sx, sy);
        score = scores[x][y];
        time = System.nanoTime() - timeStart;

        // save memory by deleting score matrix
        if (!storingScoreMatrix) {
            scores = null;
        }
    }

    // helper method that resets output fields
    private void reset() {
        pair = null;
        scores = null;
        int subLength = Math.min(query.getLength(), target.getLength()), maxLength = query.getLength()
                + target.getLength(), penalties = gapPenalty.getOpenPenalty() + gapPenalty.getExtensionPenalty();
        max = (short) (subLength * subMatrix.getMaxValue());
        score = min = (short) Math.min(subLength * subMatrix.getMinValue() + (maxLength - subLength) * penalties,
                maxLength * penalties);
        time = 0;
    }

    public static class Default {

        private static NeedlemanWunsch<?, ?> instance;

        public static NeedlemanWunsch<?, ?> getInstance() {
            return instance;
        }

        public static <S extends Sequence<C>, C extends Compound> void set(S query, S target) {
            instance = new NeedlemanWunsch<S, C>(query, target);
        }

        public static <S extends Sequence<C>, C extends Compound> void set(GapPenalty gapPenalty,
                SubstitutionMatrix<C> subMatrix) {
            instance = new NeedlemanWunsch<S, C>(gapPenalty, subMatrix);
        }

        public static <S extends Sequence<C>, C extends Compound> void set(S query, S target, GapPenalty gapPenalty,
                SubstitutionMatrix<C> subMatrix) {
            instance = new NeedlemanWunsch<S, C>(query, target, gapPenalty, subMatrix);
        }

    }

}

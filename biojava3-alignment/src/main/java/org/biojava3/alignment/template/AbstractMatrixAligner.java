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
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.Sequence;

/**
 * Implements common code for an {@link Aligner} which builds a score matrix during computation.
 *
 * @author Mark Chapman
 * @param <S> each element of the alignment {@link Profile} is of type S
 * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
 */
public abstract class AbstractMatrixAligner<S extends Sequence<C>, C extends Compound> extends AbstractScorer
        implements MatrixAligner<S, C> {

    // helper class for alignment with linear gap penalties
    public static class Cross {

        public int y;

        public Cross(int y) {
            this.y = y;
        }

    }

    // helper class for alignment with affine gap penalties
    public static class AffineCross extends Cross {

        public Last last;

        public AffineCross(int x, Last last) {
            super(x);
            this.last = last;
        }

    }

    // helper enum for alignment with affine gap penalties
    public enum Last {
        SUB,
        DEL,
        INS
    }

    // input fields
    private GapPenalty gapPenalty;
    private SubstitutionMatrix<C> subMatrix;
    private boolean local, storingScoreMatrix;

    // output fields
    protected Profile<S, C> profile;
    protected int xMax, yMax, xStart, yStart;
    protected short max, min, score, maxScore;
    private short[][][] scores;
    private String[] types;
    private long time = -1;

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
        this(gapPenalty, subMatrix, false);
    }

    /**
     * Prepares for an alignment.
     *
     * @param gapPenalty the gap penalties used during alignment
     * @param subMatrix the set of substitution scores used during alignment
     * @param local if true, find a region of similarity rather than aligning every compound
     */
    protected AbstractMatrixAligner(GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix, boolean local) {
        this.gapPenalty = gapPenalty;
        this.subMatrix = subMatrix;
        this.local = local;
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
     * Returns whether alignment finds a region of similarity rather than aligning every compound.
     *
     * @return true if alignment finds a region of similarity rather than aligning every compound
     */
    public boolean isLocal() {
        return local;
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
    public short[][][] getScoreMatrix() {
        boolean tempStoringScoreMatrix = storingScoreMatrix;
        if (scores == null) {
            storingScoreMatrix = true;
            align();
            if (scores == null) {
                return null;
            }
        }
        short[][][] copy = scores;
        if (tempStoringScoreMatrix) {
            copy = new short[scores.length][scores[0].length][];
            for (int i = 0; i < copy.length; i++) {
                for (int j = 0; j < copy[0].length; j++) {
                    copy[i][j] = Arrays.copyOf(scores[i][j], scores[i][j].length);
                }
            }
        }
        setStoringScoreMatrix(tempStoringScoreMatrix);
        return copy;
    }

    @Override
    public String getScoreMatrixAsString() {
        short[][][] scores = getScoreMatrix();
        StringBuilder s = new StringBuilder();
        CompoundSet<C> compoundSet = getCompoundSet();
        int lengthCompound = compoundSet.getMaxSingleCompoundStringLength(), lengthRest =
                Math.max(Math.max(Short.toString(min).length(), Short.toString(max).length()), lengthCompound) + 1;
        String padCompound = "%" + Integer.toString(lengthCompound) + "s",
                padRest = "%" + Integer.toString(lengthRest);
        List<C> query = getCompoundsOfQuery(), target = getCompoundsOfTarget();
        for (int type = 0; type < getScoreMatrixDimensions()[2]; type++) {
            if (type > 0) {
                s.append(String.format("%n"));
            }
            if (types[type] != null) {
                s.append(String.format("%s%n", types[type]));
            }
            s.append(String.format(padCompound, ""));
            s.append(String.format(padRest + "s", ""));
            for (C col : target) {
                s.append(String.format(padRest + "s", compoundSet.getStringForCompound(col)));
            }
            s.append(String.format("%n"));
            for (int x = 0; x < scores.length; x++) {
                s.append(String.format(padCompound, (x == 0) ? "" :
                    compoundSet.getStringForCompound(query.get(x - 1))));
                for (int y = 0; y < scores[0].length; y++) {
                    s.append(scores[x][y][type] >= min ? String.format(padRest + "d", scores[x][y][type]) :
                            String.format(padRest + "s", "-\u221E"));
                }
                s.append(String.format("%n"));
            }
        }
        return s.toString();
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

    // helper methods

    // performs alignment
    protected void align() {
        if (!isReady()) {
            return;
        }

        long timeStart = System.nanoTime();

        int xmax = getScoreMatrixDimensions()[0];
        scores = new short[xmax][][];
        Cross[][][] traceback = new Cross[xmax][][];
        List<Step> sx = new ArrayList<Step>(), sy = new ArrayList<Step>();
        if (getGapPenalty().getType() == GapPenalty.Type.LINEAR) {
            types = new String[] { null };
            if (local) {
                for (int x = 0; x < xmax; x++) {
                    setScoreVectorLinearLocal(x, traceback);
                }
            } else {
                for (int x = 0; x < xmax; x++) {
                    setScoreVectorLinearGlobal(x, traceback);
                }
            }
            setStepsLinear(traceback, sx, sy);
        } else {
            types = new String[] { "Substitution", "Deletion", "Insertion" };
            if (local) {
                for (int x = 0; x < xmax; x++) {
                    setScoreVectorLocal(x, traceback);
                }
                setStepsLocal(traceback, sx, sy);
            } else {
                for (int x = 0; x < xmax; x++) {
                    setScoreVectorGlobal(x, traceback);
                }
                setStepsGlobal(traceback, sx, sy);
            }
        }
        setOutputs(sx, sy);

        time = System.nanoTime() - timeStart;

        // save memory by deleting score matrix
        if (!storingScoreMatrix) {
            scores = null;
        }
    }

    // resets output fields; should be overridden to set max and min
    protected void reset() {
        xMax = yMax = xStart = yStart = 0;
        scores = null;
        types = null;
        time = -1;
        profile = null;
    }

    // sets output fields for a given alignment path; should be overridden to set profile
    protected void setOutputs(List<Step> sx, List<Step> sy) {
        if (local) {
            score = scores[xMax][yMax][0];
        } else {
            int x = scores.length - 1, y = scores[x].length - 1;
            for (int z = 0; z < scores[x][y].length; z++) {
                score = (short) Math.max(score, scores[x][y][z]);
            }
        }
    }

    // scores alignment for a given position in both sequences
    protected void setScorePoint(int x, int y, Cross[][][] traceback) {

        // substitution
        if (scores[x - 1][y - 1][1] >= scores[x - 1][y - 1][0] && scores[x - 1][y - 1][1] >= scores[x - 1][y - 1][2]) {
            scores[x][y][0] = (short) (scores[x - 1][y - 1][1] + getSubstitutionScore(x, y));
            traceback[x][y][0] = new AffineCross(y - 1, Last.DEL);
        } else if (scores[x - 1][y - 1][0] >= scores[x - 1][y - 1][1] &&
                scores[x - 1][y - 1][0] >= scores[x - 1][y - 1][2]) {
            scores[x][y][0] = (short) (scores[x - 1][y - 1][0] + getSubstitutionScore(x, y));
            traceback[x][y][0] = new AffineCross(y - 1, Last.SUB);
        } else {
            scores[x][y][0] = (short) (scores[x - 1][y - 1][2] + getSubstitutionScore(x, y));
            traceback[x][y][0] = new AffineCross(y - 1, Last.INS);
        }

        // deletion
        if (scores[x - 1][y][1] >= scores[x - 1][y][0] + gapPenalty.getOpenPenalty()) {
            scores[x][y][1] = (short) (scores[x - 1][y][1] + gapPenalty.getExtensionPenalty());
            traceback[x][y][1] = new AffineCross(y, Last.DEL);
        } else {
            scores[x][y][1] = (short) (scores[x - 1][y][0] + gapPenalty.getOpenPenalty() +
                    gapPenalty.getExtensionPenalty());
            traceback[x][y][1] = new AffineCross(y, Last.SUB);
        }

        // insertion
        if (scores[x][y - 1][0] + gapPenalty.getOpenPenalty() >= scores[x][y - 1][2]) {
            scores[x][y][2] = (short) (scores[x][y - 1][0] + gapPenalty.getOpenPenalty() +
                    gapPenalty.getExtensionPenalty());
            traceback[x][y][2] = new AffineCross((traceback[x][y - 1][0] == null) ? 0 : traceback[x][y - 1][0].y,
                    Last.SUB);
        } else {
            scores[x][y][2] = (short) (scores[x][y - 1][2] + gapPenalty.getExtensionPenalty());
            traceback[x][y][2] = new AffineCross((traceback[x][y - 1][2] == null) ? 0 : traceback[x][y - 1][2].y,
                    Last.INS);
        }

    }

    // scores alignment for a given position in both sequences for linear gap penalty
    protected void setScorePointLinear(int x, int y, Cross[][][] traceback, int sub, int del, int ins) {
        if (del >= sub && del >= ins) {
            scores[x][y][0] = (short) del;
            traceback[x][y][0] = new Cross(y);
        } else if (sub >= del && sub >= ins) {
            scores[x][y][0] = (short) sub;
            traceback[x][y][0] = new Cross(y - 1);
        } else {
            scores[x][y][0] = (short) ins;
            traceback[x][y][0] = traceback[x][y - 1][0];
        }
    }

    // scores global alignment for a given position in the query sequence
    protected void setScoreVectorGlobal(int x, Cross[][][] traceback) {
        int[] dim = getScoreMatrixDimensions();
        scores[x] = new short[dim[1]][dim[2]];
        traceback[x] = new AffineCross[dim[1]][dim[2]];
        short min = (short) (Short.MIN_VALUE - gapPenalty.getOpenPenalty() - gapPenalty.getExtensionPenalty());
        if (x == 0) {
            scores[0][0][1] = scores[0][0][2] = gapPenalty.getOpenPenalty();
            for (int y = 1; y < dim[1]; y++) {
                scores[0][y][0] = scores[0][y][1] = min;
                scores[0][y][2] = (short) (scores[0][y - 1][2] + gapPenalty.getExtensionPenalty());
                traceback[0][y][0] = traceback[0][y][1] = traceback[0][y][2] = new AffineCross(0, Last.INS);
            }
            xMax = scores.length - 1;
            yMax = dim[1] - 1;
        } else {
            scores[x][0][0] = scores[x][0][2] = min;
            scores[x][0][1] = (short) (scores[x - 1][0][1] + gapPenalty.getExtensionPenalty());
            traceback[x][0][0] = traceback[x][0][1] = traceback[x][0][2] = new AffineCross(0, Last.DEL);
            for (int y = 1; y < dim[1]; y++) {
                setScorePoint(x, y, traceback);
            }
            if (!storingScoreMatrix) {
                scores[x - 1] = null;
            }
        }
    }

    // scores local alignment for a given position in the query sequence
    protected void setScoreVectorLocal(int x, Cross[][][] traceback) {
        int[] dim = getScoreMatrixDimensions();
        scores[x] = new short[dim[1]][dim[2]];
        traceback[x] = new AffineCross[dim[1]][dim[2]];
        if (x > 0) {
            for (int y = 1; y < dim[1]; y++) {
                setScorePoint(x, y, traceback);
                for (int z = 0; z < dim[2]; z++) {
                    if (scores[x][y][z] <= 0) {
                        scores[x][y][z] = 0;
                        traceback[x][y][z] = null;
                    }
                }
                if (scores[x][y][0] > maxScore) {
                    xMax = x;
                    yMax = y;
                    maxScore = scores[x][y][0];
                }
            }
            if (!isStoringScoreMatrix()) {
                scores[x - 1] = null;
            }
        }
    }

    // scores global alignment for a given position in the query sequence for a linear gap penalty
    protected void setScoreVectorLinearGlobal(int x, Cross[][][] traceback) {
        int[] dim = getScoreMatrixDimensions();
        scores[x] = new short[dim[1]][dim[2]];
        traceback[x] = new Cross[dim[1]][dim[2]];
        if (x == 0) {
            for (int y = 1; y < dim[1]; y++) {
                scores[0][y][0] = (short) (scores[0][y - 1][0] + gapPenalty.getExtensionPenalty());
                traceback[0][y][0] = traceback[0][y - 1][0];
            }
            xMax = scores.length - 1;
            yMax = dim[1] - 1;
        } else {
            scores[x][0][0] = (short) (scores[x - 1][0][0] + gapPenalty.getExtensionPenalty());
            traceback[x][0][0] = new Cross(0);
            for (int y = 1; y < dim[1]; y++) {
                int del = scores[x - 1][y][0] + gapPenalty.getExtensionPenalty(),
                        ins = scores[x][y - 1][0] + gapPenalty.getExtensionPenalty(),
                        sub = scores[x - 1][y - 1][0] + getSubstitutionScore(x, y);
                setScorePointLinear(x, y, traceback, sub, del, ins);
            }
            if (!storingScoreMatrix) {
                scores[x - 1] = null;
            }
        }
    }

    // scores local alignment for a given position in the query sequence for a linear gap penalty
    protected void setScoreVectorLinearLocal(int x, Cross[][][] traceback) {
        int[] dim = getScoreMatrixDimensions();
        scores[x] = new short[dim[1]][dim[2]];
        traceback[x] = new Cross[dim[1]][dim[2]];
        GapPenalty gapPenalty = getGapPenalty();
        if (x > 0) {
            for (int y = 1; y < dim[1]; y++) {
                int del = scores[x - 1][y][0] + gapPenalty.getExtensionPenalty(),
                        ins = scores[x][y - 1][0] + gapPenalty.getExtensionPenalty(),
                        sub = scores[x - 1][y - 1][0] + getSubstitutionScore(x, y);
                if (del <= 0 && sub <= 0 && ins <= 0) {
                    scores[x][y][0] = (short) 0;
                } else {
                    setScorePointLinear(x, y, traceback, sub, del, ins);
                }
                if (scores[x][y][0] > scores[xMax][yMax][0]) {
                    xMax = x;
                    yMax = y;
                }
            }
            if (!isStoringScoreMatrix()) {
                scores[x - 1] = null;
            }
        }
    }

    // finds alignment path through traceback matrix
    protected void setSteps(Cross[][][] traceback, List<Step> sx, List<Step> sy, Last last) {
        int x = xMax, y = yMax;
        while (traceback[x][y][0] != null) {
            switch (last) {
            case DEL:
                if (traceback[x][y][1] instanceof AffineCross) {
                    last = ((AffineCross) traceback[x][y][1]).last;
                }
                sx.add(0, Step.COMPOUND);
                sy.add(0, Step.GAP);
                x--;
                break;
            case SUB:
                if (traceback[x][y][0] instanceof AffineCross) {
                    last = ((AffineCross) traceback[x][y][0]).last;
                }
                sx.add(0, Step.COMPOUND);
                sy.add(0, Step.COMPOUND);
                x--;
                y--;
                break;
            case INS:
                if (traceback[x][y][2] instanceof AffineCross) {
                    last = ((AffineCross) traceback[x][y][2]).last;
                }
                sx.add(0, Step.GAP);
                sy.add(0, Step.COMPOUND);
                y--;
            }
        }
        xStart = x;
        yStart = y;
    }

    // finds global alignment path through traceback matrix
    protected void setStepsGlobal(Cross[][][] traceback, List<Step> sx, List<Step> sy) {
        setSteps(traceback, sx, sy, (scores[xMax][yMax][1] > scores[xMax][yMax][0] && scores[xMax][yMax][1] >
                scores[xMax][yMax][2]) ? Last.DEL : (scores[xMax][yMax][0] > scores[xMax][yMax][1] &&
                scores[xMax][yMax][0] > scores[xMax][yMax][2]) ? Last.SUB : Last.INS);
    }

    // finds local alignment path through traceback matrix
    protected void setStepsLocal(Cross[][][] traceback, List<Step> sx, List<Step> sy) {
        setSteps(traceback, sx, sy, Last.SUB);
    }

    // finds alignment path through traceback matrix for linear gap penalty
    protected void setStepsLinear(Cross[][][] traceback, List<Step> sx, List<Step> sy) {
        int x = xMax, y = yMax;
        while (traceback[x][y][0] != null) {
            if (traceback[x][y][0].y == y) {
                sx.add(0, Step.COMPOUND);
                sy.add(0, Step.GAP);
                x--;
            } else if (traceback[x][y][0].y == y - 1) {
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

    // abstract methods

    // returns compound set of sequences
    protected abstract CompoundSet<C> getCompoundSet();

    // returns compounds in query sequence/profile
    protected abstract List<C> getCompoundsOfQuery();

    // returns compounds in target sequence/profile
    protected abstract List<C> getCompoundsOfTarget();

    // returns the 3 score matrix dimensions
    protected abstract int[] getScoreMatrixDimensions();

    // returns score for the alignment of two columns
    protected abstract short getSubstitutionScore(int queryColumn, int targetColumn);

    // prepares for alignment; returns true if everything is set to run the alignment
    protected abstract boolean isReady();

}

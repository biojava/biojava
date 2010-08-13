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
import java.util.Collections;
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
    protected short max, min, score;
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
        for (int type = 0; type < scores[0][0].length; type++) {
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

        int[] dim = getScoreMatrixDimensions();
        if (storingScoreMatrix) {
            scores = new short[dim[0]][dim[1]][dim[2]];
        } else {
            scores = new short[dim[0]][][];
            scores[0] = new short[dim[1]][dim[2]];
            scores[1] = new short[dim[1]][dim[2]];
        }
        List<Step> sx = new ArrayList<Step>(), sy = new ArrayList<Step>();
        if (gapPenalty.getType() == GapPenalty.Type.LINEAR) {
            types = new String[] { null };
            int[][][] traceback = new int[dim[0]][dim[1]][dim[2]];
            if (local) {
                for (int x = 1; x < dim[0]; x++) {
                    setScoreVectorLinearLocal(x, traceback);
                }
            } else {
                for (int x = 0; x < dim[0]; x++) {
                    setScoreVectorLinearGlobal(x, traceback);
                }
            }
            setStepsLinear(traceback, sx, sy);
        } else {
            types = new String[] { "Substitution", "Deletion", "Insertion" };
            Last[][][] traceback = new Last[dim[0]][dim[1]][dim[2]];
            if (local) {
                for (int x = 1; x < dim[0]; x++) {
                    setScoreVectorLocal(x, traceback);
                }
                setStepsLocal(traceback, sx, sy);
            } else {
                for (int x = 0; x < dim[0]; x++) {
                    setScoreVectorGlobal(x, traceback);
                }
                setStepsGlobal(traceback, sx, sy);
            }
        }
        if (!local) {
            for (int z = 0; z < scores[xMax][yMax].length; z++) {
                score = (short) Math.max(score, scores[xMax][yMax][z]);
            }
        }
        setProfile(sx, sy);

        time = System.nanoTime() - timeStart;

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

    // scores alignment for a given position in both sequences
    private void setScorePoint(int x, int y, Last[][][] traceback) {

        // substitution
        if (scores[x - 1][y - 1][1] >= scores[x - 1][y - 1][0] && scores[x - 1][y - 1][1] >= scores[x - 1][y - 1][2]) {
            scores[x][y][0] = (short) (scores[x - 1][y - 1][1] + getSubstitutionScore(x, y));
            traceback[x][y][0] = Last.DEL;
        } else if (scores[x - 1][y - 1][0] >= scores[x - 1][y - 1][1] &&
                scores[x - 1][y - 1][0] >= scores[x - 1][y - 1][2]) {
            scores[x][y][0] = (short) (scores[x - 1][y - 1][0] + getSubstitutionScore(x, y));
            traceback[x][y][0] = Last.SUB;
        } else {
            scores[x][y][0] = (short) (scores[x - 1][y - 1][2] + getSubstitutionScore(x, y));
            traceback[x][y][0] = Last.INS;
        }

        // deletion
        if (scores[x - 1][y][1] >= scores[x - 1][y][0] + gapPenalty.getOpenPenalty()) {
            scores[x][y][1] = (short) (scores[x - 1][y][1] + gapPenalty.getExtensionPenalty());
            traceback[x][y][1] = Last.DEL;
        } else {
            scores[x][y][1] = (short) (scores[x - 1][y][0] + gapPenalty.getOpenPenalty() +
                    gapPenalty.getExtensionPenalty());
            traceback[x][y][1] = Last.SUB;
        }

        // insertion
        if (scores[x][y - 1][0] + gapPenalty.getOpenPenalty() >= scores[x][y - 1][2]) {
            scores[x][y][2] = (short) (scores[x][y - 1][0] + gapPenalty.getOpenPenalty() +
                    gapPenalty.getExtensionPenalty());
            traceback[x][y][2] = Last.SUB;
        } else {
            scores[x][y][2] = (short) (scores[x][y - 1][2] + gapPenalty.getExtensionPenalty());
            traceback[x][y][2] = Last.INS;
        }

    }

    // scores alignment for a given position in both sequences for linear gap penalty
    private void setScorePointLinear(int x, int y, int[][][] traceback, int sub, int del, int ins) {
        if (del >= sub && del >= ins) {
            scores[x][y][0] = (short) del;
            traceback[x][y][0] = y;
        } else if (sub >= del && sub >= ins) {
            scores[x][y][0] = (short) sub;
            traceback[x][y][0] = y - 1;
        } else {
            scores[x][y][0] = (short) ins;
            traceback[x][y][0] = traceback[x][y - 1][0];
        }
    }

    // scores global alignment for a given position in the query sequence
    private void setScoreVectorGlobal(int x, Last[][][] traceback) {
        if (x == 0) {
            short min = (short) (Short.MIN_VALUE - gapPenalty.getOpenPenalty() - gapPenalty.getExtensionPenalty());
            scores[0][0][1] = scores[0][0][2] = gapPenalty.getOpenPenalty();
            for (int y = 1; y < scores[0].length; y++) {
                scores[0][y][0] = scores[0][y][1] = min;
                scores[0][y][2] = (short) (scores[0][y - 1][2] + gapPenalty.getExtensionPenalty());
                traceback[0][y][0] = traceback[0][y][1] = traceback[0][y][2] = Last.INS;
            }
            xMax = scores.length - 1;
            yMax = scores[0].length - 1;
        } else {
            if (!storingScoreMatrix && x > 1) {
                scores[x] = scores[x - 2];
            }
            scores[x][0][0] = scores[x][0][2] = (x > 1) ? scores[x - 1][0][0] : scores[0][1][0];
            scores[x][0][1] = (short) (scores[x - 1][0][1] + gapPenalty.getExtensionPenalty());
            traceback[x][0][0] = traceback[x][0][1] = traceback[x][0][2] = Last.DEL;
            for (int y = 1; y < scores[x].length; y++) {
                setScorePoint(x, y, traceback);
            }
        }
    }

    // scores local alignment for a given position in the query sequence
    private void setScoreVectorLocal(int x, Last[][][] traceback) {
        if (!storingScoreMatrix && x > 1) {
            scores[x] = scores[x - 2];
        }
        for (int y = 1; y < scores[x].length; y++) {
            setScorePoint(x, y, traceback);
            for (int z = 0; z < scores[x][y].length; z++) {
                if (scores[x][y][z] <= 0) {
                    scores[x][y][z] = 0;
                    traceback[x][y][z] = null;
                }
            }
            if (scores[x][y][0] > score) {
                xMax = x;
                yMax = y;
                score = scores[x][y][0];
            }
        }
    }

    // scores global alignment for a given position in the query sequence for a linear gap penalty
    private void setScoreVectorLinearGlobal(int x, int[][][] traceback) {
        if (x == 0) {
            for (int y = 1; y < scores[0].length; y++) {
                scores[0][y][0] = (short) (scores[0][y - 1][0] + gapPenalty.getExtensionPenalty());
            }
            xMax = scores.length - 1;
            yMax = scores[0].length - 1;
        } else {
            if (!storingScoreMatrix && x > 1) {
                scores[x] = scores[x - 2];
            }
            scores[x][0][0] = (short) (scores[x - 1][0][0] + gapPenalty.getExtensionPenalty());
            for (int y = 1; y < scores[x].length; y++) {
                int del = scores[x - 1][y][0] + gapPenalty.getExtensionPenalty(),
                        ins = scores[x][y - 1][0] + gapPenalty.getExtensionPenalty(),
                        sub = scores[x - 1][y - 1][0] + getSubstitutionScore(x, y);
                setScorePointLinear(x, y, traceback, sub, del, ins);
            }
        }
    }

    // scores local alignment for a given position in the query sequence for a linear gap penalty
    private void setScoreVectorLinearLocal(int x, int[][][] traceback) {
        if (!storingScoreMatrix && x > 1) {
            scores[x] = scores[x - 2];
        }
        for (int y = 1; y < scores[x].length; y++) {
            int del = scores[x - 1][y][0] + gapPenalty.getExtensionPenalty(),
                    ins = scores[x][y - 1][0] + gapPenalty.getExtensionPenalty(),
                    sub = scores[x - 1][y - 1][0] + getSubstitutionScore(x, y);
            if (del > 0 || sub > 0 || ins > 0) {
                setScorePointLinear(x, y, traceback, sub, del, ins);
                if (scores[x][y][0] > scores[xMax][yMax][0]) {
                    xMax = x;
                    yMax = y;
                }
            }
        }
    }

    // finds alignment path through traceback matrix
    private void setSteps(Last[][][] traceback, List<Step> sx, List<Step> sy, Last last) {
        int x = xMax, y = yMax;
        while (local ? traceback[x][y][last.ordinal()] != null : x > 0 || y > 0) {
            switch (last) {
            case DEL:
                last = traceback[x][y][1];
                sx.add(Step.COMPOUND);
                sy.add(Step.GAP);
                x--;
                break;
            case SUB:
                last = traceback[x][y][0];
                sx.add(Step.COMPOUND);
                sy.add(Step.COMPOUND);
                x--;
                y--;
                break;
            case INS:
                last = traceback[x][y][2];
                sx.add(Step.GAP);
                sy.add(Step.COMPOUND);
                y--;
            }
        }
        Collections.reverse(sx);
        Collections.reverse(sy);
        xStart = x;
        yStart = y;
    }

    // finds global alignment path through traceback matrix
    private void setStepsGlobal(Last[][][] traceback, List<Step> sx, List<Step> sy) {
        setSteps(traceback, sx, sy, (scores[xMax][yMax][1] > scores[xMax][yMax][0] && scores[xMax][yMax][1] >
                scores[xMax][yMax][2]) ? Last.DEL : (scores[xMax][yMax][0] > scores[xMax][yMax][1] &&
                scores[xMax][yMax][0] > scores[xMax][yMax][2]) ? Last.SUB : Last.INS);
    }

    // finds local alignment path through traceback matrix
    private void setStepsLocal(Last[][][] traceback, List<Step> sx, List<Step> sy) {
        setSteps(traceback, sx, sy, Last.SUB);
    }

    // finds alignment path through traceback matrix for linear gap penalty
    private void setStepsLinear(int[][][] traceback, List<Step> sx, List<Step> sy) {
        int x = xMax, y = yMax;
        while (local ? traceback[x][y][0] > 0 : x > 0 || y > 0) {
            if (traceback[x][y][0] == y) {
                sx.add(Step.COMPOUND);
                sy.add(Step.GAP);
                x--;
            } else if (traceback[x][y][0] == y - 1) {
                sx.add(Step.COMPOUND);
                sy.add(Step.COMPOUND);
                x--;
                y--;
            } else {
                sx.add(Step.GAP);
                sy.add(Step.COMPOUND);
                y--;
            }
        }
        Collections.reverse(sx);
        Collections.reverse(sy);
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

    // sets profile following the given alignment path
    protected abstract void setProfile(List<Step> sx, List<Step> sy);

}

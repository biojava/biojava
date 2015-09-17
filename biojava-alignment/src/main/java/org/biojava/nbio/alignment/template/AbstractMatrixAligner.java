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

package org.biojava.nbio.alignment.template;

import org.biojava.nbio.alignment.routines.AlignerHelper.Anchor;
import org.biojava.nbio.alignment.routines.AlignerHelper.Last;
import org.biojava.nbio.alignment.routines.AlignerHelper.Subproblem;
import org.biojava.nbio.alignment.template.AlignedSequence.Step;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.CompoundSet;
import org.biojava.nbio.core.sequence.template.Sequence;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.biojava.nbio.alignment.routines.AlignerHelper.setScoreVector;
import static org.biojava.nbio.alignment.routines.AlignerHelper.setSteps;

/**
 * Implements common code for an {@link Aligner} which builds a score matrix during computation.
 *
 * @author Mark Chapman
 * @author Daniel Cameron
 * @param <S> each element of the alignment {@link Profile} is of type S
 * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
 */
public abstract class AbstractMatrixAligner<S extends Sequence<C>, C extends Compound> extends AbstractScorer
        implements MatrixAligner<S, C> {

    // input fields
    protected GapPenalty gapPenalty;
    private SubstitutionMatrix<C> subMatrix;
    private boolean local, storingScoreMatrix;
    protected List<Anchor> anchors = new ArrayList<Anchor>();
    protected int cutsPerSection;

    // output fields
    protected Profile<S, C> profile;
    /**
     * Start position of the aligned sequence in the query and target respectively
     */
    protected int[] xyStart;
    /**
     * End position of the aligned sequence in the query and target respectively
     */
    protected int[] xyMax;
    protected int max, min, score;
    /**
     * Dynamic programming score matrix
     * The first dimension has the length of the first (query) sequence + 1
     * The second has the length of the second (target) sequence + 1
     * The third has length 1 for linear gap penalty and 3 for affine/constant gap
     * (one each for match/substitution, deletion, insertion) 
     */
    protected int[][][] scores;
    /**
     * Friendly name of each copy of the scoring matrix.
     * The number of elements must match the number of elements in third dimension of @see scores 
     */
    private String[] types;
    protected long time = -1;

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
    public int[][][] getScoreMatrix() {
        boolean tempStoringScoreMatrix = storingScoreMatrix;
        if (scores == null) {
            storingScoreMatrix = true;
            align();
            if (scores == null) {
                return null;
            }
        }
        int[][][] copy = scores;
        if (tempStoringScoreMatrix) {
            copy = new int[scores.length][scores[0].length][];
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
        int[][][] scores = getScoreMatrix();
        return scoreMatrixToString(scores);
    }
	private String scoreMatrixToString(int[][][] scores) {
		StringBuilder s = new StringBuilder();
        CompoundSet<C> compoundSet = getCompoundSet();
        int lengthCompound = compoundSet.getMaxSingleCompoundStringLength(), lengthRest =
                Math.max(Math.max(Integer.toString(min).length(), Integer.toString(max).length()), lengthCompound) + 1;
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
    public double getMaxScore() {
        if (profile == null) {
            align();
        }
        return max;
    }

    @Override
    public double getMinScore() {
        if (profile == null) {
            align();
        }
        return min;
    }

    @Override
    public double getScore() {
        if (profile == null) {
            align();
        }
        return score;
    }

    // helper methods

    /**
     * Performs alignment
     */
    protected void align() {
        if (!isReady()) {
            return;
        }

        long timeStart = System.nanoTime();

        int[] dim = getScoreMatrixDimensions();
        if (storingScoreMatrix) {
            scores = new int[dim[0]][dim[1]][dim[2]];
        } else {
            scores = new int[dim[0]][][];
            scores[0] = new int[dim[1]][dim[2]];
            scores[1] = new int[dim[1]][dim[2]];
        }
        boolean linear = (gapPenalty.getType() == GapPenalty.Type.LINEAR);
        Last[][][] traceback = new Last[dim[0]][][];
        List<Step> sx = new ArrayList<Step>(), sy = new ArrayList<Step>();

        if (!local) {
        	xyMax = new int[] { dim[0] - 1, dim[1] - 1 };
        	xyStart = new int[] { 0, 0 };
            score = 0;
            List<Subproblem> problems = Subproblem.getSubproblems(anchors, xyMax[0], xyMax[1]);
            assert problems.size() == anchors.size() + 1;
            for (int i = 0; i < problems.size(); i++) {
            	Subproblem subproblem = problems.get(i);
            	for (int x = subproblem.getQueryStartIndex(); x <= subproblem.getQueryEndIndex(); x++) {
                	
            		traceback[x] = 
                		
            		  linear ? 
                		setScoreVector(x, subproblem, gapPenalty.getExtensionPenalty(), getSubstitutionScoreVector(x, subproblem), storingScoreMatrix, scores) : 
                		
                		setScoreVector(x, subproblem, gapPenalty.getOpenPenalty(), gapPenalty.getExtensionPenalty(), getSubstitutionScoreVector(x, subproblem), storingScoreMatrix, scores);
                }
            }
            setSteps(traceback, scores, sx, sy);
            score = Integer.MIN_VALUE;
            int[] finalScore = scores[xyMax[0]][xyMax[1]];
            for (int z = 0; z < finalScore.length; z++) {
            	score = Math.max(score, finalScore[z]);
            }
        } else {
            for (int x = 0; x < dim[0]; x++) {
            	
            	traceback[x] = 
            			
            	  linear ? 
           			setScoreVector(x, gapPenalty.getExtensionPenalty(), getSubstitutionScoreVector(x), storingScoreMatrix, scores, xyMax, score) :
                    
           			setScoreVector(x, gapPenalty.getOpenPenalty(), gapPenalty.getExtensionPenalty(), getSubstitutionScoreVector(x), storingScoreMatrix, scores, xyMax, score);
           				
                if (xyMax[0] == x) {
                    score = scores[x][xyMax[1]][0];
                }
            }
            xyStart = local ? setSteps(traceback, xyMax, sx, sy) : setSteps(traceback, scores, sx, sy);
        }

        setProfile(sx, sy);

        if (!storingScoreMatrix) {
            scores = null;
        }

        time = System.nanoTime() - timeStart;
    }

    /**
     * Returns score for the alignment of the query column to all target columns
     * @param queryColumn
     * @return
     */
    protected int[] getSubstitutionScoreVector(int queryColumn) {
        return getSubstitutionScoreVector(queryColumn, new Subproblem(0, 0, scores.length - 1, scores[0].length - 1));
    }

    /**
     * Returns score for the alignment of the query column to all target columns
     * @param queryColumn
     * @param subproblem
     * @return
     */
    protected int[] getSubstitutionScoreVector(int queryColumn, Subproblem subproblem) {
        int[] subs = new int[subproblem.getTargetEndIndex() + 1];
        if (queryColumn > 0) {
            for (int y = Math.max(1, subproblem.getTargetStartIndex()); y <= subproblem.getTargetEndIndex(); y++) {
                subs[y] = getSubstitutionScore(queryColumn, y);
            }
        }
        return subs;
    }

    /**
     * Resets output fields; should be overridden to set max and min
     */
    protected void reset() {
        xyMax = new int[] {0, 0};
        xyStart = new int[] {0, 0};
        scores = null;
        types = (gapPenalty == null || gapPenalty.getType() == GapPenalty.Type.LINEAR) ? new String[] { null } :
                new String[] { "Substitution", "Deletion", "Insertion" };
        time = -1;
        profile = null;
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
    protected abstract int getSubstitutionScore(int queryColumn, int targetColumn);

    // prepares for alignment; returns true if everything is set to run the alignment
    protected abstract boolean isReady();

    // sets profile following the given alignment path
    protected abstract void setProfile(List<Step> sx, List<Step> sy);

}

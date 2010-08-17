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
 * Created on August 13, 2010
 * Author: Mark Chapman
 */

package org.biojava3.alignment.routines;

import java.util.Collections;
import java.util.List;

import org.biojava3.alignment.template.AlignedSequence.Step;

/**
 * Static utility to construct alignment routines from a common library of methods.
 *
 * @author Mark Chapman
 */
public class AlignerHelper {

    // types

    /**
     * Defines a traceback pointer for the three edit operations: substitution (match/replacement of a query compound
     * with a target compound), deletion (removal of a query compound leaving a gap in the target sequence), and
     * insertion (addition of a target compound opening a gap in the query sequence).
     */
    public enum Last {
        SUBSTITUTION,
        DELETION,
        INSERTION
    }

    /**
     * Defines a 'cut' row for divide-and-conquer alignment in which a new anchor is found.
     */
    public static class Cut {

        private int queryIndex;
        private int[][] targetIndices, tiLast, ti1, ti2;
        private Last[][] pointers, pLast, p1, p2;

        public Cut(int queryIndex, int[] dim) {
            this.queryIndex = queryIndex;
            targetIndices = ti1 = new int[dim[1]][dim[2]];
            ti2 = new int[dim[1]][dim[2]];
            pointers = p1 = new Last[dim[1]][dim[2]];
            p2 = new Last[dim[1]][dim[2]];
        }

        public Last getPointer(int z) {
            return pointers[pointers.length - 1][z];
        }

        public int getQueryIndex() {
            return queryIndex;
        }

        public int getTargetIndex(int z) {
            return targetIndices[targetIndices.length - 1][z];
        }

        public void update(int x, int[] subproblem, Last[][] pointerUpdates) {
            if (queryIndex == x - 1) {
                // TODO create
                for (int y = subproblem[1]; y <= subproblem[3]; y++) {
                    for (int z = 0; z < pointers[y].length; z++) {
                        pointers[y][z] = pointerUpdates[y][z];
                        if (pointers[y][z] != null) {
                            switch (pointers[y][z]) {
                            case DELETION:
                                targetIndices[y][z] = y;
                                break;
                            case SUBSTITUTION:
                                targetIndices[y][z] = y - 1;
                                break;
                            case INSERTION:
                                targetIndices[y][z] = targetIndices[y - 1][z];
                            }
                        }
                    }
                }
            } else if (queryIndex < x) {
                // TODO advance
                tiLast = targetIndices;
                targetIndices = (targetIndices == ti2) ? ti1 : ti2;
                pLast = pointers;
                pointers = (pointers == p2) ? p1 : p2;
                for (int y = subproblem[1]; y <= subproblem[3]; y++) {
                    for (int z = 0; z < pointers[y].length; z++) {
                        switch (pointerUpdates[y][z]) {
                        case DELETION:
                            targetIndices[y][z] = tiLast[y][z];
                            pointers[y][z] = pLast[y][z];
                            break;
                        case SUBSTITUTION:
                            targetIndices[y][z] = tiLast[y - 1][z];
                            pointers[y][z] = pLast[y - 1][z];
                            break;
                        case INSERTION:
                            targetIndices[y][z] = targetIndices[y - 1][z];
                            pointers[y][z] = pointers[y - 1][z];
                        }
                    }
                }
            }
        }

    }

    /**
     * Implements a data structure for linking query compounds and target compounds.
     */
    public static class Anchor {

        private int targetIndex;
        private Last pointer;

        public Anchor() {
        }

        public Anchor(Cut c, int z) {
            this(c.getTargetIndex(z), c.getPointer(z));
        }

        public Anchor(int y, Last p) {
            targetIndex = y;
            pointer = p;
        }

        public Last getPointer() {
            return pointer;
        }

        public int getTargetIndex() {
            return targetIndex;
        }

    }

    // methods

    public static void addAnchors(Cut[] cuts, Anchor[] anchors, int z) {
        for (Cut c : cuts) {
            anchors[c.getQueryIndex()] = new Anchor(c, z);
        }
    }

    public static Cut[] getCuts(int k, int[] subproblem, int[] dim) {
        Cut[] cuts;
        int m = subproblem[2] - subproblem[0] - 1;
        if (k < m) {
            cuts = new Cut[k];
            for (int i = 0; i < k; i++) {
                cuts[i] = new Cut(subproblem[0] + ((i + 1) * (m + 1) * k / (k + 1)), dim);
            }
        } else {
            cuts = new Cut[m];
            for (int i = 0, x = subproblem[0] + 1; i < m; i++, x++) {
                cuts[i] = new Cut(x, dim);
            }
        }
        return cuts;
    }

    /**
     * Returns the coordinates for the next subproblem.  The first element is the starting index in the query sequence.
     * The second element is the starting index in the target sequence.  Similarly, the third and fourth elements are
     * the ending indices of the query and target sequences, respectively.  When the alignment is finished, there is no
     * next subproblem and this method returns {@code null}.
     *
     * @param anchors current list of anchors
     * @param targetLength length of the target sequence
     * @return the coordinates for the next subproblem
     */
    public static int[] getNextSubproblem(Anchor[] anchors, int targetLength) {
        int[] subproblem = new int[4];
        int x = 0;

        // find first unanchored x (query sequence index)
        while (x < anchors.length && anchors[x] != null) {
            x++;
        }
        if (x == anchors.length) {
            return null; // no unanchored x, therefore alignment is complete
        }

        // save last anchor as starting point
        subproblem[0] = x - 1;
        subproblem[1] = anchors[x - 1].getTargetIndex();

        // find next anchored x or reach the end
        while (x < anchors.length && anchors[x] == null) {
            x++;
        }
        if (x == anchors.length) {
            subproblem[2] = x - 1;
            subproblem[3] = targetLength;
        } else {
            subproblem[2] = x;
            subproblem[3] = anchors[x].getTargetIndex();
        }

        return subproblem;
    }

    // updates cut rows given the latest row of traceback pointers
    public static void setCuts(int x, int[] subproblem, Last[][] pointers, Cut[]cuts) {
        for (Cut c : cuts) {
            c.update(x, subproblem, pointers);
        }
    }

    // scores alignment for a given position in both sequences
    public static Last[] setScorePoint(int x, int y, short gop, short gep, short sub, short[][][] scores) {
        Last[] pointers = new Last[3];

        // substitution
        if (scores[x - 1][y - 1][1] >= scores[x - 1][y - 1][0] && scores[x - 1][y - 1][1] >= scores[x - 1][y - 1][2]) {
            scores[x][y][0] = (short) (scores[x - 1][y - 1][1] + sub);
            pointers[0] = Last.DELETION;
        } else if (scores[x - 1][y - 1][0] >= scores[x - 1][y - 1][2]) {
            scores[x][y][0] = (short) (scores[x - 1][y - 1][0] + sub);
            pointers[0] = Last.SUBSTITUTION;
        } else {
            scores[x][y][0] = (short) (scores[x - 1][y - 1][2] + sub);
            pointers[0] = Last.INSERTION;
        }

        // deletion
        if (scores[x - 1][y][1] >= scores[x - 1][y][0] + gop) {
            scores[x][y][1] = (short) (scores[x - 1][y][1] + gep);
            pointers[1] = Last.DELETION;
        } else {
            scores[x][y][1] = (short) (scores[x - 1][y][0] + gop + gep);
            pointers[1] = Last.SUBSTITUTION;
        }

        // insertion
        if (scores[x][y - 1][0] + gop >= scores[x][y - 1][2]) {
            scores[x][y][2] = (short) (scores[x][y - 1][0] + gop + gep);
            pointers[2] = Last.SUBSTITUTION;
        } else {
            scores[x][y][2] = (short) (scores[x][y - 1][2] + gep);
            pointers[2] = Last.INSERTION;
        }

        return pointers;
    }

    // scores alignment for a given position in both sequences for linear gap penalty
    public static Last setScorePoint(int x, int y, short gep, short sub, short[][][] scores) {
        int d = scores[x - 1][y][0] + gep, i = scores[x][y - 1][0] + gep, s = scores[x - 1][y - 1][0] + sub;
        if (d >= s && d >= i) {
            scores[x][y][0] = (short) d;
            return Last.DELETION;
        } else if (s >= i) {
            scores[x][y][0] = (short) s;
            return Last.SUBSTITUTION;
        } else {
            scores[x][y][0] = (short) i;
            return Last.INSERTION;
        }
    }

    // scores global alignment for a given position in the query sequence
    public static Last[][] setScoreVector(int x, short gop, short gep, short[] subs, boolean storing,
            short[][][] scores) {
        return setScoreVector(x, 0, 0, scores[0].length - 1, gop, gep, subs, storing, scores);
    }

    // scores global alignment for a given position in the query sequence
    public static Last[][] setScoreVector(int x, int[] subproblem, short gop, short gep, short[] subs, boolean storing,
            short[][][] scores) {
        return setScoreVector(x, subproblem[0], subproblem[1], subproblem[3], gop, gep, subs, storing, scores);
    }

    // scores global alignment for a given position in the query sequence
    public static Last[][] setScoreVector(int x, int xb, int yb, int ye, short gop, short gep, short[] subs,
            boolean storing, short[][][] scores) {
        Last[][] pointers = new Last[ye + 1][];
        short min = (short) (Short.MIN_VALUE - gop - gep);
        if (x == xb) {
            scores[xb][yb][1] = scores[xb][yb][2] = gop;
            pointers[yb] = new Last[] {null, null, null};
            Last[] insertions = new Last[] { Last.INSERTION, Last.INSERTION, Last.INSERTION };
            for (int y = yb + 1; y <= ye; y++) {
                scores[xb][y][0] = scores[xb][y][1] = min;
                scores[xb][y][2] = (short) (scores[xb][y - 1][2] + gep);
                pointers[y] = insertions;
            }
        } else {
            if (!storing && x > xb + 1) {
                scores[x] = scores[x - 2];
            }
            scores[x][yb][0] = scores[x][yb][2] = min;
            scores[x][yb][1] = (short) (scores[x - 1][yb][1] + gep);
            pointers[yb] = new Last[] { Last.DELETION, Last.DELETION, Last.DELETION };
            for (int y = yb + 1; y <= ye; y++) {
                pointers[y] = setScorePoint(x, y, gop, gep, subs[y], scores);
            }
        }
        return pointers;
    }

    // scores global alignment for a given position in the query sequence for a linear gap penalty
    public static Last[][] setScoreVector(int x, short gep, short[] subs, boolean storing, short[][][] scores) {
        return setScoreVector(x, 0, scores[0].length - 1, gep, subs, storing, scores);
    }

    // scores global alignment for a given position in the query sequence for a linear gap penalty
    public static Last[][] setScoreVector(int x, int[] subproblem, short gep, short[] subs, boolean storing,
            short[][][] scores) {
        return setScoreVector(x, subproblem[1], subproblem[3], gep, subs, storing, scores);
    }

    // scores global alignment for a given position in the query sequence for a linear gap penalty
    public static Last[][] setScoreVector(int x, int yb, int ye, short gep, short[] subs, boolean storing,
            short[][][] scores) {
        Last[][] pointers = new Last[ye + 1][1];
        if (x == 0) {
            for (int y = yb + 1; y <= ye; y++) {
                scores[0][y][0] = (short) (scores[0][y - 1][0] + gep);
                pointers[y][0] = Last.INSERTION;
            }
        } else {
            if (!storing && x > 1) {
                scores[x] = scores[x - 2];
            }
            scores[x][yb][0] = (short) (scores[x - 1][yb][0] + gep);
            pointers[yb][0] = Last.DELETION;
            for (int y = yb + 1; y <= ye; y++) {
                pointers[y][0] = setScorePoint(x, y, gep, subs[y], scores);
            }
        }
        return pointers;
    }

    // scores local alignment for a given position in the query sequence
    public static Last[][] setScoreVector(int x, short gop, short gep, short[] subs, boolean storing,
            short[][][] scores, int[] xyMax, int score) {
        return setScoreVector(x, 0, scores[0].length - 1, gop, gep, subs, storing, scores, xyMax, score);
    }

    // scores local alignment for a given position in the query sequence
    public static Last[][] setScoreVector(int x, int yb, int ye, short gop, short gep, short[] subs, boolean storing,
            short[][][] scores, int[] xyMax, int score) {
        Last[][] pointers;
        if (x == 0) {
            pointers = new Last[ye + 1][scores[0][0].length];
        } else {
            pointers = new Last[ye + 1][];
            pointers[0] = new Last[scores[0][0].length];
            if (!storing && x > 1) {
                scores[x] = scores[x - 2];
            }
            for (int y = 1; y < scores[0].length; y++) {
                pointers[y] = setScorePoint(x, y, gop, gep, subs[y], scores);
                for (int z = 0; z < scores[0][0].length; z++) {
                    if (scores[x][y][z] <= 0) {
                        scores[x][y][z] = 0;
                        pointers[y][z] = null;
                    }
                }
                if (scores[x][y][0] > score) {
                    xyMax[0] = x;
                    xyMax[1] = y;
                    score = scores[x][y][0];
                }
            }
        }
        return pointers;
    }

    // scores local alignment for a given position in the query sequence for a linear gap penalty
    public static Last[][] setScoreVector(int x, short gep, short[] subs, boolean storing, short[][][] scores,
            int[] xyMax, int score) {
        return setScoreVector(x, 0, scores[0].length - 1, gep, subs, storing, scores, xyMax, score);
    }

    // scores local alignment for a given position in the query sequence for a linear gap penalty
    public static Last[][] setScoreVector(int x, int yb, int ye, short gep, short[] subs, boolean storing,
            short[][][] scores, int[] xyMax, int score) {
        Last[][] pointers;
        if (x == 0) {
            pointers = new Last[ye + 1][1];
        } else {
            pointers = new Last[ye + 1][];
            pointers[0] = new Last[1];
            if (!storing && x > 1) {
                scores[x] = scores[x - 2];
            }
            for (int y = 1; y < scores[x].length; y++) {
                pointers[y][0] = setScorePoint(x, y, gep, subs[y], scores);
                if (scores[x][y][0] <= 0) {
                    scores[x][y][0] = 0;
                    pointers[y][0] = null;
                } else if (scores[x][y][0] > score) {
                    xyMax[0] = x;
                    xyMax[1] = y;
                    score = scores[x][y][0];
                }
            }
        }
        return pointers;
    }

    /**
     * Sets the alignment path following the list of anchors.
     *
     * @param anchors a complete list of anchors (input)
     * @param sx steps for the query sequence (output)
     * @param sy steps for the target sequence (output)
     */
    public static void setSteps(Anchor[] anchors, List<Step> sx, List<Step> sy) {
        for (int x = 1; x < anchors.length; x++) {
            switch (anchors[x].getPointer()) {
            case DELETION:
                sx.add(Step.COMPOUND);
                sy.add(Step.GAP);
                break;
            case SUBSTITUTION:
                sx.add(Step.COMPOUND);
                sy.add(Step.COMPOUND);
                break;
            case INSERTION:
                for (int dy = anchors[x].getTargetIndex() - anchors[x - 1].getTargetIndex(); dy > 1; dy--) {
                    sx.add(Step.GAP);
                    sy.add(Step.COMPOUND);
                }
                sx.add(Step.COMPOUND);
                sy.add(Step.COMPOUND);
            }
        }
    }

    // finds alignment path through traceback matrix
    public static int[] setSteps(Last[][][] traceback, boolean local, int[] xyMax, Last last, List<Step> sx,
            List<Step> sy) {
        int x = xyMax[0], y = xyMax[1];
        boolean linear = (traceback[x][y].length == 1);
        while (local ? (linear ? last : traceback[x][y][last.ordinal()]) != null : x > 0 || y > 0) {
            switch (last) {
            case DELETION:
                sx.add(Step.COMPOUND);
                sy.add(Step.GAP);
                last = linear ? traceback[--x][y][0] : traceback[x--][y][1];
                break;
            case SUBSTITUTION:
                sx.add(Step.COMPOUND);
                sy.add(Step.COMPOUND);
                last = linear ? traceback[--x][--y][0] : traceback[x--][y--][0];
                break;
            case INSERTION:
                sx.add(Step.GAP);
                sy.add(Step.COMPOUND);
                last = linear ? traceback[x][--y][0] : traceback[x][y--][2];
            }
        }
        Collections.reverse(sx);
        Collections.reverse(sy);
        return new int[] {x, y};
    }

    // finds global alignment path through traceback matrix
    public static int[] setSteps(Last[][][] traceback, short[][][] scores, List<Step> sx, List<Step> sy) {
        int xMax = scores.length - 1, yMax = scores[xMax].length - 1;
        boolean linear = (traceback[xMax][yMax].length == 1);
        Last last = linear ? traceback[xMax][yMax][0] : (scores[xMax][yMax][1] > scores[xMax][yMax][0] &&
                scores[xMax][yMax][1] > scores[xMax][yMax][2]) ? Last.DELETION : (scores[xMax][yMax][0] >
                scores[xMax][yMax][2]) ? Last.SUBSTITUTION : Last.INSERTION;
        return setSteps(traceback, false, new int[] {xMax, yMax}, last, sx, sy);
    }

    // finds local alignment path through traceback matrix
    public static int[] setSteps(Last[][][] traceback, int[] xyMax, List<Step> sx, List<Step> sy) {
        return setSteps(traceback, true, xyMax, Last.SUBSTITUTION, sx, sy);
    }

}

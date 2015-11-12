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

package org.biojava.nbio.alignment.routines;

import org.biojava.nbio.core.alignment.template.AlignedSequence.Step;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
//import org.slf4j.Logger;
//import org.slf4j.LoggerFactory;


/**
 * Static utility to construct alignment routines from a common library of methods.
 *
 * @author Mark Chapman
 * @author Daniel Cameron
 */
public class AlignerHelper {

	//private static final Logger logger = LoggerFactory.getLogger(AlignerHelper.class);
	
    // types

    /**
     * Define a traceback pointer for the three edit operations: substitution (match/replacement of a query compound
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

        public Cut(int queryIndex, int[] dim) {
            this.queryIndex = queryIndex;
            targetIndices = ti1 = new int[dim[1]][dim[2]];
            ti2 = new int[dim[1]][dim[2]];
        }

        public int getQueryIndex() {
            return queryIndex;
        }

        public int getTargetIndex(int z) {
            return targetIndices[targetIndices.length - 1][z];
        }

        public void update(int x, Subproblem subproblem, Last[][] pointers) {
            if (pointers[subproblem.getTargetStartIndex()].length == 1) {
                if (queryIndex == x - 1) {
                    updateLinearInitial(subproblem, pointers);
                } else if (queryIndex < x) {
                    updateLinearAdvance(subproblem, pointers);
                }
            } else {
                if (queryIndex == x - 1) {
                    updateInitial(subproblem, pointers);
                } else if (queryIndex < x) {
                    updateAdvance(subproblem, pointers);
                }
            }
        }

        private void updateAdvance(Subproblem subproblem, Last[][] pointers) {
            tiLast = targetIndices;
            targetIndices = (targetIndices == ti2) ? ti1 : ti2;
            for (int y = subproblem.getTargetStartIndex(); y <= subproblem.getTargetEndIndex(); y++) {
                if (pointers[y][0] != null) {
                    targetIndices[y][0] = tiLast[y - 1][pointers[y][0].ordinal()];
                }
                if (pointers[y][1] != null) {
                    targetIndices[y][1] = tiLast[y][pointers[y][1].ordinal()];
                }
                if (pointers[y][2] != null) {
                    targetIndices[y][2] = targetIndices[y - 1][pointers[y][2].ordinal()];
                }
            }
        }

        private void updateInitial(Subproblem subproblem, Last[][] pointers) {
            for (int y = subproblem.getTargetStartIndex(); y <= subproblem.getTargetEndIndex(); y++) {
                if (pointers[y][0] != null) {
                    targetIndices[y][0] = y - 1;
                }
                if (pointers[y][1] != null) {
                    targetIndices[y][1] = y;
                }
                if (pointers[y][2] != null) {
                    targetIndices[y][2] = targetIndices[y - 1][2];
                }
            }
        }

        private void updateLinearAdvance(Subproblem subproblem, Last[][] pointers) {
            tiLast = targetIndices;
            targetIndices = (targetIndices == ti2) ? ti1 : ti2;
            for (int y = subproblem.getTargetStartIndex(); y <= subproblem.getTargetEndIndex(); y++) {
                switch (pointers[y][0]) {
                case DELETION:
                    targetIndices[y][0] = tiLast[y][0];
                    break;
                case SUBSTITUTION:
                    targetIndices[y][0] = tiLast[y - 1][0];
                    break;
                case INSERTION:
                    targetIndices[y][0] = targetIndices[y - 1][0];
                }
            }
        }

        private void updateLinearInitial(Subproblem subproblem, Last[][] pointers) {
            for (int y = subproblem.getTargetStartIndex(); y <= subproblem.getTargetEndIndex(); y++) {
                if (pointers[y][0] != null) {
                    switch (pointers[y][0]) {
                    case DELETION:
                        targetIndices[y][0] = y;
                        break;
                    case SUBSTITUTION:
                        targetIndices[y][0] = y - 1;
                        break;
                    case INSERTION:
                        targetIndices[y][0] = targetIndices[y - 1][0];
                    }
                }
            }
        }

    }
    public static int addAnchors(Cut[] cuts, int[] scores, boolean addScore, int[] anchors) {
        int zMax = 0, subscore = scores[0];
        for (int z = 1; z < scores.length; z++) {
            if (scores[z] > subscore) {
                zMax = z;
                subscore = scores[z];
            }
        }
        for (Cut c : cuts) {
            anchors[c.getQueryIndex()] = c.getTargetIndex(zMax);
        }
        return addScore ? (int) subscore : 0;
    }

    public static Cut[] getCuts(int k, Subproblem subproblem, int[] dim, boolean anchor0) {
        Cut[] cuts;
        int m = subproblem.getQueryEndIndex() - subproblem.getQueryStartIndex() - (anchor0 ? 1 : 0);
        if (k < m) {
            cuts = new Cut[k];
            int firstCutIndex = subproblem.getQueryStartIndex() + (anchor0 ? 1 : 0);
            for (int i = 0; i < k; i++) {            	
            	cuts[i] = new Cut(firstCutIndex + i * (m - 1) / (k - 1), dim);
            }
        } else {
            cuts = new Cut[m];
            for (int i = 0, x = subproblem.getQueryStartIndex() + (anchor0 ? 1 : 0); i < m; i++, x++) {
                cuts[i] = new Cut(x, dim);
            }
        }
        return cuts;
    }
    /**
     * Compounds in query and target sequences that must align
     * @author Daniel Cameron
     */
    public static class Anchor {
    	public int getQueryIndex() {
			return queryIndex;
		}
		public int getTargetIndex() {
			return targetIndex;
		}
		private final int queryIndex;
    	private final int targetIndex;
    	public Anchor(int queryIndex, int targetIndex) {
    		this.queryIndex = queryIndex;
    		this.targetIndex = targetIndex;
    	}
    	public static class QueryIndexComparator implements Comparator<Anchor> {
			@Override
			public int compare(Anchor o1, Anchor o2) {
				return o1.getQueryIndex() - o2.getQueryIndex();
			}
    	}
    }
    /**
     * Alignment subproblem. The bounds of the subproblem are the
     * indicies representing the inclusive bounds of the dynamic programming
     * alignment problem.
     * @author Daniel Cameron
     */
    public static class Subproblem {
    	public int getTargetStartIndex() {
			return targetStartIndex;
		}
		public int getQueryEndIndex() {
			return queryEndIndex;
		}
		public int getTargetEndIndex() {
			return targetEndIndex;
		}
		public int getQueryStartIndex() {
			return queryStartIndex;
		}
		/**
		 * Indicates whether the start query and start target index compounds
		 * are anchored to each other
		 * @return true if the compounds are anchored in the alignment, false otherwise
		 */
		public boolean isStartAnchored() {
			return isAnchored;
		}
		private int queryStartIndex; // [0]
    	private int targetStartIndex; // [1]
    	private int queryEndIndex; // [2]
    	private int targetEndIndex; // [3]
    	private boolean isAnchored;
    	public Subproblem(int queryStartIndex, int targetStartIndex, int queryEndIndex, int targetEndIndex) {
    		this(queryStartIndex, targetStartIndex, queryEndIndex, targetEndIndex, false);
    	}
    	public Subproblem(int queryStartIndex, int targetStartIndex, int queryEndIndex, int targetEndIndex, boolean isAnchored) {
    		this.queryStartIndex = queryStartIndex;
    		this.targetStartIndex = targetStartIndex;
    		this.queryEndIndex = queryEndIndex;
    		this.targetEndIndex = targetEndIndex;
    		this.isAnchored = isAnchored;
    	}
    	/**
    	 * Convert a list of anchors into a subproblem list.
    	 * @param anchors anchored read pairs
    	 * @param querySequenceLength length of query sequence
    	 * @param targetSequenceLength length of target sequence
    	 * @return list alignment subproblems
    	 */
    	public static List<Subproblem> getSubproblems(List<Anchor> anchors, int querySequenceLength, int targetSequenceLength) {
    		Collections.sort(anchors, new Anchor.QueryIndexComparator());
    		List<Subproblem> list = new ArrayList<Subproblem>();
    		Anchor last = new Anchor(-1, -1); // sentinal anchor
    		boolean isAnchored = false;
    		for (int i = 0; i < anchors.size(); i++) {
    			if (anchors.get(i).targetIndex <= last.targetIndex ||
					anchors.get(i).queryIndex <= last.queryIndex) {
    				throw new IllegalArgumentException("Anchor set must allow at least one possible alignment.");
    			}
    			list.add(new Subproblem(
    					last.queryIndex + 1,
						last.targetIndex + 1,
						anchors.get(i).queryIndex,
						anchors.get(i).targetIndex,
						isAnchored));
    			last = anchors.get(i);
    			isAnchored = true;
    		}
    		list.add(new Subproblem(
    				last.queryIndex + 1,
					last.targetIndex + 1,
					querySequenceLength,
					targetSequenceLength,
					isAnchored));
    		return list;
    	}
    }
    // updates cut rows given the latest row of traceback pointers
    public static void setCuts(int x, Subproblem subproblem, Last[][] pointers, Cut[]cuts) {
        for (Cut c : cuts) {
            c.update(x, subproblem, pointers);
        }
    }
    /**
     * Calculate the optimal alignment score for the given sequence positions with an affine or constant gap penalty
     * @param x position in query
     * @param y position in target
     * @param gop gap opening penalty
     * @param gep gap extension penalty
     * @param sub compound match score
     * @param scores dynamic programming score matrix to fill at the given position
     * @return traceback direction for substitution, deletion and insertion
     */
    public static Last[] setScorePoint(int x, int y, int gop, int gep, int sub, int[][][] scores) {
        Last[] pointers = new Last[3];

        // substitution
        if (scores[x - 1][y - 1][1] >= scores[x - 1][y - 1][0] && scores[x - 1][y - 1][1] >= scores[x - 1][y - 1][2]) {
            scores[x][y][0] = scores[x - 1][y - 1][1] + sub;
            pointers[0] = Last.DELETION;
        } else if (scores[x - 1][y - 1][0] >= scores[x - 1][y - 1][2]) {
            scores[x][y][0] = scores[x - 1][y - 1][0] + sub;
            pointers[0] = Last.SUBSTITUTION;
        } else {
            scores[x][y][0] = scores[x - 1][y - 1][2] + sub;
            pointers[0] = Last.INSERTION;
        }

        // deletion
        if (scores[x - 1][y][1] >= scores[x - 1][y][0] + gop) {
            scores[x][y][1] = scores[x - 1][y][1] + gep;
            pointers[1] = Last.DELETION;
        } else {
            scores[x][y][1] = scores[x - 1][y][0] + gop + gep;
            pointers[1] = Last.SUBSTITUTION;
        }

        // insertion
        if (scores[x][y - 1][0] + gop >= scores[x][y - 1][2]) {
            scores[x][y][2] = scores[x][y - 1][0] + gop + gep;
            pointers[2] = Last.SUBSTITUTION;
        } else {
            scores[x][y][2] = scores[x][y - 1][2] + gep;
            pointers[2] = Last.INSERTION;
        }

        return pointers;
    }
    /**
     * Calculates the optimal alignment score for the given sequence positions and a linear gap penalty
     * @param x position in query
     * @param y position in target
     * @param gep gap extension penalty
     * @param sub compound match score
     * @param scores dynamic programming score matrix to fill at the given position
     * @return traceback directions for substitution, deletion and insertion respectively
     */
    public static Last setScorePoint(int x, int y, int gep, int sub, int[][][] scores) {
        int d = scores[x - 1][y][0] + gep;
        int i = scores[x][y - 1][0] + gep;
        int s = scores[x - 1][y - 1][0] + sub;
        if (d >= s && d >= i) {
            scores[x][y][0] = (int) d;
            return Last.DELETION;
        } else if (s >= i) {
            scores[x][y][0] = (int) s;
            return Last.SUBSTITUTION;
        } else {
            scores[x][y][0] = (int) i;
            return Last.INSERTION;
        }
    }

    /**
     * Score global alignment for a given position in the query sequence
     * @param x
     * @param subproblem
     * @param gop
     * @param gep
     * @param subs
     * @param storing
     * @param scores
     * @return
     */
    public static Last[][] setScoreVector(int x, Subproblem subproblem, int gop, int gep, int[] subs, boolean storing,
            int[][][] scores) {
        return setScoreVector(x, subproblem.getQueryStartIndex(), subproblem.getTargetStartIndex(), subproblem.getTargetEndIndex(), gop, gep, subs, storing, scores, subproblem.isStartAnchored());
    }

    /**
     * Score global alignment for a given position in the query sequence
     * @param x
     * @param xb
     * @param yb
     * @param ye
     * @param gop
     * @param gep
     * @param subs
     * @param storing
     * @param scores
     * @param startAnchored
     * @return
     */
    public static Last[][] setScoreVector(int x, int xb, int yb, int ye, int gop, int gep, int[] subs,
            boolean storing, int[][][] scores, boolean startAnchored) {
        Last[][] pointers = new Last[ye + 1][];
        int min = Integer.MIN_VALUE - gop - gep;
        ensureScoringMatrixColumn(x, storing, scores);
        if (x == xb) {
            scores[xb][yb][1] = scores[xb][yb][2] = gop;
            pointers[yb] = new Last[] {null, null, null};
            if (startAnchored) {
        		assert (xb > 0 && yb > 0);
        		int subproblemStartingScore = scores[xb - 1][yb - 1][0] + subs[yb];
        		scores[xb][yb][0] = subproblemStartingScore;
        		scores[xb][yb][1] = subproblemStartingScore + gop;
        		scores[xb][yb][2] = subproblemStartingScore + gop;
        		pointers[yb] = new Last[] {Last.SUBSTITUTION, Last.SUBSTITUTION, Last.SUBSTITUTION};
        	}
            Last[] insertion = new Last[] { null, null, Last.INSERTION };
            for (int y = yb + 1; y <= ye; y++) {
                scores[xb][y][0] = scores[xb][y][1] = min;
                scores[xb][y][2] = (int) (scores[xb][y - 1][2] + gep);
                pointers[y] = insertion;
            }
        } else {
            scores[x][yb][0] = scores[x][yb][2] = min;
            scores[x][yb][1] = scores[x - 1][yb][1] + gep;
            pointers[yb] = new Last[] { null, Last.DELETION, null };
            for (int y = yb + 1; y <= ye; y++) {
                pointers[y] = setScorePoint(x, y, gop, gep, subs[y], scores);
            }
        }
        return pointers;
    }

    /**
     * Score global alignment for a given position in the query sequence for a linear gap penalty
     * @param x
     * @param subproblem
     * @param gep
     * @param subs
     * @param storing
     * @param scores
     * @return
     */
    public static Last[][] setScoreVector(int x, Subproblem subproblem, int gep, int[] subs, boolean storing,
            int[][][] scores) {
        return setScoreVector(x, subproblem.getQueryStartIndex(), subproblem.getTargetStartIndex(), subproblem.getTargetEndIndex(), gep, subs, storing, scores, subproblem.isStartAnchored());
    }

    /**
     * Score global alignment for a given position in the query sequence for a linear gap penalty
     * @param x
     * @param xb
     * @param yb
     * @param ye
     * @param gep
     * @param subs
     * @param storing
     * @param scores
     * @param startAnchored
     * @return
     */
    public static Last[][] setScoreVector(int x, int xb, int yb, int ye, int gep, int[] subs, boolean storing,
            int[][][] scores, boolean startAnchored) {
        Last[][] pointers = new Last[ye + 1][1];
        ensureScoringMatrixColumn(x, storing, scores);
        if (x == xb) {
        	if (startAnchored) {
        		assert (xb > 0 && yb > 0);
        		scores[xb][yb][0] = (int) (scores[xb - 1][yb - 1][0] + subs[yb]);
        		pointers[yb][0] = Last.SUBSTITUTION;
        	}
            for (int y = yb + 1; y <= ye; y++) {
                scores[xb][y][0] = (int) (scores[xb][y - 1][0] + gep);
                pointers[y][0] = Last.INSERTION;
            }
        } else {
            scores[x][yb][0] = (int) (scores[x - 1][yb][0] + gep);
            pointers[yb][0] = Last.DELETION;
            for (int y = yb + 1; y <= ye; y++) {
                pointers[y][0] = setScorePoint(x, y, gep, subs[y], scores);
            }
        }
        return pointers;
    }

    /**
     * Score local alignment for a given position in the query sequence
     * @param x
     * @param gop
     * @param gep
     * @param subs
     * @param storing
     * @param scores
     * @param xyMax
     * @param score
     * @return
     */
    public static Last[][] setScoreVector(int x, int gop, int gep, int[] subs, boolean storing,
            int[][][] scores, int[] xyMax, int score) {
        return setScoreVector(x, 0, 0, scores[0].length - 1, gop, gep, subs, storing, scores, xyMax, score);
    }

    /**
     * Score local alignment for a given position in the query sequence
     * @param x
     * @param xb
     * @param yb
     * @param ye
     * @param gop
     * @param gep
     * @param subs
     * @param storing
     * @param scores
     * @param xyMax
     * @param score
     * @return
     */
    public static Last[][] setScoreVector(int x, int xb, int yb, int ye, int gop, int gep, int[] subs,
            boolean storing, int[][][] scores, int[] xyMax, int score) {
        Last[][] pointers;
        ensureScoringMatrixColumn(x, storing, scores);
        if (x == xb) {
            pointers = new Last[ye + 1][scores[0][0].length];
        } else {
            pointers = new Last[ye + 1][];
            pointers[0] = new Last[scores[0][0].length];   
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
    
    /**
     * Score local alignment for a given position in the query sequence for a linear gap penalty
     * @param x
     * @param gep
     * @param subs
     * @param storing
     * @param scores
     * @param xyMax
     * @param score
     * @return
     */
    public static Last[][] setScoreVector(int x, int gep, int[] subs, boolean storing, int[][][] scores,
            int[] xyMax, int score) {
        return setScoreVector(x, 0, 0, scores[0].length - 1, gep, subs, storing, scores, xyMax, score);
    }

    /**
     * Score local alignment for a given position in the query sequence for a linear gap penalty
     * @param x
     * @param xb
     * @param yb
     * @param ye
     * @param gep
     * @param subs
     * @param storing
     * @param scores
     * @param xyMax
     * @param score
     * @return
     */
    public static Last[][] setScoreVector(int x, int xb, int yb, int ye, int gep, int[] subs, boolean storing,
            int[][][] scores, int[] xyMax, int score) {
        Last[][] pointers;
        ensureScoringMatrixColumn(x, storing, scores);
        if (x == xb) {
            pointers = new Last[ye + 1][1];
        } else {
            pointers = new Last[ye + 1][];
            pointers[0] = new Last[1];
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
    
    private static void ensureScoringMatrixColumn(int x, boolean storingFullMatrix, int[][][] scores) {
    	if (!storingFullMatrix && x > 1) {
            scores[x] = scores[x - 2];
        }
    }
    
    /**
     * Find alignment path through traceback matrix
     * @param traceback
     * @param local
     * @param xyMax
     * @param last
     * @param sx
     * @param sy
     * @return
     */
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

    /**
     * Find global alignment path through traceback matrix
     * @param traceback
     * @param scores
     * @param sx
     * @param sy
     * @return
     */
    public static int[] setSteps(Last[][][] traceback, int[][][] scores, List<Step> sx, List<Step> sy) {
        int xMax = scores.length - 1, yMax = scores[xMax].length - 1;
        boolean linear = (traceback[xMax][yMax].length == 1);
        
        Last last = 
        		
        	linear ? 
        		traceback[xMax][yMax][0] : 
        			
        		(scores[xMax][yMax][1] > scores[xMax][yMax][0] &&
        		 scores[xMax][yMax][1] > scores[xMax][yMax][2] ) ? 
        						
        				Last.DELETION : 
        					(scores[xMax][yMax][0] > scores[xMax][yMax][2]) ? 
        							Last.SUBSTITUTION : 
        							Last.INSERTION;
        
        		
        return setSteps(traceback, false, new int[] {xMax, yMax}, last, sx, sy);
    }

    /**
     * Find local alignment path through traceback matrix
     * @param traceback
     * @param xyMax
     * @param sx
     * @param sy
     * @return
     */
    public static int[] setSteps(Last[][][] traceback, int[] xyMax, List<Step> sx, List<Step> sy) {
        return setSteps(traceback, true, xyMax, Last.SUBSTITUTION, sx, sy);
    }
    
    public static String tracebackToString(Last[][][] traceback) {
    	StringBuilder sb = new StringBuilder();
    	for (int z = 0; z < 3; z++) {
    		for (int i = 0; i < traceback.length; i++) {
    			if (traceback[i] != null) {
	    			for (int j = 0; j < traceback[i].length; j++) {
	    				if (traceback[i][j] == null || z >= traceback[i][j].length || traceback[i][j][z] == null) {
	    					sb.append('.');
	    				} else {
		    				switch (traceback[i][j][z]) {
		    	            case DELETION:
		    	            	sb.append('^');
		    	            	break;
		    	            case SUBSTITUTION:
		    	            	sb.append('\\');
		    	                break;
		    	            case INSERTION:
		    	            	sb.append('<');
		    	            	break;
		    	            }
	    				}
	    			}
    			}
    			sb.append('\n');
    		}
    		sb.append("\n\n");
    	}
    	return sb.toString();
    }
}

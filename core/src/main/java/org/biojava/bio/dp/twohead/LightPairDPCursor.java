/*
 * BioJava development code
 * 
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 * 
 * http://www.gnu.org/copyleft/lesser.html
 * 
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 * 
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 * 
 * http://www.biojava.org
 *
 */

package org.biojava.bio.dp.twohead;

import java.util.Arrays;
import java.util.List;

import org.biojava.bio.dp.BackPointer;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;

/**
 * A LIGHT implementation of PairDPCursor.
 * <p>
 * This object manages memory that is linear on the length of the shortest
 * sequence. It does not maintain any data beyond that necessary for the next
 * round of calcCell invocations.
 *
 * @author     Matthew Pocock 
 * @author     David Huen (fixes for magical state)
 */
public class LightPairDPCursor implements PairDPCursor {
  //    protected CrossProductAlphabet alpha;
  private int[] pos;
  private boolean flip;
  private SymbolList[] seqs;
  private double[][][] columns;
  private double[][][] emissions;
  /**
   *  Description of the Field 
   */
  protected BackPointer[][][] bPointers;
  /**
   *  Description of the Field 
   */
  protected int numStates;
  /**
   *  Description of the Field 
   */
  protected double[] zeroCol;
  /**
   *  Description of the Field 
   */
  protected BackPointer[] emptyBP;
  /**
   *  Description of the Field 
   */
  protected int[] depth;
  /**
   *  Description of the Field 
   */
  protected EmissionCache eCache;


  /**
   *  Constructor for the LightPairDPCursor object 
   *
   * @param  seq1                        First sequence in this twohead DP.
   * @param  seq2                        Second sequence in this twohead DP.
   * @param  depth1                      The number of bases of context required in 
   *                                     first sequence to compute DP matrix at a point (= max advance in first sequence + 1).
   * @param  depth2                      The number of bases of context required in 
   *                                     second sequence to compute DP matrix at a point (= max advance in second sequence + 1).
   * @param  numStates                   Total number of states in model.
   * @param  eCache                      Emission cache to be used with this run.
   * @exception  IllegalSymbolException  Description of Exception 
   */
  public LightPairDPCursor(
      SymbolList seq1, 
      SymbolList seq2, 
  //      CrossProductAlphabet alpha,
      int depth1, 
      int depth2, 
      int numStates, 
      EmissionCache eCache
  ) throws IllegalSymbolException {
    this.numStates = numStates;
    //      this.alpha = alpha;
    this.zeroCol = new double[this.numStates];
    // don't touch this, please...
    for (int i = 0; i < zeroCol.length; ++i) {
      this.zeroCol[i] = Double.NaN;
    }
    this.emptyBP = new BackPointer[numStates];
    this.pos = new int[2];
    this.pos[0] = 0;
    this.pos[1] = 0;
    this.seqs = new SymbolList[2];
    this.seqs[0] = seq1;
    this.seqs[1] = seq2;
    this.depth = new int[2];
    this.depth[0] = depth1;
    this.depth[1] = depth2;
    this.eCache = eCache;

    this.flip = this.seqs[1].length() > this.seqs[0].length();
    //System.out.println("flip=" + this.flip);
    if (flip) {
      this.columns = 
          new double[depth[0]][seqs[1].length() + 2][numStates];
      this.emissions = 
          new double[depth[0]][seqs[1].length() + 2][];
      this.bPointers = 
          new BackPointer[depth[0]][seqs[1].length() + 2][numStates];
    }
    else {
      this.columns = 
          new double[depth[1]][seqs[0].length() + 2][numStates];
      this.emissions = 
          new double[depth[1]][seqs[0].length() + 2][];
      this.bPointers = 
          new BackPointer[depth[1]][seqs[0].length() + 2][numStates];
    }

    for (int i = 0; i < columns.length; i++) {
      double[][] ci = columns[i];
      for (int j = 0; j < ci.length; j++) {
        double[] cj = ci[j];
        for (int k = 0; k < cj.length; k++) {
          cj[k] = Double.NaN;
        }
      }
    }

    // compute the emissions vector at origin of matrix
    calcEmissions(emissions[0]);
  }


  /**
   *  Gets the Depth attribute of the LightPairDPCursor object 
   *
   * @return    The Depth value 
   */
  public int[] getDepth() {
    return depth;
  }


  /**
   *  Are there further Cells to be computed? 
   *
   * @return    Description of the Returned Value 
   */
  public boolean hasNext() {
    int i = flip ? 0 : 1;
    return 
        pos[i] <= seqs[i].length() + 1;
  }


  /**
   * Returns the minimal context of the DP matrix
   * necessary to compute the value of a single point
   * in that matrix.
   * <p>
   * The Cell [][] array has the origin as the
   * current point to be evaluated and successive
   * rows/columns represent rows/columns backwards
   * within the DP matrix. 
   *
   * @return    An array representing the immediate context around the element to be computed. 
   */
  public Cell[][] press() {
    Cell[][] cells = new Cell[depth[0]][depth[1]];
    for (int i = 0; i < cells.length; i++) {
      Cell[] ci = cells[i];
      for (int j = 0; j < ci.length; j++) {
        ci[j] = new Cell();
      }
    }
    return cells;
  }


  /**
   *  Description of the Method 
   *
   * @param  cells                       Description of Parameter 
   * @exception  IllegalSymbolException  Description of Exception 
   */
  public void next(Cell[][] cells) throws IllegalSymbolException {
    //System.out.println("Pos=" + pos[0] + ", " + pos[1] + " " + flip);
    if (flip) {
      for (int i = 0; i < depth[0]; i++) {
        int ii = pos[0] - i;
        boolean outI = (ii < 0) || (ii > seqs[0].length() + 1);
        Cell[] cellI = cells[i];
        // lucky in this case - can pre-fetch cells
        double[][] columnsI = columns[i];
        double[][] emissionsI = emissions[i];
        BackPointer[][] bPointersI = bPointers[i];
        for (int j = 0; j < depth[1]; j++) {
          int jj = pos[1] - j;
          boolean outJ = (jj < 0) || (jj > seqs[1].length() + 1);
          //System.out.println("at " + i + "->" + ii + ", " + j + "->" + jj);
          Cell c = cellI[j];
          if (outI || outJ) {
            c.scores = zeroCol;
            c.emissions = zeroCol;
            c.backPointers = emptyBP;
          }
          else {
            c.scores = columnsI[jj];
            c.emissions = emissionsI[jj];
            c.backPointers = bPointersI[jj];
          }
        }
      }
      if (pos[1] <= seqs[1].length()) {
        pos[1]++;
      }
      else {
        pos[1] = 0;
        pos[0]++;

        // advance arrays
        int depth = this.depth[0];
        double[][] tempC = columns[depth - 1];
        double[][] tempE = emissions[depth - 1];
        BackPointer[][] tempBP = bPointers[depth - 1];
        for (int i = 1; i < depth; i++) {
          columns[i] = columns[i - 1];
          emissions[i] = emissions[i - 1];
          bPointers[i] = bPointers[i - 1];
        }
        columns[0] = tempC;
        emissions[0] = tempE;
        bPointers[0] = tempBP;

        calcEmissions(tempE);
      }
    }
    else {
      // flip == false
      for (int i = 0; i < depth[1]; i++) {
        int ii = pos[1] - i;
        boolean outI = (ii < 0) || (ii > seqs[1].length() + 1);
        //Cell[] cellI = cells[i]; // can't pre-fetch as loop wrong way
        double[][] columnsI = columns[i];
        double[][] emissionsI = emissions[i];
        BackPointer[][] bPointersI = bPointers[i];
        for (int j = 0; j < depth[0]; j++) {
          int jj = pos[0] - j;
          boolean outJ = (jj < 0) || (jj > seqs[0].length() + 1);
          //System.out.println("at " + i + "->" + ii + ", " + j + "->" + jj);
          Cell c = cells[j][i];
          if (outI || outJ) {
            c.scores = zeroCol;
            c.emissions = zeroCol;
            c.backPointers = emptyBP;
          }
          else {
            c.scores = columnsI[jj];
            c.emissions = emissionsI[jj];
            c.backPointers = bPointersI[jj];
          }
        }
      }
      if (pos[0] <= seqs[0].length()) {
        pos[0]++;
      }
      else {
        pos[0] = 0;
        pos[1]++;

        // advance arrays
        int depth = this.depth[1];
        double[][] tempC = columns[depth - 1];
        double[][] tempE = emissions[depth - 1];
        BackPointer[][] tempBP = bPointers[depth - 1];
        for (int i = 1; i < depth; i++) {
          columns[i] = columns[i - 1];
          emissions[i] = emissions[i - 1];
          bPointers[i] = bPointers[i - 1];
        }
        columns[0] = tempC;
        emissions[0] = tempE;
        bPointers[0] = tempBP;

        calcEmissions(tempE);
      }
    }
  }


  /**
   *  Description of the Method 
   *
   * @param  em                          Description of Parameter 
   * @exception  IllegalSymbolException  Description of Exception 
   */
  private void calcEmissions(double[][] em)
       throws IllegalSymbolException {
    if (flip) {
      //System.out.println("Calculating emissions at " + pos[0] + ", " + pos[1]);
      Symbol[] symL = new Symbol[2];
      List symList = Arrays.asList(symL);
      int i = pos[0];
      symL[0] = (i < 1 || i > seqs[0].length())
           ? AlphabetManager.getGapSymbol()
           : seqs[0].symbolAt(i);
      for (int j = 0; j <= seqs[1].length() + 1; j++) {
        symL[1] = (j < 1 || j > seqs[1].length())
             ? AlphabetManager.getGapSymbol()
             : seqs[1].symbolAt(j);
        em[j] = eCache.getEmissions(symList, !(( i < 1 && j < 1 ) || ((i > seqs[0].length()) && (j > seqs[1].length())) ) );
//        System.out.println("symbol " + symL[0].getName() + ", " + symL[1].getName() + "->" + em[j] + " " + (( i < 1 && j < 1 ) || ((i > seqs[0].length()) && (j > seqs[1].length())) ) );
      }
    }
    else {
      //System.out.println("Calculating emissions at " + pos[0] + ", " + pos[1]);
      Symbol[] symL = new Symbol[2];
      List symList = Arrays.asList(symL);
      int j = pos[1];
      symL[1] = (j < 1 || j > seqs[1].length())
           ? AlphabetManager.getGapSymbol()
           : seqs[1].symbolAt(j);
      for (int i = 0; i <= seqs[0].length() + 1; i++) {
        symL[0] = (i < 1 || i > seqs[0].length())
             ? AlphabetManager.getGapSymbol()
             : seqs[0].symbolAt(i);
//        System.out.println("symbol " + symL[0].getName() + ", " + symL[1].getName() + " " + (( i < 1 && j < 1 ) || ((i > seqs[0].length()) && (j > seqs[1].length())) ) );
        em[i] = eCache.getEmissions(symList, !(( i < 1 && j < 1 ) || ((i > seqs[0].length()) && (j > seqs[1].length())) ) );
      }
    }
  }
}


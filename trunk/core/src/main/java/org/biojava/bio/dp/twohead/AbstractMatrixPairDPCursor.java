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
 * @author Matthew Pocock
 */
public abstract class AbstractMatrixPairDPCursor
  implements PairDPCursor {
    protected int[] pos;
    protected SymbolList[] seqs;
    protected double[][][] columns;
    protected BackPointer[][][] bPointers;
    protected double[][][] emissions;
    protected int numStates;
    protected double[] zeroCol;
    protected BackPointer[] emptyBP;
    protected int[] depth;
    protected double[][][] sMatrix;
    protected EmissionCache eCache;

    public AbstractMatrixPairDPCursor(
      SymbolList seq1,
      SymbolList seq2,
      int start1,
      int start2,
      int depth1,
      int depth2,
      PairDPMatrix matrix,
      EmissionCache eCache
    ) throws IllegalSymbolException {
      this.numStates = matrix.states().length;

      this.zeroCol = new double[this.numStates]; // don't touch this, please...
      for (int i = 0; i < zeroCol.length; ++i) {
        this.zeroCol[i] = Double.NaN;
      }
      this.emptyBP = new BackPointer[numStates];
      
      this.sMatrix = matrix.getScoreArray();

      this.pos = new int[2];
      this.pos[0] = start1;
      this.pos[1] = start2;
      this.seqs = new SymbolList[2];
      this.seqs[0] = seq1;
      this.seqs[1] = seq2;
      this.depth = new int[2];
      this.depth[0] = depth1;
      this.depth[1] = depth2;
      this.bPointers = new BackPointer[seq1.length()+2][seq2.length()+2][numStates];
      this.emissions = new double[seq1.length()+2][seq2.length()+2][];
      this.eCache = eCache;
      
      Symbol [] symArray = new Symbol[2];
      List symList = Arrays.asList(symArray);
      for(int i = 0; i <= seq1.length()+1; i++) {
        symArray[0] = (i < 1 || i > seq1.length())
          ? AlphabetManager.getGapSymbol()
          : seq1.symbolAt(i);
        double [][] ei = emissions[i];
        for(int j = 0; j <= seq2.length()+1; j++) {
          symArray[1] = (j < 1 || j > seq2.length())
            ? AlphabetManager.getGapSymbol()
            : seq2.symbolAt(j);
          ei[j] = eCache.getEmissions(symList, !((i < 1 && j < 1) || (i > seq1.length() && j <= seq2.length())) );
        }
      }
    }
    
    public int [] getDepth() {
      return depth;
    }

    public Cell[][] press() {
      Cell [][] cells = new Cell[depth[0]][depth[1]];
      for(int i = 0; i < cells.length; i++) {
        Cell [] ci = cells[i];
        for(int j = 0; j < ci.length; j++) {
          ci[j] = new Cell();
        }
      }
      return cells;
    }
  }

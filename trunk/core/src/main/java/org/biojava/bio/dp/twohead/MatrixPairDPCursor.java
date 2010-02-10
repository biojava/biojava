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

import org.biojava.bio.dp.BackPointer;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.SymbolList;

  /**
   * @author Matthew Pocock
   */
  public class MatrixPairDPCursor
  extends AbstractMatrixPairDPCursor {
    public MatrixPairDPCursor(
      SymbolList seq1,
      SymbolList seq2,
      int depth1,
      int depth2,
      PairDPMatrix matrix,
      EmissionCache eCache
    ) throws IllegalSymbolException {
      super(seq1, seq2, 0, 0, depth1, depth2, matrix, eCache);
    }
    
    public boolean hasNext() {
      return
        pos[1] <= (seqs[1].length()+1);
    }
    
    public void next(Cell[][] cells) {
      for(int i = 0; i < depth[0]; i++) {
        Cell[] cellI = cells[i];
        int ii = pos[0] - i;
        boolean outI = ii < 0 || ii > seqs[0].length()+1;
        if(outI) {
          for(int j = 0; j < depth[1]; j++) {
            Cell c = cellI[j];
            c.scores = zeroCol;
            c.emissions = zeroCol;
            c.backPointers = emptyBP;
          }
        } else {
          double[][] sMatI = this.sMatrix[ii];
          double[][] emisI = this.emissions[ii];
          BackPointer[][] bPI = this.bPointers[ii];
          for(int j = 0; j < depth[1]; j++) {
            int jj = pos[1] - j;
            boolean outJ = jj < 0 || jj > seqs[1].length()+1;
            Cell c = cellI[j];
            if(outJ) {
              c.scores = zeroCol;
              c.emissions = zeroCol;
              c.backPointers = emptyBP;
            } else {
              c.scores = sMatI[jj];
              c.emissions = emisI[jj];
              c.backPointers = bPI[jj];
            }
          }
        }
      }
      if(pos[0] <= seqs[0].length()) {
        pos[0]++;
      } else {
        pos[0] = 0;
        pos[1]++;
      }
    }
  }


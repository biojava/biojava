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
 */

package org.biojava.bio.dp.twohead;

import org.biojava.bio.symbol.IllegalSymbolException;

/**
 * A cursor over a DP matrix.
 *
 * @author Matthew Pocock
 */
public interface PairDPCursor {
  /** test wether the cursor can be advanced further */
  boolean hasNext();
  /** retrieve the next block of cells */
  void next(Cell [][] cells) throws IllegalSymbolException;
  /** press out a new correctly sized cell array */
  Cell [][] press();
  /** retrieve the depth of this cursor */
  int [] getDepth();
}

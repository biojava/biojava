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

import org.biojava.bio.dp.BackPointer;

/**
 * <p>
 * A single cell in the DP matrix.
 * </p>
 *
 * <p>
 * The cell is divided into parallel arrays for the cell scores, backpointers
 * and emission probabilities. The scores and backpointer arrays are as long
 * as the number of states in the model. The emissions array is as long as
 * the number of emitting states in the model.
 * </p>
 *
 * @author Matthew Pocock
 * @since 1.2
 */
public final class Cell {
  public double [] scores;
  public BackPointer [] backPointers;
  public double [] emissions;
}

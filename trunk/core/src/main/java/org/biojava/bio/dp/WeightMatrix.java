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


package org.biojava.bio.dp;

import org.biojava.bio.dist.Distribution;
import org.biojava.bio.symbol.Alphabet;

/**
 * A log odds weight matrix.
 * <p>
 * The weight matrix uses computer-coordinates. Thus, a 10 column weight matrix
 * has columns (0 - 9). I guess that if you try to access columns outside the
 * logical range, the implementation may throw an IndexOutOfBoundsException.
 *
 * @author Matthew Pocock
 */
public interface WeightMatrix {
  /**
   * The alphabet for the sequences that this weight matrix models.
   *
   * @return  the Alphabet
   */
  Alphabet getAlphabet();
  
  /**
   * The number of columns modeled by the weight matrix.
   *
   * @return the number of columns
   */
  int columns();
  
  /**
   * Retrieve a column as an EmissionState.
   * <p>
   * To find the emission probability for Symbol sym at column col use:
   * <code>wm.getColumn(col).getWeight(sym)</code>.
   *
   * @param column  the weight matrix column to retrieve
   * @throws IndexOutOfBoundsException if column is not between 0 and
   *         columns()-1
   * @return the EmissionState that represents the individual column
   */
  Distribution getColumn(int column) throws IndexOutOfBoundsException;
}

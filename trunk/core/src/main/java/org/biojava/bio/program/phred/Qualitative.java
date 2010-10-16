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

package org.biojava.bio.program.phred;

import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;

/**
 * <p>Qualitative is an interface for classes wanting to hold quality
 * data in symbolic form such as Phred scores.<p>
 *
 * <p>Copyright (c) 2001</p>
 * <p>Company:      AgResearch</p>
 *
 * @author Mark Schreiber
 * @since 1.1
 */

public interface Qualitative {

  /**
   * Retreives the list of quality symbols from the underlying object.
   */
  SymbolList getQuality();

  /**
   * Retreives the quality symbol for the specified index.
   * @param index - Must be greater than zero.
   * @throws IndexOutOfBoundsException if index is outside of the quality symbol list.
   */
  Symbol getQualityAt(int index) throws IndexOutOfBoundsException;
}

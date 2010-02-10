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
import org.biojava.bio.dp.IllegalTransitionException;
import org.biojava.bio.dp.ScoreType;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;

/**
 * @author Matthew Pocock
 */
public interface CellCalculatorFactory {
  CellCalculator forwards(ScoreType scoreType)
  throws IllegalSymbolException, IllegalAlphabetException, IllegalTransitionException;

  CellCalculator backwards(ScoreType scoreType)
  throws IllegalSymbolException, IllegalAlphabetException, IllegalTransitionException;

  CellCalculator viterbi(ScoreType scoreType, BackPointer terminal)
  throws IllegalSymbolException, IllegalAlphabetException, IllegalTransitionException;
}

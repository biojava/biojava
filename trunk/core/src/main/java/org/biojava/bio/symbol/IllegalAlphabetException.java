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


package org.biojava.bio.symbol;

import org.biojava.bio.BioException;

/**
 * <p>
 * The exception to indicate that an invalid alphabet has been used.
 * </p>
 *
 * <p>
 * The usual reason for throwing an IllegalAlphabetException is that you are
 * trying to parse a SymbolList into a method that only works for some
 * alphabets, but not for the alphabet associated with that SymbolList.
 * </p>
 *
 * @author Matthew Pocock
 */
public class IllegalAlphabetException extends BioException {
  /**
   * Just make the exception.
   */
  public IllegalAlphabetException() { super(); }

  /**
   * Make the exception with a message.
   */
  public IllegalAlphabetException(String message) { super(message); }

  public IllegalAlphabetException(Throwable t) { super(t); }

  public IllegalAlphabetException(Throwable t, String message) { super( message, t); }
}

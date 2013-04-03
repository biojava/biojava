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
 * Created on DATE
 *
 */
package org.biojava3.core.exceptions;

/**
 * Thrown from AbstractCompundTranslator
 * @author Andy Yates
 *
 */

public class TranslationException extends RuntimeException {

  private static final long serialVersionUID = -3017433758219757440L;

  public TranslationException(String m) {
    super(m);
  }

  public TranslationException(Exception t) {
    super(t);
  }

  public TranslationException(String m, Exception t) {
    super(m, t);
  }

}

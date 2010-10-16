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


/*
 * BioException.java
 *
 * Coppied from AceError By Thomas Down <td2@sanger.ac.uk>
 */

package org.biojava.bio;

/**
 * A nestable biological error.
 *

 * Catch this whenever it, or one of it's sub-classes are thrown and you know
 * what to do once you've got it. Note: in general, you should not be catching
 * errors. However, there are cases where it is necisary e.g. for logging. You
 * will nearly always want to either re-throw the Error, throw a new Error
 * or exit the current thread.
 *
 * 
 * Throw this when something has gone wrong and in general people should not be
 * handeling it.
 * @author Matthew Pocock
 * @since 1.0
 */
public class BioError extends Error {
  /**
   * Create a new BioError with a message.
   *
   * @param message  the message
   */
  public BioError(String message) {
	  super(message);
  }

  /**
   * Create a new BioError with a cause.
   *
   * @param ex  the Throwable that caused this BioError
   */
  public BioError(Throwable ex) {
    super(ex);
  }

  /**
   * Create a new BioError with a cause and a message.
   *
   * @param ex  the Throwable that caused this BioError
   * @param message  the message
   * @deprecated Use BioError(message, ex) instead.
   */
  public BioError(Throwable ex, String message) {
    this(message, ex);
  }

  /**
   * Create a new BioError with a cause and a message.
   *
   * @param message  the message
   * @param ex  the Throwable that caused this BioError
   */
  public BioError(String message, Throwable ex) {
    super(message, ex);
  }
  
  /**
   * Create a new BioError.
   */
  public BioError() {
	  super();
  }
}

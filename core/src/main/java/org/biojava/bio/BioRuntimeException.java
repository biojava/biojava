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
package org.biojava.bio;

/**
 * A nestable biological exception.
 *
 * <p>
 * In BioJava, checked exceptions are generally preferred to RuntimeExceptions,
 * but RuntimeExceptions can be used as a fall-back if you are implementing
 * an interface which doesn't support checked exceptions.  If you do this,
 * please document this clearly in the implementing class.
 * </p>
 * 
 * Occasionaly methods will document that the object they return may emit
 * BioRuntimeExceptions. If you know what you are doing, you can catch these
 * and do something sensible with them. An example would be an iterator() method
 * that returns an Iterator where next() could fail due to a database connection
 * flaking out. This could raise a BioRuntimeException wrapping the real cause.
 * It would then not be unreasonable for the code calling next() to catch this
 * BioRuntimeException and extract the cause.
 *
 * 
 * Occasionaly it is necisary to abuse the exception system by throwing BioError
 * to get arround limitations in APIs that you do not control. For example, when
 * implementing Iterator.next(), you may need to call a method that can fail.
 * Catch the failure and throw it as a BioError. If you want people to handle
 * this gracefully in their scripts, then document this behavior and hope they
 * catch the error.
 *
 * @author Matthew Pocock
 * @author Thomas Down
 *
 * 
 * 
 *
 * @since 1.0
 */

public class BioRuntimeException extends RuntimeException {
  /**
   * Create a new BioRuntimeException with a message.
   *
   * @param message  the message
   */
  public BioRuntimeException(String message) {
	  super(message);
  }


  /**
   * Create a new BioRuntimeException with a cause.
   *
   * @param ex  the Throwable that caused this BioRuntimeException
   */
  public BioRuntimeException(Throwable ex) {
    super(ex);
  }

  /**
   * Create a new BioRuntimeException with a cause and a message.
   *
   * @param ex  the Throwable that caused this BioRuntimeException
   * @param message  the message
   * @deprecated use new BioRuntimeException(message, ex) instead
   */
  public BioRuntimeException(Throwable ex, String message) {
    this(message, ex);
  }
  
  /**
   * Create a new BioRuntimeException with a cause and a message.
   *
   * @param message  the message
   * @param ex  the Throwable that caused this BioRuntimeException
   */
  public BioRuntimeException(String message, Throwable ex) {
    super(message, ex);
  }
  
  /**
   * Create a new BioRuntimeException.
   */
  public BioRuntimeException() {
	  super();
  }
}

package org.biojava.bibliography;

/**
 * An exception raised when communciation with the BibRef APIs fails.
 *
 * @author Matthew Pocock
 * @since 1.4
 */
public class BibRefException
extends Exception {
  /**
   * Create a new BibRefException with a message.
   *
   * @param message  the message of the exception
   */
  public BibRefException(String message) {
    super(message);
  }

  /**
   * Create a new BibRefException with a root cause.
   *
   * @param cause  the unerlying cause of this exception
   */
  public BibRefException(Throwable cause) {
    super(cause);
  }

  /**
   * Create a nw BibRefException with a message and a root cause.
   *
   * @param message   the message for the exception
   * @param cause     the underlying cuase of this exception
   */
  public BibRefException(String message, Throwable cause) {
    super(message, cause);
  }
}

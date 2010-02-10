package org.biojava.bio.program.ssaha;

/**
 * There has been some failure that prevents a search from completing.
 *
 * @author Matthew Pocock
 * @since 1.4
 */
public class SearchException
extends Exception {
  public SearchException(String message) {
    super(message);
  }
  
  public SearchException(Throwable cause) {
    super(cause);
  }
  
  public SearchException(String message, Throwable cause) {
    super(message, cause);
  }
}

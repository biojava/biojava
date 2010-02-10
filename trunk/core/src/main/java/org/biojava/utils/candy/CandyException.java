package org.biojava.utils.candy;

public class CandyException
extends Exception {
  public CandyException(String message) {
    super(message);
  }
  
  public CandyException(Throwable cause) {
    super(cause);
  }
  
  public CandyException(String message, Throwable cause) {
    super(message, cause);
  }
}

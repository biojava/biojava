package org.biojava.utils;

public class CommitFailure
extends Exception {
  public CommitFailure(String message) {
    super(message);
  }
  
  public CommitFailure(Throwable cause) {
    super(cause);
  }
  
  public CommitFailure(String message, Throwable cause) {
    super(message, cause);
  }
}

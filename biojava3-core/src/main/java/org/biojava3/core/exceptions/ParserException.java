package org.biojava3.core.exceptions;

public class ParserException extends RuntimeException {

  private static final long serialVersionUID = -4101924035353204493L;

  public ParserException(String message) {
    super(message);
  }

  public ParserException(Exception e) {
    super(e);
  }

  public ParserException(String message, Exception e) {
    super(message, e);
  }

}

package org.biojava3.core.exceptions;

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

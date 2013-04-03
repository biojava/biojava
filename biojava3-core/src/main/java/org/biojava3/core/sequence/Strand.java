package org.biojava3.core.sequence;

/**
 * Provides a way of representing the strand of a sequence or location.
 *
 * @author ayates
 */
public enum Strand {

  POSITIVE("+", 1), NEGATIVE("-", -1), UNKNOWN("?", 0);

  private final String stringRepresentation;
  private final int    numericRepresentation;

  private Strand(String stringRepresentation, int numericRepresentation) {
    this.stringRepresentation = stringRepresentation;
    this.numericRepresentation = numericRepresentation;
  }

  public int getNumericRepresentation() {
    return numericRepresentation;
  }

  public String getStringRepresentation() {
    return stringRepresentation;
  }
}

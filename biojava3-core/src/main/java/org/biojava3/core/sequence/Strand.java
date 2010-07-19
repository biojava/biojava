package org.biojava3.core.sequence;

/**
 * Provides a way of representing the strand of a sequence, location
 * hit or feature.
 *
 * @author ayates
 */
public enum Strand {

    POSITIVE("+", 1), NEGATIVE("-", -1), UNDEFINED(".", 0);
    private final String stringRepresentation;
    private final int numericRepresentation;

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

    public Strand getReverse() {
        switch (this) {
            case POSITIVE:
                return NEGATIVE;
            case NEGATIVE:
                return POSITIVE;
            default:
                return UNDEFINED;
        }
    }
}

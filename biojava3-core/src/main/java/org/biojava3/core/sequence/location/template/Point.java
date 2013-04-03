package org.biojava3.core.sequence.location.template;

/**
 * Holds a single point part of a location
 */
public interface Point extends Comparable<Point> {

    /**
     * Used to resolve a position about a point
     */
    public interface Resolver<T extends Point> {
        int resolve(T point);
    }

    /**
     * Returns the position held by this object
     */
    Integer getPosition();

    /**
     * Returns true if the current position is unknown but is
     * beyond the position encoded for. This is the same as the position
     * <pre>&gt;80</pre> as encoded by UniProt.
     */
    boolean isUnknown();

    /**
     * Returns a true if the exact point is unknown. Equivalent position
     * from UniProt is <pre>?80</pre>.
     */
    boolean isUncertain();

    /**
     * Returns the equivalent position on the reverse strand
     *
     * @param length Length of the sequence to translate to
     */
    Point reverse(int length);

    /**
     * Returns a new point offset by the given distance
     */
    Point offset(int distance);

    /**
     * Returns true if the current point is at a lower position than the
     * point given.
     */
    boolean isLower(Point point);

    /**
     * Returns true if the point is higher in value to the current point
     */
    boolean isHigher(Point point);

    /**
     * Returns a copy of this point
     */
    Point clonePoint();
}

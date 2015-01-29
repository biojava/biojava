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
 * Created on 01-21-2010
 */
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

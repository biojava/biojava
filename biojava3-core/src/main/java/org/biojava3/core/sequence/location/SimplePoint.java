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
package org.biojava3.core.sequence.location;

import org.biojava3.core.sequence.location.template.Point;
import org.biojava3.core.util.Hashcoder;
import org.biojava3.core.util.Equals;

/**
 * Basic implementation of the Point interface.
 *
 * @author ayates
 */
public class SimplePoint implements Point {

    private int position;
    private boolean unknown;
    private boolean uncertain;

    protected SimplePoint() {
        super();
    }

    public SimplePoint(int position) {
        this.position = position;
    }

    public SimplePoint(int position, boolean unknown, boolean uncertain) {
        this.position = position;
        this.unknown = unknown;
        this.uncertain = uncertain;
    }

    @Override
    public Integer getPosition() {
        return position;
    }

    protected void setPosition(int position) {
        this.position = position;
    }

    @Override
    public boolean isUnknown() {
        return unknown;
    }

    protected void setUnknown(boolean unknown) {
        this.unknown = unknown;
    }

    @Override
    public boolean isUncertain() {
        return uncertain;
    }

    protected void setUncertain(boolean uncertain) {
        this.uncertain = uncertain;
    }

    @Override
    public Point reverse(int length) {
        int translatedPosition = reverse(getPosition(), length);
        return new SimplePoint(translatedPosition, isUnknown(), isUncertain());
    }

    @Override
    public Point offset(int distance) {
        int offsetPosition = getPosition() + distance;
        return new SimplePoint(offsetPosition, isUnknown(), isUncertain());
    }

    protected int reverse(int position, int length) {
        return (length - position) + 1;
    }

    @Override
    @SuppressWarnings("EqualsWhichDoesntCheckParameterClass")
    public boolean equals(Object obj) {
        boolean equals = false;
        if (Equals.classEqual(this, obj)) {
            SimplePoint p = (SimplePoint) obj;
            equals = (Equals.equal(getPosition(), p.getPosition())
                    && Equals.equal(isUncertain(), p.isUncertain())
                    && Equals.equal(isUnknown(), p.isUnknown()));
        }
        return equals;
    }

    @Override
    public int hashCode() {
        int r = Hashcoder.SEED;
        r = Hashcoder.hash(r, getPosition());
        r = Hashcoder.hash(r, isUncertain());
        r = Hashcoder.hash(r, isUnknown());
        return r;
    }

    @Override
    public String toString() {
        return Integer.toString(getPosition());
    }

    @Override
    public int compareTo(Point o) {
        return getPosition().compareTo(o.getPosition());
    }

    @Override
    public boolean isLower(Point point) {
        return (compareTo(point) < 0);
    }

    @Override
    public boolean isHigher(Point point) {
        return (compareTo(point) > 0);
    }

    @Override
    public Point clonePoint() {
        return this.offset(0);
    }
}

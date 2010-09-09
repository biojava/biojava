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

    
    public Integer getPosition() {
        return position;
    }

    protected void setPosition(int position) {
        this.position = position;
    }

    
    public boolean isUnknown() {
        return unknown;
    }

    protected void setUnknown(boolean unknown) {
        this.unknown = unknown;
    }

    
    public boolean isUncertain() {
        return uncertain;
    }

    protected void setUncertain(boolean uncertain) {
        this.uncertain = uncertain;
    }

    
    public Point reverse(int length) {
        int translatedPosition = reverse(getPosition(), length);
        return new SimplePoint(translatedPosition, isUnknown(), isUncertain());
    }

    
    public Point offset(int distance) {
        int offsetPosition = getPosition() + distance;
        return new SimplePoint(offsetPosition, isUnknown(), isUncertain());
    }

    protected int reverse(int position, int length) {
        return (length - position) + 1;
    }

    
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

    
    public int hashCode() {
        int r = Hashcoder.SEED;
        r = Hashcoder.hash(r, getPosition());
        r = Hashcoder.hash(r, isUncertain());
        r = Hashcoder.hash(r, isUnknown());
        return r;
    }

    
    public String toString() {
        return Integer.toString(getPosition());
    }

    
    public int compareTo(Point o) {
        return getPosition().compareTo(o.getPosition());
    }

    
    public boolean isLower(Point point) {
        return (compareTo(point) < 0);
    }

    
    public boolean isHigher(Point point) {
        return (compareTo(point) > 0);
    }
}

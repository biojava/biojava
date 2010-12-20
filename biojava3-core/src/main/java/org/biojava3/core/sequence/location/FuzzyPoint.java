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
import org.biojava3.core.util.Equals;
import org.biojava3.core.util.Hashcoder;

/**
 * Implementation for resolving fuzzy locations. Caches the calculated
 * value.
 *
 * @author ayates
 */
public class FuzzyPoint extends SimplePoint {

    /**
     * Always uses the min point to resolve a location
     */
    public static final Resolver<FuzzyPoint> MIN_RESOLVER = new Resolver<FuzzyPoint>() {
        @Override
        public int resolve(FuzzyPoint point) {
            return point.getMin();
        }
    };

    /**
     * Always uses the max point to resolve a location
     */
    public static final Resolver<FuzzyPoint> MAX_RESOLVER = new Resolver<FuzzyPoint>() {
        @Override
        public int resolve(FuzzyPoint point) {
            return point.getMax();
        }
    };

    /**
     * Combines min and max and then gets the mean of it
     */
    public static final Resolver<FuzzyPoint> MEAN_RESOLVER = new Resolver<FuzzyPoint>() {
        @Override
        public int resolve(FuzzyPoint point) {
            return (point.getMin() + point.getMax()) / 2;
        }
    };

    private final int min;
    private final int max;
    private final Resolver<FuzzyPoint> resolver;

    public FuzzyPoint(int minPoint, int maxPoint) {
        this(minPoint, maxPoint, MEAN_RESOLVER, false, false);
    }

    public FuzzyPoint(int minPoint, int maxPoint, Resolver<FuzzyPoint> resolver) {
        this(minPoint, maxPoint, resolver, false, false);
    }

    public FuzzyPoint(int minPoint, int maxPoint, Resolver<FuzzyPoint> resolver, boolean unknown, boolean uncertain) {
        this.min = minPoint;
        this.max = maxPoint;
        this.resolver = resolver;
        setUncertain(uncertain);
        setUnknown(unknown);
        setPosition(-1); //Means we have not resolved this position yet
    }

    @Override
    public Integer getPosition() {
        if(super.getPosition() == -1) {
            super.setPosition(getResolver().resolve(this));
        }
        return super.getPosition();
    }

    protected Integer getMax() {
        return max;
    }

    protected Integer getMin() {
        return min;
    }

    protected Resolver<FuzzyPoint> getResolver() {
        return resolver;
    }

    @Override
    public Point reverse(int length) {
        int revMin = reverse(getMin(), length);
        int revMax = reverse(getMax(), length);
        return new FuzzyPoint(revMin, revMax, getResolver(), isUnknown(), isUncertain());
    }

    @Override
    public Point offset(int distance) {
        int offMin = getMin() + distance;
        int offMax = getMax() + distance;
        return new FuzzyPoint(offMin, offMax, getResolver(), isUnknown(), isUncertain());
    }


    @Override
    @SuppressWarnings("EqualsWhichDoesntCheckParameterClass")
    public boolean equals(Object obj) {
        boolean equals = false;
        if (Equals.classEqual(this, obj)) {
            FuzzyPoint p = (FuzzyPoint) obj;
            equals = (Equals.equal(getMin(), p.getMin())
                    && Equals.equal(getMax(), p.getMax())
                    && Equals.equal(isUnknown(), p.isUnknown())
                    && Equals.equal(isUncertain(), p.isUncertain())
                    );
        }
        return equals;
    }

    @Override
    public int hashCode() {
        int r = Hashcoder.SEED;
        r = Hashcoder.hash(r, getMin());
        r = Hashcoder.hash(r, getMax());
        r = Hashcoder.hash(r, isUncertain());
        r = Hashcoder.hash(r, isUnknown());
        return r;
    }

    @Override
    public int compareTo(Point point) {
        //If we can assign this to a FuzzyPoint then work with a bit more info
        if(FuzzyPoint.class.isAssignableFrom(point.getClass())) {
            FuzzyPoint fuzzy = (FuzzyPoint)point;
            int minComparison = getMin().compareTo(fuzzy.getMin());
            if(minComparison != 0)
                return minComparison;
            return getMax().compareTo(fuzzy.getMax());
        }
        //If not fuzzy then compare on position as normal
        return super.compareTo(point);
    }
}

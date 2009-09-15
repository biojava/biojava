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
 */

package org.biojava.bio.symbol;

import java.io.Serializable;

/**
 * A 'fuzzy' location a-la Embl fuzzy locations.
 * <p>
 * Fuzzy locations have propreties that indicate that they may start before min
 * and end after max. However, this in no way affects how they interact with
 * other locations.
 * </p>
 *
 * @author Matthew Pocock
 * @author Thomas Down
 * @author Greg Cox
 */
public class FuzzyLocation
extends AbstractRangeLocation
implements Serializable {
    /**
     * Always use the `inner' values.
     */

    public final static RangeResolver RESOLVE_INNER;

    /**
     * Use the `outer' values, unless they are unbounded in which case the
     * `inner' values are used.
     */

    public final static RangeResolver RESOLVE_OUTER;

    /**
     * Use the arithmetic mean of the `inner' and `outer' values, unless the
     * outer value is unbounded.
     */

    public final static RangeResolver RESOLVE_AVERAGE;

    static {
        RESOLVE_INNER = new InnerRangeResolver();
        RESOLVE_OUTER = new OuterRangeResolver();
        RESOLVE_AVERAGE = new AverageRangeResolver();
    }

    private int outerMin;
    private int innerMin;
    private int innerMax;
    private int outerMax;
    private boolean mIsMinFuzzy;
    private boolean mIsMaxFuzzy;
    private RangeResolver resolver;

  /**
   * Create a new FuzzyLocation with endpoints (outerMin.innerMin) and (innerMax.outerMax).
   *
   * @param outerMin the lower bound on the location's min value.
   *				 Integer.MIN_VALUE indicates unbounded.
   * @param outerMax the upper bound on the location's max value.
   * 				 Integer.MAX_VALUE indicates unbounded.
   * @param innerMin the upper bound on the location's min value.
   * @param innerMax the lower bound on the location's max value.
   * @param resolver a RangeResolver object which defines the policy used to calculate
   *                 the location's min and max properties.
   */

  public FuzzyLocation(
    int outerMin, int outerMax,
    int innerMin, int innerMax,
    RangeResolver resolver
  ) {
        boolean isMinFuzzy = false;
        boolean isMaxFuzzy = false;
        if (outerMin != innerMin)
        {
            isMinFuzzy = true;
        }
        if (outerMax != innerMax)
        {
            isMaxFuzzy = true;
        }
        this.initializeVariables(outerMin, outerMax, innerMin, innerMax, isMinFuzzy, isMaxFuzzy, resolver);
  }

    /**
     * Create a new FuzzyLocation with endpoints (outerMin.innerMin) and
     * (innerMax.outerMax).  This constructor allows you to explicitly mark an
     * endpoint as fuzzy, even if there is no other information about it.  For
     * example, a valid swissprot location "?5 10" would be a fuzzy location 5
     * to 10 where the min is fuzzy and the max is not.
     * <p>
     * Note that it is not logical to specify inner and outer values that
     * clearly denote fuzzy boundaries and the set the <code>isMinFuzzy</code> or
     * <code>isMaxFuzzy</code> value to false. This object makes
     * no specific check of your logic so be careful.
     *
     * @param outerMin the lower bound on the location's min value.
     *				 Integer.MIN_VALUE indicates unbounded.
     * @param outerMax the upper bound on the location's max value.
     * 				 Integer.MAX_VALUE indicates unbounded.
     * @param innerMin the upper bound on the location's min value.
     * @param innerMax the lower bound on the location's max value.
     * @param isMinFuzzy Explictly state if the minimum is fuzzy
     * @param isMaxFuzzy Explictly state if the maximum is fuzzy
     * @param resolver a RangeResolver object which defines the policy used to
     * 				   calculate the location's min and max properties.
     */

    public FuzzyLocation(int outerMin, int outerMax,
                         int innerMin, int innerMax,
                         boolean isMinFuzzy, boolean isMaxFuzzy,
                         RangeResolver resolver)
    {
        this.initializeVariables(outerMin, outerMax, innerMin, innerMax, isMinFuzzy, isMaxFuzzy, resolver);
    }

  public Location translate(int dist) {
    return new FuzzyLocation(
      outerMin + dist,
      outerMax + dist,
      innerMin + dist,
      innerMax + dist,
      resolver
    );
  }

  /**
   * Retrieve the Location that this decorates.
   *
   * @return the Location instance that stores all of the Location interface
   *         data
   */

  public RangeResolver getResolver() {
    return resolver;
  }

  public int getOuterMin() {
    return outerMin;
  }


  public int getOuterMax() {
    return outerMax;
  }

  public int getInnerMin() {
    return innerMin;
  }


  public int getInnerMax() {
    return innerMax;
  }

  public int getMin() {
    return resolver.resolveMin(this);
  }

  public int getMax() {
    return resolver.resolveMax(this);
  }

  public boolean hasBoundedMin() {
    return outerMin != Integer.MIN_VALUE;
  }

  public boolean hasBoundedMax() {
    return outerMax != Integer.MAX_VALUE;
  }

    public String toString()
    {
        return "["
            + (hasBoundedMin() ? Integer.toString(getMin()) : "<" + Integer.toString(getMin()))
            + ", "
            + (hasBoundedMax() ? Integer.toString(getMax()) : ">" + Integer.toString(getMax()))
            + "]";
    }

    /**
     * Determines how a <code>FuzzyLocation</code> should be treated when used
     * as a normal <code>Location</code>.
     *
     * Use one of the implementations of this interface when creating a <code>FuzzyLocation</code>
     * to specify how the fuzzy (inner/outer) properties are translated into the standard
     * Location min and max properties.
     *
     * It is possible to write custom implementations of this to create <code>FuzzyLocations</code>
     * with exotic behaviour.
     */

    public static interface RangeResolver extends Serializable {
        /**
         * Delegate for the getMin() method.
         * @param loc The Location to resolve
         * @return the resolved Min
         */

        public int resolveMin(FuzzyLocation loc);

        /**
         * Delegate for the getMax() method.
         * @param loc The Location to resolve
         * @return the resolved Max
         */

        public int resolveMax(FuzzyLocation loc);
    }

    private static class InnerRangeResolver implements RangeResolver {
        public int resolveMin(FuzzyLocation loc) {
            return loc.getInnerMin();
        }

        public int resolveMax(FuzzyLocation loc) {
            return loc.getInnerMax();
        }
    }

    private static class OuterRangeResolver implements RangeResolver {
        public int resolveMin(FuzzyLocation loc) {
            if (loc.hasBoundedMin()) {
                return loc.getOuterMin();
            } else {
                return loc.getInnerMin();
            }
        }

        public int resolveMax(FuzzyLocation loc) {
            if (loc.hasBoundedMax()) {
                return loc.getOuterMax();
            } else {
                return loc.getInnerMax();
            }
        }
    }

    private static class AverageRangeResolver implements RangeResolver {
        public int resolveMin(FuzzyLocation loc) {
            if (loc.hasBoundedMin()) {
                return (loc.getOuterMin() + loc.getInnerMin()) / 2;
            } else {
                return loc.getInnerMin();
            }
        }

        public int resolveMax(FuzzyLocation loc) {
            if (loc.hasBoundedMax()) {
                return (loc.getOuterMax() + loc.getInnerMax()) / 2;
            } else {
                return loc.getInnerMax();
            }
        }
    }

    public boolean isMinFuzzy()
    {
        return mIsMinFuzzy;
    }

    public boolean isMaxFuzzy()
    {
        return mIsMaxFuzzy;
    }

    /**
     * Refactored initialization code from the constructors.
     *
     * @param outerMin the lower bound on the location's min value.
     *				 Integer.MIN_VALUE indicates unbounded.
     * @param outerMax the upper bound on the location's max value.
     * 				 Integer.MAX_VALUE indicates unbounded.
     * @param innerMin the upper bound on the location's min value.
     * @param innerMax the lower bound on the location's max value.
     * @param resolver a RangeResolver object which defines the policy used to calculate
     *                 the location's min and max properties.
     */
    protected void initializeVariables(int outerMin, int outerMax,
                                       int innerMin, int innerMax,
                                       boolean isMinFuzzy, boolean isMaxFuzzy,
                                       RangeResolver resolver)
    {
        this.outerMin = outerMin;
        this.outerMax = outerMax;
        this.innerMin = innerMin;
        this.innerMax = innerMax;
        this.resolver = resolver;
        this.mIsMinFuzzy = isMinFuzzy;
        this.mIsMaxFuzzy = isMaxFuzzy;
    }
}

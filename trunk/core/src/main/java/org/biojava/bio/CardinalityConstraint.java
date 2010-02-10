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
package org.biojava.bio;

import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.RangeLocation;

/**
 * A constraint on the number of values a property can have.
 *
 * @author Matthew Pocock
 * @since 1.3
 *
 * Usefull constants for whenever you need one of the common
 * cardinalitites. Otherwise, build a Location using the normal Location
 * APIs.:
 */
public final class CardinalityConstraint {
  /**
   * This cardinality contains no intengers, not even zero. It means that there
   * is no way to fulfill this cardinality constraint. It's like Double.NaN
   */
  public static final Location NONE
    = Location.empty;
  /**
   * The property should have zero values. This means that it should be absent.
   */
  public static final Location ZERO
    = new RangeLocation(0, 0);
  /**
   * The property should have zero or one values. This means that it is optional
   * but if present must have exactly one value.
   */
  public static final Location ZERO_OR_ONE
    = new RangeLocation(0, 1);
  /**
   * The property can have any number of values, including none.
   */
  public static final Location ANY
    = new RangeLocation(0, Integer.MAX_VALUE);
  /**
   * The property should have exactly one value.
   */
  public static final Location ONE
    = new RangeLocation(1, 1);
  /**
   * The property should have one or more values. It can not be absent.
   */
  public static final Location ONE_OR_MORE
    = new RangeLocation(1, Integer.MAX_VALUE);
  
  private CardinalityConstraint() {}
}

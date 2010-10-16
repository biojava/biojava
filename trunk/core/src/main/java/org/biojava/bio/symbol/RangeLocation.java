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
 * A simple implementation of Location that contains all points between
 * getMin and getMax inclusive.
 * <p>
 * This will in practice be the most commonly used pure-java implementation.
 *
 * @author Matthew Pocock
 */
public class RangeLocation
extends AbstractRangeLocation
implements Serializable {
  /**
   * The minimum point contained.
   */
  private int min;

  /**
   * The maximum point contained.
   */
  private int max;

  public int getMin() {
    return min;
  }

  public int getMax() {
    return max;
  }

  /**
   * Construct a new RangeLocation from <code>min</code> to <code>max</code>.
   */
  
  public RangeLocation(int min, int max) throws IndexOutOfBoundsException {
    if(max < min) {
      throw new IndexOutOfBoundsException(
        "max must exceed min: min=" + min + ", max=" + max
      );
    }
    this.min = min;
    this.max = max;
  }

  public Location translate(int dist) {
    return new RangeLocation(min + dist, max + dist);
  }

  public String toString() {
    return "[" + getMin() + "," + getMax() + "]";
  }
}

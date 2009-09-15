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
 * A location representing a single point.  This can be considered as
 * the singleton set of one integer.
 * <p>
 * min and max are always equal for this implementation
 * </p>
 *
 * @author Matthew Pocock
 */
public class PointLocation
extends AbstractRangeLocation
implements Location, Serializable {
  /**
   * The actual index contained.
   */
  private int point;

  public int getMin()	{ return point; }
  public int getMax()	{ return point; }
  public boolean contains(int p)	{ return this.point == p; }
  
  public Location translate(int dist) {
      if (dist == 0)
	  return this;
      return new PointLocation(this.point + dist);
  }
  
  public PointLocation(int point) {
    this.point = point;
  }
  
  public String toString() {
    return String.valueOf(point);
  }
}

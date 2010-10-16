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
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import org.biojava.bio.BioError;

/**
 * <p>
 * A complex location. It is made up from multiple sub-locations and is essentially
 * the point-wise union of the child locations.
 * </p>
 *
 * <p>
 * <strong>NOTE:</strong> It is no longer possible to directly construct
 * CompoundLocations.  Use LocationTools.union instead.
 * </p>
 *
 * @author Matthew Pocock
 * @author Thomas Down
 * @author Greg Cox
 */

class CompoundLocation extends AbstractLocation
      implements Location, Serializable 
{
  /**
   * The list of child locations in order. Should contain only RangeLocation
   * instances.
   */
  private List locations;

  /**
   * Minimum index contained.
   */
  private int min = Integer.MAX_VALUE;

  /**
   * Maximum index contained.
   */
  private int max = Integer.MIN_VALUE;

  /**
   * <p>
   * Generate a new CompoundLocation from a list of locations.
   * </p>
   *
   * <p>
   * You will nearly always want to generate these beasts using the
   * Location.union method.
   * </p>
   *
   * <p>
   * The locations list should contain contiguous locations, sorted by getMin()
   * and guaranteed to be non-overlapping.
   * </p>
   *
   * @param locations a list of Location instances to combine into a single
   *        compound location
   */
  CompoundLocation(List locs) {
    this.locations = new ArrayList(locs);
    Collections.sort(locations, Location.naturalOrder);
    locations = Collections.unmodifiableList(locations);

    Location minL = (Location) locations.get(0);
    Location maxL = (Location) locations.get(locs.size() - 1);
    this.min = minL.getMin();
    this.max = maxL.getMax();
  }


  public int getMin() {
    return min;
  }

  public int getMax() {
    return max;
  }

  public boolean contains(int p) {
    if(p < min || p > max)
      return false;

    int m = locations.size();
    int n = 0;
    while (m >= n) {
	int i = (m + n) / 2;
	Location l = (Location) locations.get(i);

	if (p < l.getMin()) {
	    m = i - 1;
	} else if (p > l.getMax()) {
	    n = i + 1;
	} else {
	    return true;
	}
    }

    return false;
  }

  public boolean contains(Location l) {
    return LocationTools.contains(this, l);
  }

  public boolean overlaps(Location l) {
    return LocationTools.overlaps(this, l);
  }

  public Location union(Location loc) {
    return LocationTools.union(this, loc);
  }

  public Location intersection(Location loc) {
    return LocationTools.intersection(this, loc);
  }

  public boolean equals(Object o) {
    if(!(o instanceof Location)) {
      return false;
    } else {
      return LocationTools.areEqual(this, (Location) o);
    }
  }

  public SymbolList symbols(SymbolList seq) {
      if (isContiguous())
	  return seq.subList(min, max);

    List res = new ArrayList();
    for (Iterator i = blockIterator(); i.hasNext(); ) {
	Location l = (Location) i.next();
	res.addAll(l.symbols(seq).toList());
    }

    try {
      return new SimpleSymbolList(seq.getAlphabet(), res);
    } catch (IllegalSymbolException ex) {
      throw new BioError(ex);
    }
  }

  public Location translate(int dist) {
      if (dist == 0) {
	  return this;
      }

    List res = new ArrayList();

    for(Iterator i = locations.iterator(); i.hasNext(); ) {
      res.add( ((Location) i.next()).translate(dist) );
    }

    return new CompoundLocation(res);
  }

  public boolean isContiguous() {
    return locations.size() <= 1;
  }

    List getBlockList() {
	return locations;
    }

  public Iterator blockIterator() {
    return locations.iterator();
  }

  public int hashCode() {
    return getMin() ^ getMax();
  }

  public String toString() {
    StringBuffer sb = new StringBuffer();
    sb.append(getMin() + ", " + getMax() + " {");
    Iterator i = locations.iterator();
    if(i.hasNext())
      sb.append("(" + i.next() + ")");
    while(i.hasNext())
      sb.append(", (" + i.next() + ")");
    sb.append("}");

    return sb.substring(0);
  }
}

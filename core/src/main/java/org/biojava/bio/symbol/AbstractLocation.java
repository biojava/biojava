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
 * An abstract implementation of <code>Location</code>.
 *
 * This provides implementations of the binary operators which delegate to
 *                the <code>LocationTools</code> class.
 * @author Matthew Pocock
 */
public abstract class AbstractLocation
implements Location, Serializable {
  public Location getDecorator(Class decoratorClass) {
    if(decoratorClass.isInstance(this)) {
      return this;
    } else {
      return null;
    }
  }
  
  public Location newInstance(Location loc) {
    return loc;
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
  
  public int hashCode() {
    return getMin() ^ getMax();
  }
}

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

import java.io.NotSerializableException;
import java.io.ObjectStreamException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;

import org.biojava.bio.BioError;
import org.biojava.utils.StaticMemberPlaceHolder;

/**
 * Package-private location iml for empty locations.
 *
 * Location class should ensure that only one of these exist.
 *
 * @author Matthew Pocock
 * @author Thomas Down
 */
class EmptyLocation implements Location, Serializable {
    public Location getDecorator(Class decoratorClass) {
      if(decoratorClass.isInstance(this)) {
        return this;
      } else {
        return null;
      }
    }
    public Location newInstance(Location loc) { return loc; }
    public int getMin() { return Integer.MAX_VALUE; }
    public int getMax() { return Integer.MIN_VALUE; }
    public boolean overlaps(Location l) { return false; }
    public boolean contains(Location l) { return false; }
    public boolean contains(int p) { return false; }
    public boolean equals(Object o) {
      if(o instanceof Location) {
        if(o instanceof EmptyLocation) {
          return true;
        } else {
          return Location.naturalOrder.areEqual(this, (Location) o);
        }
      } else {
        return false;
      }
    }
    public Location intersection(Location l) { return empty; }
    public Location union(Location l) { return l; }
    public SymbolList symbols(SymbolList seq) {
      try {
        return new SimpleSymbolList(seq.getAlphabet(), new ArrayList());
      } catch (IllegalSymbolException ex) {
        throw new BioError(ex);
      }
    }
    public Location translate(int dist) { return this; }
    public boolean isContiguous() { return true; }
    public Iterator blockIterator() { return Collections.EMPTY_SET.iterator(); }
    private Object writeReplace() throws ObjectStreamException {
      try {
        return new StaticMemberPlaceHolder(Location.class.getField("empty"));
      } catch (NoSuchFieldException nsfe) {
        throw new NotSerializableException(nsfe.getMessage());
      }
    }
    public String toString() {
      return "{}";
    }
}


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
import java.util.Comparator;
import java.util.Iterator;

import org.biojava.utils.StaticMemberPlaceHolder;

/**
 * Comparator class for locations.
 *
 * Location should ensure that only one instance exists.
 *
 * @author Matthew Pocock
 * @author Thomas Down
 */
class LocationComparator implements Comparator, Serializable {
    public int compare(Object o1, Object o2) {
      int d = 0;

      Location l1 = (Location) o1;
      Location l2 = (Location) o2;

      Iterator i1 = l1.blockIterator();
      Iterator i2 = l2.blockIterator();

      while(i1.hasNext() && i2.hasNext()) {
        Location li1 = (Location) i1.next();
        Location li2 = (Location) i2.next();

        d = li1.getMin() - li2.getMin();
        if(d != 0) {
          return d;
        }
        d = li1.getMax() - li2.getMax();
        if(d != 0) {
          return d;
        }
      }
      if(i2.hasNext()) {
        return 1;
      } else if(i1.hasNext()) {
        return -1;
      }

      return 0;
    }

    public boolean equals(Object obj) {
      return obj == this;
    }

    /**
    *Test whether two locations are equal or not
    */
    public boolean areEqual(Location l1, Location l2) {
      Iterator i1 = l1.blockIterator();
      Iterator i2 = l2.blockIterator();

      while(i1.hasNext() && i2.hasNext()) {
        if(! i1.next().equals(i2.next()) ) {
          return false;
        }
      }

      if(!i1.hasNext() && !i2.hasNext()) {
        return false;
      }

      return true;
    }

    private Object writeReplace() throws ObjectStreamException {
      try {
        return new StaticMemberPlaceHolder(Location.class.getField("naturalOrder"));
      } catch (NoSuchFieldException nsfe) {
        throw new NotSerializableException(nsfe.getMessage());
      }
    }
}

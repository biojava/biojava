package org.biojava3.core.sequence.location.template;

import static org.biojava3.core.util.Equals.classEqual;
import static org.biojava3.core.util.Equals.equal;

import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Stack;

import org.biojava3.core.sequence.Strand;
import org.biojava3.core.util.Hashcoder;

public abstract class AbstractLocation implements Location {

  private final int            min;
  private final int            max;
  private final Strand         strand;
  private final List<Location> subLocations;

  public AbstractLocation(int min, int max, Strand strand) {
    this.min = min;
    this.max = max;
    this.strand = strand;
    this.subLocations = Collections.emptyList();
  }

  public AbstractLocation(int min, int max, Strand strand,
      List<Location> subLocations) {
    this.min = min;
    this.max = max;
    this.strand = strand;
    this.subLocations = Collections.unmodifiableList(subLocations);
  }

  public int getMax() {
    return max;
  }

  public int getMin() {
    return min;
  }

  public Strand getStrand() {
    return strand;
  }

  public List<Location> getSubLocations() {
    return subLocations;
  }

  public boolean isComplex() {
    return getSubLocations().size() > 0;
  }

  public Iterator<Location> iterator() {
    if(!isComplex()) {
      //Basic iterator
      return new Iterator<Location>() {
        private int index = 0;
        public boolean hasNext() {
          return index++ < 1;
        }
        public Location next() {
          return AbstractLocation.this;
        }
        public void remove() {
          throw new UnsupportedOperationException("Cannot remove from a Location");
        };
      };
    }
    else {
      //more complex one which decends down each location in turn
      return new Iterator<Location>() {
        private Stack<Iterator<Location>> iteratorStack = new Stack<Iterator<Location>>();
        private Iterator<Location> currentIterator = AbstractLocation.this.getSubLocations().iterator();
        {
          iteratorStack.push(currentIterator);
        }

        public boolean hasNext() {
          boolean hasNext = currentIterator.hasNext();
          return hasNext;
        }

        public Location next() {
          Location nextLocation = currentIterator.next();
          if(nextLocation == null) {
          }
          else if(nextLocation.isComplex()) {
            //TODO How do we descend down from here?!? Need to jump into the sub-locations
          }
          return nextLocation;
        }

        public void remove() {
          throw new UnsupportedOperationException("Cannot remove from a Location");
        }
      };
    }
  }

  public boolean equals(Object obj) {
    boolean equals = false;
    if (classEqual(this, obj)) {
      AbstractLocation l = (AbstractLocation) obj;
      equals = (
          equal(getMin(), l.getMin()) &&
          equal(getMax(), l.getMax()) &&
          equal(getStrand(), l.getStrand()) &&
          equal(getSubLocations(), l.getSubLocations())
      );
    }
    return equals;
  }

  public int hashCode() {
    int r = Hashcoder.SEED;
    r = Hashcoder.hash(r, getMin());
    r = Hashcoder.hash(r, getMax());
    r = Hashcoder.hash(r, getStrand());
    r = Hashcoder.hash(r, getSubLocations());
    return r;
  }
}

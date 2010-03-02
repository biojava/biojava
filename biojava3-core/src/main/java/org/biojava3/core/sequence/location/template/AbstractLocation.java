package org.biojava3.core.sequence.location.template;

import static org.biojava3.core.util.Equals.classEqual;
import static org.biojava3.core.util.Equals.equal;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import org.biojava3.core.sequence.Strand;
import org.biojava3.core.util.Hashcoder;

/**
 * Basic implementation of a location which has a start, end and a location.
 *
 * @author ayates
 */
public abstract class AbstractLocation implements Location {

  private final int             start;
  private final int             end;
  private final Strand          strand;
  private final List<Location>  subLocations;

  public AbstractLocation(int start, int end, Strand strand) {
    this.start = start;
    this.end = end;
    this.strand = strand;
    this.subLocations = Collections.emptyList();
  }

  public AbstractLocation(int start, int end, Strand strand,
      List<Location> subLocations) {
    this.start = start;
    this.end = end;
    this.strand = strand;
    this.subLocations = Collections.unmodifiableList(subLocations);
  }

  public AbstractLocation(int start, int end, Strand strand,
      Location... subLocations) {
    this(start, end, strand, Arrays.asList(subLocations));
  }

  public int getEnd() {
    return end;
  }

  public int getStart() {
    return start;
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
    List<Location> locations;
    if(isComplex()) {
      locations = new ArrayList<Location>(getSubLocations());
    }
    else {
      locations = new ArrayList<Location>(Arrays.asList(AbstractLocation.this));
    }
    return locations.iterator();
  }

  public List<Location> getAllSubLocations() {
    return getAllSubLocations(this);
  }

  private List<Location> getAllSubLocations(Location location) {
    List<Location> subLocations = new ArrayList<Location>();
    for(Location l: location) {
      subLocations.add(l);
      if(l.isComplex()) {
        subLocations.addAll(getAllSubLocations(l));
      }
    }
    return subLocations;
  }

  public boolean equals(Object obj) {
    boolean equals = false;
    if (classEqual(this, obj)) {
      AbstractLocation l = (AbstractLocation) obj;
      equals = (
          equal(getStart(), l.getStart()) &&
          equal(getEnd(), l.getEnd()) &&
          equal(getStrand(), l.getStrand()) &&
          equal(getSubLocations(), l.getSubLocations())
      );
    }
    return equals;
  }

  public int hashCode() {
    int r = Hashcoder.SEED;
    r = Hashcoder.hash(r, getStart());
    r = Hashcoder.hash(r, getEnd());
    r = Hashcoder.hash(r, getStrand());
    r = Hashcoder.hash(r, getSubLocations());
    return r;
  }

  /**
   * Always returns false
   */
  public boolean isCircular() {
    return false;
  }

  public String toString() {
    return getStart()+":"+getEnd()+"("+getStrand().getStringRepresentation()+")";
  }
}

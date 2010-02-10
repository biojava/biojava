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

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.BioRuntimeException;

/**
 * <p>
 * A class to hide the messiness of CircularLocs from LocationTools
 * </p>
 *
 * <p>
 * <b>WARNING</b> All binary operations on CircularLocations are currently assumed to
 * be performed on Locations from the same sequence, or sequences of the same
 * length. To enforce this binary operations will currently only accept circular
 * locations of the same length.
 * </p>
 *
 * @author Mark Schreiber
 * @author Greg Cox
 * @author Thomas Down
 * @version 1.0
 */

final class CircularLocationTools {
  /**
   * Translates coordinates from circular into linear where base 1 is the first
   * base and base 0 is the base before 1 or the last base of the sequence base
   * -1 is the penultimate base of the sequence etc.
   *
   * @param val The circular coordinate
   * @param length The length of the circular molecule
   * @return the transformed coordinate
   */
  private static int realValue(int val, int length){
    val = (val % length);
    if(val < 1) val = length + val;
    return val;
  }

  /**
   * Makes a circular location, called by LocationTools.
   *
   * @param min  minimum index of location
   * @param max  maximum index of Location
   * @param seqLength  used to work out if the feature wraps arround the origin
   */
  protected static CircularLocation makeCircLoc(int min, int max, int seqLength){
      if(min == max)
        return new CircularLocation(new PointLocation(min), seqLength);

      boolean overlap = false;
      int rmin = realValue(min, seqLength);
      int rmax = realValue(max, seqLength);

      if(rmin > rmax){
        overlap = true;
        rmax+= seqLength;
      }

      if(!overlap){
        RangeLocation loc = new RangeLocation(rmin,rmax);
        return new CircularLocation(loc,seqLength);
      }else{
        Location locA;
        Location locB;
        Location compound;

        if(rmin == seqLength){
          locB = new PointLocation(rmin);
        }else{
          locB = new RangeLocation(rmin, seqLength);
        }
        if(rmax - seqLength == 1){
          locA = new PointLocation(1);
        }else{
          locA = new RangeLocation(1,rmax - seqLength);
        }
        compound = LocationTools.union(locA,locB);
        return new CircularLocation(compound,seqLength, rmin);
      }
  }

    static Location intersection(Location locA, Location locB) {
        int circularityA, circularityB;
        Location rawA, rawB;
        if (locA instanceof CircularLocation) {
            circularityA = ((CircularLocation) locA).getLength();
            rawA = ((CircularLocation) locA).getWrapped();
        } else if (locA == Location.empty) {
            return Location.empty;
        } else {
            throw new BioError("Assertion failure: not circular");
        }

        if (locA instanceof CircularLocation) {
            circularityB = ((CircularLocation) locB).getLength();
            rawB = ((CircularLocation) locB).getWrapped();
        } else if (locB == Location.empty) {
            return Location.empty;
        } else {
            throw new BioError("Assertion failure: not circular");
        }

        if (circularityA != circularityB) {
            throw new BioRuntimeException("Can't find intersection of locations on circular sequences of non-equal length");
        }

        Location intersect = LocationTools.intersection(rawA, rawB);
        if (intersect != Location.empty) {
            intersect = new CircularLocation(intersect, circularityA);
        }
        return intersect;
    }


    /**
     * Creates a union of two <code>Locations</code> in a Circular context. This
     * is effectively a join operation. If <code>locA</code> and <code>locB</code> overlap they will be merged.
     * <p> The order in which <code>Location</code>s are joined is important in a
     * circular context. It will be assumed that the 5' end of the resulting join
     * will be the 5' end of <code>locA</code>
     * <p> If the union event merges two locations such that it is not sensible to
     * maintain the original 5' end of <code>locA</code> the 5' end will become the
     * minimum value of the resulting merged Location eg the value returned by <code>getMin()</code>.
     *
     * @param locA the first Location (also determines the 5' End)
     * @param locB the second Location
     * @return the joined location.
     */
  protected static CircularLocation union(CircularLocation locA, CircularLocation locB){
    int length = locA.getLength();
    int fivePrimeEnd = locA.get5PrimeEnd();
    Location temp;


    if(
      locA.isContiguous() &&
      locB.isContiguous() &&
      locA.overlaps(locB)
    ) {
      // the simple case
      try {
        temp = MergeLocation.mergeLocations(locA,locB);
      }
      catch (BioException ex) {
        //this shouldn't happen as conditions have been checked above
        throw new BioError("Assertion Error, cannot build MergeLocation", ex);
      }

      try{
        return new CircularLocation(temp, length, fivePrimeEnd);
      }catch(IllegalArgumentException ae){
        return new CircularLocation(temp, length);
      }

    } else {
      // either may be compound. They may not overlap. We must build the
      // complete list of blocks, merge overlapping blocks and then create the
      // appropriate implementation of Location for the resulting list.

      // list of all blocks
      List locList = new ArrayList();

      // add all blocks in locA
      for(Iterator i = locA.blockIterator(); i.hasNext(); ) {
        locList.add(i.next());
      }

      // add all blocks in locB
      for(Iterator i = locB.blockIterator(); i.hasNext(); ) {
        locList.add(i.next());
      }

      temp = LocationTools._union(locList);
      try{
        return new CircularLocation(temp, length, fivePrimeEnd);
      }catch(IllegalArgumentException iae){
        return new CircularLocation(temp, length);
      }
    }
  }

  /**
   * Tests a location to see if it overlaps the origin of the circular
   * molecule to which it belongs. The origin is taken to be position
   * one.
   *
   * @return true if it does overlap or if it contains the entire sequence
   */
  protected static boolean overlapsOrigin(CircularLocation loc){
    if(loc.getWrapped() instanceof CompoundLocation){
      boolean min = false;
      boolean max = false;
      for(Iterator i = loc.blockIterator(); i.hasNext();){
        Location l = (Location)(i.next());
        if(l.getMax() == loc.getLength()) max = true;
        if(l.getMin() == 1 ) min = true;
      }
      if(min && max){
        return true;
      }else{
        return false;
      }
    }
    else if(loc.getMin() ==1 && loc.getMax() == loc.getLength())return true;
    else return false;
  }


  /**
   * Tests if the specified location is an instance of CircularLocation.
   *
   * @param loc  Locatoin to test
   * @return true if it is a circular location
   */
  protected static boolean isCircular(Location loc){
        boolean toReturn = false;
    try
    {
            if(loc.getDecorator(Class.forName("org.biojava.bio.symbol.CircularLocation")) != null)
            {
                toReturn = true;
            }
        }
        catch(Exception e)
        {
                throw new org.biojava.bio.BioError("class org.biojava.bio.symbol.BetweenLocation could not be loaded");
        }
    return toReturn;
  }


  public static boolean safeOperation(Location locA, Location locB){
    if(isCircular(locA) && isCircular(locB)){
      if (((CircularLocation)locA).getLength() ==
            ((CircularLocation)locB).getLength()) {
        return true;
      }
      else {
        return false;
      }

    }else{
      return false;
    }
  }
}

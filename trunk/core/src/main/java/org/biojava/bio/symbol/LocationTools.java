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
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;

/**
 * Tools class containing a number of operators for working with <code>Location</code> objects.
 *
 * <p>
 * Most of the methods in this class are simple set-wise binary operators: for example,
 * calculate the intersection of two locations.
 * </p>
 *
 * This class provides helpful methods for set-wise manipulation of Location objects.
 *
 * @author Matthew Pocock
 * @author Greg Cox
 * @author Thomas Down
 * @author Mark Schreiber
 * @author Francois Pepin
 * @since 1.2
 */
final public class LocationTools {
    /**
     * Nobody needs one of these.
     */

    private LocationTools() {
    }

  /**
   * Return the union of two locations.
   *
   * <p>
   * The union will be a Location instance that contains every index contained
   * by either locA or locB.
   * </p>
   *
   * @param locA  the first Location
   * @param locB  the second Location
   * @return a Location that is the union of locA and locB
   */
  public static Location union(Location locA, Location locB) {
        if(isDecorated(locA) || isDecorated(locB))
        {
          handleDecorations(locA, locB);
          if(locA instanceof CircularLocation && locB instanceof CircularLocation){
            return CircularLocationTools.union((CircularLocation)locA,(CircularLocation)locB);
          }
          if(BetweenLocationTools.isBetween(locA) || (BetweenLocationTools.isBetween(locB)))
            {
               return BetweenLocationTools.union(locA, locB);
            }
        }

    if(
      locA.isContiguous() &&
      locB.isContiguous() &&
      locA.overlaps(locB)
    ) {
      // the simple case
      Location mloc = null;
      try {
        mloc = MergeLocation.mergeLocations(locA, locB);
      }
      catch (BioException ex) {
        //this shouldn't happen as conditions have been checked above
        throw new BioError("Assertion Error, cannot build MergeLocation",ex);
      }
      return mloc;

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

      return _union(locList);
    }
  }

  /**
   * Return the intersection of two locations.
   * <p>
   * The intersection will be a Location instance that contains every index
   * contained by both locA and locB.
   * </p>
   *
   * @param locA  the first Location
   * @param locB  the second Location
   * @return a Location that is the intersection of locA and locB
   */
  public static Location intersection(Location locA, Location locB) {

    if(isDecorated(locA) || isDecorated(locB))
    {
        handleDecorations(locA, locB);
        if(BetweenLocationTools.isBetween(locA) || (BetweenLocationTools.isBetween(locB)))
            {
                return BetweenLocationTools.intersection(locA, locB);
            }
        if (CircularLocationTools.isCircular(locA) || CircularLocationTools.isCircular(locB)) {
            return CircularLocationTools.intersection(locA, locB);
        }
    }

    if(locA.isContiguous() && locB.isContiguous()) {
        // handle easy case of solid locations
        // Ought to make this bit more efficient --THOMASD
        if (LocationTools.contains(locA, locB)) {
            return locB;
        } else if (LocationTools.contains(locB, locA)) {
            return locA;
        } if(LocationTools.overlaps(locA, locB)) {
            int min = Math.max(locA.getMin(), locB.getMin());
            int max = Math.min(locA.getMax(), locB.getMax());
            return makeLocation(min, max);
        } else {
            return Location.empty;
        }
    } else {
        // One or other of the locations is compound. Build a list of all
        // locations created by finding intersection of all pairwise combinations
        // of blocks in locA and locB. Ignore all Location.empty. Create the
        // appropriate Location instance.
        List locList = new ArrayList();

        List blA = getBlockList(locA);
        {
            int minBlock = blockContainingOrFollowingPoint(blA, locB.getMin());
            int maxBlock = blockContainingOrPreceedingPoint(blA, locB.getMax());
            if (minBlock == -1 || maxBlock == -1) {
                return Location.empty;
            }
            blA = blA.subList(minBlock, maxBlock + 1);
        }
        List blB = getBlockList(locB);
        {
            int minBlock = blockContainingOrFollowingPoint(blB, locA.getMin());
            int maxBlock = blockContainingOrPreceedingPoint(blB, locA.getMax());
            if (minBlock == -1 || maxBlock == -1) {
                return Location.empty;
            }
            blB = blB.subList(minBlock, maxBlock + 1);
        }

        if (blA.size() > blB.size()) {
            List temp = blA;
            blA = blB;
            blB = temp;
        }

        for (Iterator i = blA.iterator(); i.hasNext(); ) {
            Location blocA = (Location) i.next();
            int minHitIndex = blockContainingOrFollowingPoint(blB, blocA.getMin());
            int maxHitIndex = blockContainingOrPreceedingPoint(blB, blocA.getMax());
            for (int hitIndex = minHitIndex; hitIndex <= maxHitIndex; ++hitIndex) {
                Location blocB = (Location) blB.get(hitIndex);
                locList.add(LocationTools.intersection(blocA, blocB));
            }
        }

        return buildLoc(locList);
    }
  }

    /**
     * Return an ordered list of non-overlapping blocks of a location
     * (singleton list for contiguous locations, empty list for empty locations)
     *
     * Note that from this point of view, contiguous circular locations
     * aren't necessarily contiguous :-(.
     *
     * @param l the <code>Location</code> to get the blocks from
     * @return and ordered <code>List</code> of <code>Locations</code>
     */

    private static List getBlockList(Location l) {
        if (l == Location.empty) {
            return Collections.EMPTY_LIST;
        } else if (l instanceof CompoundLocation) {
            return ((CompoundLocation )l).getBlockList();
        } else if (l.isContiguous() && !CircularLocationTools.isCircular(l)) {
            return Collections.nCopies(1, l);
        } else {
            List bl = new ArrayList();
            for (Iterator bi = l.blockIterator(); bi.hasNext(); ) {
                bl.add(bi.next());
            }
            Collections.sort(bl, Location.naturalOrder);
            return bl;
        }
    }

    /**
     * Return the index into <code>bl</code> of a block containing
     * <code>point</code>, or -1 if it is not found.
     */

    private static int blockContainingPoint(List bl, int point) {
        int start = 0;
        int end = bl.size() - 1;

        while (start <= end) {
            int mid = (start + end) / 2;
            Location block = (Location) bl.get(mid);
            if (block.contains(point)) {
                return mid;
            } else if (point < block.getMin()) {
                end = mid - 1;
            } else if (point > block.getMax()) {
                start = mid + 1;
            } else {
                throw new BioError("Assertion failed: overrun in binary search");
            }
        }
        return -1;
    }

    /**
     * Return the index of the block containing <code>point</code>, if it exists,
     * or the block preceeding it.
     */

    private static int blockContainingOrPreceedingPoint(List bl, int point) {
        // System.err.println("COPP");
        int start = 0;
        int length = bl.size();
        int end = length - 1;

        while (start <= end) {
            int mid = (start + end) / 2;
            // System.err.println("Start=" + start + " mid=" + mid + " end=" + end);
            Location block = (Location) bl.get(mid);
            if (block.contains(point)) {
                return mid;
            } else if (point < block.getMin()) {
                end = mid - 1;
            } else if (point > block.getMax()) {
                start = mid + 1;
            } else {
                throw new BioError("Assertion failed: overrun in binary search");
            }
        }

        if (end < length) {
            return end;
        } else {
            return -1;
        }
    }

    /**
     * Return the index of the block containing <code>point</code>, if it exists,
     * or the block preceeding it.
     */

    private static int blockContainingOrFollowingPoint(List bl, int point) {
        // System.err.println("COFP");
        int start = 0;
        int length = bl.size();
        int end = length - 1;

        while (start <= end) {
            int mid = (start + end) / 2;
            // System.err.println("Start=" + start + " mid=" + mid + " end=" + end);
            Location block = (Location) bl.get(mid);
            if (block.contains(point)) {
                return mid;
            } else if (point < block.getMin()) {
                end = mid - 1;
            } else if (point > block.getMax()) {
                start = mid + 1;
            } else {
                throw new BioError("Assertion failed: overrun in binary search");
            }
        }
        if (start >= 0) {
            return start;
        } else {
            return -1;
        }
    }

    /**
     * Determines whether the locations are touching or not (if they could be
     * merged in a single Location.
     * <p>
     * Two locations can merge if they contain at least one index of one
     * beside one index of the other.
     * </p>
     *
     * @param locA  the first Location
     * @param locB  the second Location
     * @return <code>true</code> if they can merge, <code>false</code> otherwise
     */
  public static boolean canMerge(Location locA, Location locB){
    if (overlaps(locA, locB)||
        overlaps(locA.translate(1), locB)||
        overlaps(locA.translate(-1),locB))
      return true;
    return false;
  }
  

  
    /**
     * Determines whether the locations overlap or not.
     * <p>
     * Two locations overlap if they contain at least one index in common.
     * </p>
     *
     * @param locA  the first Location
     * @param locB  the second Location
     * @return <code>true</code> if they overlap, <code>false</code> otherwise
     */

    public static boolean overlaps(Location locA, Location locB) {
        if(isDecorated(locA) || isDecorated(locB))
        {
            handleDecorations(locA, locB);
            if(BetweenLocationTools.isBetween(locA) || (BetweenLocationTools.isBetween(locB)))
                {
                    return BetweenLocationTools.overlaps(locA, locB);
                }
        }

        if(locA.isContiguous() && locB.isContiguous()) {
            // if they are both solid, return whether the extents overlap
            return !(
                     (locA.getMax() < locB.getMin()) ||
                     (locA.getMin() > locB.getMax())
                     );
        } else {
            // System.err.println("Doing complex overlap stuff");

            List blA = getBlockList(locA);
            {
                // System.err.println("In A restriction");
                int minBlock = blockContainingOrFollowingPoint(blA, locB.getMin());
                int maxBlock = blockContainingOrPreceedingPoint(blA, locB.getMax());
                if (minBlock == -1 || maxBlock == -1) {
                    // System.err.println("blA empty: minBlock=" + minBlock +", maxBlock=" + maxBlock);
                    return false;
                }
                blA = blA.subList(minBlock, maxBlock + 1);
            }
            List blB = getBlockList(locB);
            {
                // System.err.println("In B restriction");
                int minBlock = blockContainingOrFollowingPoint(blB, locA.getMin());
                int maxBlock = blockContainingOrPreceedingPoint(blB, locA.getMax());
                if (minBlock == -1 || maxBlock == -1) {
                    // System.err.println("blB empty: minBlock=" + minBlock +", maxBlock=" + maxBlock);
                    return false;
                }
                blB = blB.subList(minBlock, maxBlock + 1);
            }

            // System.err.println("Restricted lists");

            for(Iterator aI = blA.iterator(); aI.hasNext(); ) {
                Location a = (Location) aI.next();
                for(Iterator bI = blB.iterator(); bI.hasNext(); ) {
                    Location b = (Location) bI.next();
                    if(LocationTools.overlaps(a, b)) {
                        return true;
                    }
                }
            }

            return false;
        }
    }

    /**
     * Return <code>true</code> iff all indices in <code>locB</code> are also contained
     * by <code>locA</code>.
     *
     * @param locA The containing location
     * @param locB The contained location
     * @return <code>true</code> is locA contains locB
     */

    public static boolean contains(Location locA, Location locB) {
        if(isDecorated(locA) || isDecorated(locB))
        {
            handleDecorations(locA, locB);
            if(BetweenLocationTools.isBetween(locA) || (BetweenLocationTools.isBetween(locB)))
                {
                    return BetweenLocationTools.contains(locA, locB);
                }
        }

        if (locA.getMin() <= locB.getMin() && locA.getMax() >= locB.getMax()) {
            if (locA.isContiguous()) {
                return true;
            } else {
                List blA = getBlockList(locA);
                for (Iterator biB = locB.blockIterator(); biB.hasNext(); ) {
                    Location bloc = (Location) biB.next();
                    int hitIndex = blockContainingPoint(blA, bloc.getMin());
                    if (hitIndex < 0) {
                        return false;
                    } else {
                        Location hitBloc = (Location) blA.get(hitIndex);
                        if (bloc.getMax() > hitBloc.getMax()) {
                            return false;
                        }
                    }
                }
                return true;
            }
        } else {
            return false;
        }
    }

  /**
   * Return whether two locations are equal.
   * <p>
   * They are equal if both a contains b and b contains a. Equivalently, they
   * are equal if for every point p, locA.contains(p) == locB.contains(p).
   * </p>
   *
   * @param locA the first Location
   * @param locB the second Location
   * @return true if they are equivalent, false otherwise
   */
  public static boolean areEqual(Location locA, Location locB) {
        if(isDecorated(locA) || isDecorated(locB))
        {
                handleDecorations(locA, locB);
                if(BetweenLocationTools.isBetween(locA) || (BetweenLocationTools.isBetween(locB)))
                {
                        return BetweenLocationTools.areEqual(locA, locB);
                }
        }
    // simple check - if one is broken and the other isn't, they aren't equal.
    if(locA.isContiguous() != locB.isContiguous()) {
      return false;
    }

    // both contiguous if one is - check extent only
    if(locA.isContiguous()) {
      return
        (locA.getMin() == locB.getMin()) &&
        (locA.getMax() == locB.getMax());
    }

    // ok - both compound. The blocks returned from blockIterator should each be
    // equivalent.
    Iterator i1 = locA.blockIterator();
    Iterator i2 = locB.blockIterator();

    // while there are more pairs to check...
    while(i1.hasNext() && i2.hasNext()) {
      // check that this pair is equivalent
      Location l1 = (Location) i1.next();
      Location l2 = (Location) i2.next();

      if(
        (l1.getMin() != l2.getMin()) ||
        (l1.getMax() != l2.getMax())
      ) {
        // not equivalent blocks so not equal
        return false;
      }
    }

    // One of the locations had more blocks than the other
    if(i1.hasNext() || i2.hasNext()) {
      return false;
    }

    // Same number of blocks, all equivalent. Must be equal.
    return true;
  }

  /**
   * Create a Location instance from the list of contiguous locations in
   * locList.
   * <p>
   * If the list is empty then Location.empty will be produced. If it is just
   * one element long, then this will be returned as-is. If there are multiple
   * locations then they will be sorted and then used to construct a
   * CompoundLocation.
   *
   * @param locList a List<Location>, where each element is a contiguous location.
   * @return a new Location instance
   */
  static Location buildLoc(List locList) {
    if(locList.size() == 0) {
      return Location.empty;
    } else if(locList.size() == 1) {
      return (Location) locList.get(0);
    } else {
      Collections.sort(locList, Location.naturalOrder);
      return new CompoundLocation(locList);
    }
  }
  
  /**
   * Create a Location instance from the sorted list of contiguous locations in
   * locList.
   * <p>
   * If the list is empty then Location.empty will be produced. If it is just
   * one element long, then this will be returned as-is. If there are multiple
   * locations then they will be sorted and then used to construct a
   * CompoundLocation.
   *
   * @param locList a List<Location>, where each element is a contiguous location.
   * @return a new Location instance
   */
  
  static Location buildLocSorted(List locList) {
      if(locList.size() == 0) {
        return Location.empty;
      } else if(locList.size() == 1) {
        return (Location) locList.get(0);
      } else {
        return new CompoundLocation(locList);
      }
    }

    /**
     * The n-way union of a Collection of locations.  Returns a Location
     * which covers every point covered by at least one of the locations
     * in <code>locs</code>
     *
     * @param locs A collection of locations.
     * @return A union location
     * @throws ClassCastException if the collection contains non-Location objects.
     */

    public static Location union(Collection locs) {
        boolean circular = false;
        List locList = new ArrayList();
        for (Iterator li = locs.iterator(); li.hasNext(); ) {
            Location loc = (Location) li.next();
            if((loc instanceof CircularLocation)){
              circular = true;
              locList.add(loc);
            }else{
              for (Iterator bi = loc.blockIterator(); bi.hasNext(); ) {
                locList.add(bi.next());
              }
            }
        }
        if(circular){
          //need to add these one at a time
          ListIterator li = locList.listIterator();
          CircularLocation loc = (CircularLocation)li.next();
          while(li.hasNext()){
            loc = (CircularLocation)union(loc, (Location)li.next());
          }
          return loc;
        }else{
          return _union(locList);
        }
    }


    static Location _union(List locList) {
      // sort these blocks
      Collections.sort(locList, Location.naturalOrder);

      // merge into this list...
      List joinList = new ArrayList();

      // start iterating over sorted list.
      // last is used as loop variable. We must be careful about zero lengthed
      // lists and also careful to merge overlaps before adding to joinList.
      Iterator i = locList.iterator();
      Location last = Location.empty;


      // prime last
      if(i.hasNext()) {
        last = (Location) i.next();
      }

      // merge or add last with next location
      while(i.hasNext()) {
        Location cur = (Location) i.next();
        if(canMerge(last,cur)) {
          try {
            last = MergeLocation.mergeLocations(last,cur);
          }
          catch (BioException ex) {
            throw new BioError("Cannot make MergeLocation",ex);
          }
        } else {
          joinList.add(last);
          last = cur;
        }
      }

      // handle the end of the loop
      if(last == Location.empty) {
        return Location.empty;
      } else {
        joinList.add(last);
      }

      // now make the appropriate Location instance
      return buildLocSorted(joinList);
    }

  /**
   * Return a contiguous Location from min to max.
   * <p>
   * If min == max then a PointLocation will be made, otherwise, a RangeLocation
   * will be returned.
   *
   * @param min  the Location min value
   * @param max  the Location max value
   * @return a new Location from min to max
   */
  public static Location makeLocation(int min, int max) {
    if(min == max) {
      return new PointLocation(min);
    } else {
      return new RangeLocation(min, max);
    }
  }

    /**
     * A simple method to generate a RangeLocation wrapped
     * in a CircularLocation. The method will cope with situtations
     * where the min is greater than the max. Either of min or max can
     * be negative, or greater than the underlying sequence length. If min and
     * max are equal a wrapped point location will be made.
     *
     * @param min the "left" end of the location
     * @param max the "right" end of the location
     * @param seqLength the lenght of the sequence that the location will
     * be applied to (for purposes of determining origin).
     * @return the new <code>CircularLocation</code>
     */
    public static CircularLocation makeCircularLocation(int min, int max, int seqLength){
        return CircularLocationTools.makeCircLoc(min,max,seqLength);
    }

  /**
   * Checks if the location has a decorator.
   *
   * @param theLocation The location to test for decorators
   * @return True if the location has a decorator and false otherwise
   */
  static boolean isDecorated(Location theLocation)
  {
        return (theLocation instanceof AbstractLocationDecorator);
  }

  /**
   * Currently CircularLocations are handled if and only if
   * CircularLocationTools.safeOperation returns true. For this to be true
   * the locations must have the same sequence length.
   *
   * BetweenLocations are handled in all cases.  Overlap cases, such as
   * between, circular locations have indeterminate behavior.
   *
   * The default behavior is to not handle decorators.  Any new decorators
   * will have to re-visit this method
   */
  private static void handleDecorations(Location locA, Location locB){
    if(CircularLocationTools.isCircular(locA)|| CircularLocationTools.isCircular(locB)){
        if(CircularLocationTools.safeOperation(locA, locB) == false){
              throw new UnsupportedOperationException(
                "Binary operations between Circular and"+
                " Non-Circular locations, or CircularLocations"+
                " from sequences of differing length are currently unsupported.");
        }
      }
    else if(BetweenLocationTools.isBetween(locA) || (BetweenLocationTools.isBetween(locB)))
    {
        // Intentionally blank
    }
    else{
        throw new ClassCastException("Decorated locations are not handled in this version: " + locA.getClass().getName() + ", " + locB.getClass().getName());
      }
  }

  /**
   * Flips a location relative to a length.
   *
   * <p>It is very common in biological sequences to represent locations on a sequence and then reverse that
   * sequence. This method allows locations in the original coordinate space to be transformed int
   * locations in the reverse one.</p>
   *
   * @param loc  the Location to flip
   * @param len  the length of the region to flip within
   * @return  a flipped view of the location
   */
  public static Location flip(Location loc, int len) {
      if(loc instanceof PointLocation) {
          return new PointLocation(len - loc.getMin() + 1);
      } else if(loc instanceof RangeLocation) {
          return new RangeLocation(
            len - loc.getMax() + 1,
            len - loc.getMin() + 1
          );
      } else {
          Iterator bi = loc.blockIterator();
          Location res = (Location) bi.next();
          while(bi.hasNext()) {
              res = LocationTools.union(res, (Location) bi.next());
          }
          return res;
      }
  }

    /**
     * Subtract one location from another.  This methods calculates the set of
     * points which are contains in location <code>keep</code> but not in
     * <code>remove</code>.
     *
     * @param keep A location
     * @param remove A location
     * @return a location containing all points which are in x but not y
     * @since 1.3
     */

    public static Location subtract(Location keep, Location remove) {
        // fixme: documentation hack - I should re-factor this method
        Location x = keep;
        Location y = remove;
        
        if (isDecorated(x) || isDecorated(y)) {
            handleDecorations(x, y);
        }

        List spans = new ArrayList();
        for (Iterator i = x.blockIterator(); i.hasNext(); ) {
            Location xb = (Location) i.next();
            Location yb = LocationTools.intersection(xb, y);
            int pos = xb.getMin();
            for (Iterator j = yb.blockIterator(); j.hasNext(); ) {
                Location sb = (Location) j.next();
                if (sb.getMin() > pos) {
                    spans.add(new RangeLocation(pos, sb.getMin() - 1));
                }
                pos = sb.getMax() + 1;
            }
            if (pos <= xb.getMax()) {
                spans.add(new RangeLocation(pos, xb.getMax()));
            }
        }
        return LocationTools.union(spans);
    }
    
    /**
     * Return the number of positions which are covered by a <code>Location</code>
     *
     * @param loc A location
     * @return the number of distinct points contained by that location
     * @since 1.4
     */
    
    public static int coverage(Location loc) {
        int cov = 0;
        for (Iterator i = loc.blockIterator(); i.hasNext(); ) {
            Location bloc = (Location) i.next();
            cov += (bloc.getMax() - bloc.getMin() + 1);
        }
        return cov;
    }
    
    /**
     * Return a contiguous location running from the minimum to the maximum points of
     * the specified location.
     * 
     * @param loc a location
     * @return a corresponding contiguous location
     */
    
    public static Location shadow(Location loc) {
        if (loc.isContiguous()) {
            return loc;
        } else {
            return new RangeLocation(loc.getMin(), loc.getMax());
        }
    }
    
    /**
     * Return the number of contiguous blocks in a location.
     * 
     * @param loc a location
     * @return the number of blocks
     * @since 1.4
     */
    
    public static int blockCount(Location loc) {
        if (loc.isContiguous()) {
            if (loc instanceof EmptyLocation) {
                return 0;
            } else {
                return 1;
            }
        } else {
            int cnt = 0;
            for (Iterator bi = loc.blockIterator(); bi.hasNext(); ) {
                bi.next(); ++cnt;
            }
            return cnt;
        }
    }
}

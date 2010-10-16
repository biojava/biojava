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

import java.util.Iterator;

/**
 * A set of integers, often used to represent positions on biological sequences.
 *
 * <p>
 * The location will contain some indices between getMin and getMax inclusive.
 * It is not required to contain all indices within this range. It is meant
 * to contain the indices returned by the getMin or getMax. In the event that
 * an operation would produce an
 * invalid or nonsensical range, <code>Location.empty</code> should be returned.
 * </p>
 *
 * <p>
 * Location objects are <strong>always</strong> immutable.
 * </p>
 *
 * <h2>Working with locations</h2>
 *
 * <p>
 * Locations can be constructed in a number of ways:
 * </p>
 *
 * <pre>
 * Location l1 = LocationTools.makeLocation(10, 20);  // Makes a RangeLocation
 * Location l2 = LocationTools.makeLocation(25, 25);  // Makes a PointLocation
 * Location l3 = LocationTools.union(l1, l2); // Construct a non-contiguous
 *                                            // location containing the
 *                                            // points from l1 and l2
 * </pre>
 *
 * @author Matthew Pocock
 * @author Thomas Down
 */

public interface Location {
  /**
   * Create a new instance of Location with all of the same decorators as this
   * instance but with the data stored in <code>loc</code>.
   * <p>
   * The default behavior is to return <loc>loc</loc> unchanged. If the class is
   * a location decorator then it should instantiate an instance of the same
   * type that decorates <code>loc</code>.
   *
   * @param loc  the Location to use as template
   * @return a Location instance based on loc with the same decorators as the
   *         current instance
   */
  Location newInstance(Location loc);

  /**
   * Checks the decorator chain for an instance of <class>decoratorClass</class>
   * and return it if found.
   * <p>
   * The default behavior is to return null. If the current object is a
   * decorator and is an instance of <class>decoratorClass</class> it should
   * return itself. Otherwise, the decorator should chain this method onto the
   * instance it wraps.
   *
   * @param decoratorClass  the Class to check
   * @return a Location if an instance of this class is present in the decorator
   *         chain and null otherwise.
   */
  Location getDecorator(Class decoratorClass);
  /**
   * The minimum position contained.
   * <p>
   * <b>WARNING:</b> The location will <b>not</b> contain every point between <code>getMin()</code>
   * and <code>getMax()</code> if <code>isContiguous()</code> is false. If <code>isContiguous()</code>
   * does return false you should use the <code>Iterator</code> returned by <code>blockIterator()</code>
   * to iterate over the minimum set of contiguous blocks that make up this <code>Location</code>
   *
   * @return	the minimum position contained
   */
  int getMin();
  /**
   * The maximum position contained.
   * <p>
   * <b>WARNING:</b> The location will <b>not</b> contain every point between <code>getMin()</code>
   * and <code>getMax()</code> if <code>isContiguous()</code> is false. If <code>isContiguous()</code>
   * does return false you should use the <code>Iterator</code> returned by <code>blockIterator()</code>
   * to iterate over the minimum set of contiguous blocks that make up this <code>Location</code>
   *
   * @return	the maximum position contained
   */
  int getMax();

  /**
   * Checks if these two locations overlap, using this location's
   * concept of overlapping.
   * <p>
   * Abstractly, two locations overlap if they both contain any point.
   *
   * @param l	the Location to check
   * @return	true if they overlap, otherwise false
   *
   */
  boolean overlaps(Location l);
  /**
   * Checks if this location contains the other.
   * <p>
   * Abstractly, a location contains another if every point in the
   * other location is contained within this one.
   *
   * @param l	the Location to check
   * @return	true if this contains l, otherwise false
   *
   */
  boolean contains(Location l);
  /**
   * Checks if this location contains a point.
   *
   * @param p	the point to check
   * @return	true if this contains p, otherwise false
   */
  boolean contains(int p);

  /**
   * Checks if this location is equivalent to the other.
   * <p>
   * Abstractly, a location is equal to another if for every point in one
   * it is also in the other. This is equivalent to
   * a.contains(b) && b.contains(a). You should call LocationTools.areEqual
   * after casting l to Location.
   *
   * @param l	the Object to check
   * @return	true if this equals l, otherwise false
   */
  boolean equals(Object l);

  /**
   * Returns a Location that contains all points common to both ranges.
   *
   * @param l	the Location to intersect with
   * @return	a Location containing all points common to both, or
   *              the empty range if there are no such points
   *
   */
  Location intersection(Location l);
  /**
   * Return a Location containing all points in either ranges.
   *
   * @param l	the Location to union with
   * @return	a Location representing the union
   *
   */
  Location union(Location l);

  /**
   * Return the symbols in a sequence that fall within this range.
   *
   * @param seq	the SymbolList to process
   * @return	the SymbolList containing the symbols in seq in this range
   *
   */
  SymbolList symbols(SymbolList seq);

  /**
   * Create a location that is a translation of this location.
   *
   * @param dist  the distance to translate (to the right)
   */
  Location translate(int dist);

  /**
   * Determine if a Location is contiguous.
   *
   * @return <code>true</code> if and only if this Location
   *         contains every point from <code>min</code> to
   *         <code>max</code> inclusive.
   */
  boolean isContiguous();

  /**
   * Return an Iterator over the set of maximal contiguous sub-locations.
   * <p>
   * Given any location, it can be considered to contain zero or more
   * maximal contiguous blocks of width 1 or greater. The empty location is
   * composed from nothing. A contiguous location is composed from itself.
   * A non-contiguous location is composed from contiguous blocks seperated by
   * gaps.
   * <p>
   * This method should return an Iterator over these maximally contiguous blocks
   * starting with the left-most block, and finishing at the right-most block.
   *
   * @return an Iterator over Location objects that are the maximally contiguous
   *         set of locations contained within this location
   */
  Iterator<Location> blockIterator();

  /**
   * The <code>Location</code> which contains no points.
   *
   * <p>
   * This object contains nothing. Its minimum value is Integer.MAX_VALUE.
   * Its maximum value is Integer.MIN_VALUE. It overlaps nothing. It is
   * equal to nothing. Intersection results in the empty range. Union
   * results in the argument range. Symbols returns an empty array.
   * <p>
   * Every day, in every way, empty becomes more and more boring.
   */
  public static final Location empty = new EmptyLocation();
  
  /**
   * The <code>Location</code> which contains all points.
   *
   * <p>
   * This object contains every point. It's minimum value is Integer.MIN_VALUE,
   * and it's maximum value is Integer.MAX_VALUE. It overlaps and contains
   * everything.
   * </p>
   */
  public static final Location full = new RangeLocation(Integer.MIN_VALUE, Integer.MAX_VALUE);

  /**
   * Comparator which orders Locations naturally.  Locations
   * are sorted primarily on the basis of their <code>getMin()</code>
   * value.  In cases where that is equal, they are secondarily sorted
   * by <code>getMax()</code> value.
   */

  static final LocationComparator naturalOrder = new LocationComparator();
}


package org.biojava3.core.sequence.location.template;

import java.util.List;

import org.biojava3.core.sequence.Strand;
import org.biojava3.core.sequence.location.SimpleLocation;

/**
 * Sets of integers used to represent the location of features on sequence. A
 * location can be a single set of bounds or composed of multiple
 * sub-locations. Each sub-location is a Location and therefore subject to the
 * same rules.
 *
 * Rather than the normal concept of min and max here we define start and
 * stop which denotes nothing about the order of the location bounds. This
 * leans towards the idea of circular locations where min & max still apply
 * they convey incorrect ideas and assumptions about the relationship between
 * these two integer values.
 *
 * @author ayates
 */
public interface Location extends Iterable<Location> {

  /**
   * Basic location which is set to the minimum and maximum bounds of
   * {@link Integer}. {@link Strand} is set to {@link Strand#UNDEFINED}.
   */
  public static final Location EMPTY =
    new SimpleLocation(Integer.MIN_VALUE, Integer.MAX_VALUE, Strand.UNDEFINED);

  /**
   * Start of the location; not necessarily the min position
   */
  int getStart();

  /**
   * Start of the location; not necessarily the max position
   */
  int getEnd();

  /**
   * Strand which the location is located on
   */
  Strand getStrand();

  /**
   * Gives access to the sub locations for this location. However this does
   * not return sub-locations of sub-locations. For that functionality use
   * {@link #getAllSubLocations()}.
   *
   * @return A list of a single level of sub-locations
   */
  List<Location> getSubLocations();

  /**
   * An extension to {@link #getSubLocations()} which returns sub-locations
   * of sub-locations; this will continue until it runs out of those locations.
   *
   * @return List of all sub locations including sub-locations of sub locations
   */
  List<Location> getAllSubLocations();

  /**
   * Returns true if the location is considered to be complex; normally this
   * means the location is actually composed of sub-locations.
   */
  boolean isComplex();

  /**
   * Indicates if this location is circular. We do not capture how many times
   * we are circular just that we are.
   */
  boolean isCircular();
}

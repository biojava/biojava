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
 * Created on 01-21-2010
 */
package org.biojava3.core.sequence.location.template;

import java.util.ArrayList;
import java.util.List;

import org.biojava3.core.sequence.Strand;
import org.biojava3.core.sequence.location.SimpleLocation;
import org.biojava3.core.sequence.location.SimplePoint;
import org.biojava3.core.sequence.template.Accessioned;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.Sequence;

/**
 * Sets of integers used to represent the location of features on sequence. A
 * location can be a single set of bounds or composed of multiple
 * sub-locations. Each sub-location is a Location and therefore subject to the
 * same rules.
 * 
 * @author ayates
 */
public interface Location extends Iterable<Location>, Accessioned {

    /**
     * Basic location which is set to the minimum and maximum bounds of
     * {@link Integer}. {@link Strand} is set to {@link Strand#UNDEFINED}.
     */
    public static final Location EMPTY =
            new SimpleLocation(Integer.MIN_VALUE, Integer.MAX_VALUE, Strand.UNDEFINED);

    /**
     * Start of the location
     */
    Point getStart();

    /**
     * End of the location
     */
    Point getEnd();

    /**
     * Returns the length of the outer bounds of this location
     */
    int getLength();

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
    List<Location> getRelevantSubLocations();

    /**
     * Returns true if the location is considered to be complex meaning
     * the location is actually composed of sub-locations.
     */
    boolean isComplex();

    /**
     * Indicates if this location is circular.
     */
    boolean isCircular();

    /**
     * Returns true if the position is meant to represent a point between
     * two points such as 78^79. Only valid if start and stop are next to
     * each other.
     */
    boolean isBetweenCompounds();

    /**
     * Will return a SequenceReader object which represents the outer bounds
     * of this Location
     *
     * @param &lt;C&gt; The type of compound to use
     * @param sequence The sequence object to work with
     * @return The sequence
     */
    <C extends Compound> Sequence<C> getSubSequence(Sequence<C> sequence);

    /**
     * Will return a SequenceReader object which offers a view of all resolved
     * locations i.e. those locations which are not complex and define the
     * true Sequence represented
     *
     * @param &lt;C&gt; The type of compound to use
     * @param sequence The sequence object to work with
     * @return The full assembled sequence
     */
    <C extends Compound> Sequence<C> getRelevantSubSequence(Sequence<C> sequence);

    /**
     * Helper methods for use with the Location classes. Taking its
     * inspiration from the RichSequence.Tools class from the old BioJava
     */
    public static class Tools {

        /**
         * Used for building a location from a series of sub-locations
         */
        public static Location location(List<Location> locations, Integer sequenceLength, String type) {
            type = (type == null) ? "join" : type;
            sequenceLength = (sequenceLength == null) ? -1 : sequenceLength;



            return null;
        }

        /**
         * Returns a location object which unlike the location constructors
         * allows you to input reverse coordinates and will convert
         * these into the right location on the positive strand.
         */
        public static Location location(int start, int end, Strand strand, int length) {
            int min = Math.min(start, end);
            //if this is true then we have a coord on the +ve strand even though Strand could be negative
            boolean isReverse = (min != start);
            if (isReverse) {
                return new SimpleLocation(
                        new SimplePoint(start).reverse(length),
                        new SimplePoint(end).reverse(length),
                        strand);
            }
            return new SimpleLocation(start, end, strand);
        }

        /**
         * Converts a location which defines the outer bounds of a circular
         * location and splits it into the required portions. Unlike any
         * other location builder this allows you to express your input
         * location on the reverse strand
         *
         * @param location The location which currently expresses the outer
         * bounds of a circular location.
         * @param length The length of the circular genomic unit
         * @return The circular location; can optionally return a normal non
         * circular location if the one you give is within the bounds of
         * the length
         */
        public static Location circularLocation(int start, int end, Strand strand, int length) {

            int min = Math.min(start, end);
            int max = Math.max(start, end);
            //Tells us we're dealing with something that's not _right_
            boolean isReverse = (min != start);

            if (min > length) {
                throw new IllegalArgumentException("Cannot process a "
                        + "location whose lowest coordinate is less than "
                        + "the given length " + length);
            }

            //If max positon was less than length the return a normal location
            if (max <= length) {
                return location(start, end, strand, length);
            }

            //Fine for forward coords (i..e start < end)
            int modStart = modulateCircularIndex(start, length);
            int modEnd = modulateCircularIndex(end, length);
            int numberOfPasses = completeCircularPasses(Math.max(start, end), length);

            if (isReverse) {
                int reversedModStart = new SimplePoint(modStart).reverse(length).getPosition();
                int reversedModEnd = new SimplePoint(modEnd).reverse(length).getPosition();
                modStart = reversedModStart;
                modEnd = reversedModEnd;
                start = reversedModStart;
                //+1 to number of passes to skip the run encoded by the start
                end = (length * (numberOfPasses + 1)) + modEnd;
            }

            List<Location> locations = new ArrayList<Location>();
            locations.add(new SimpleLocation(modStart, length, strand));
            for (int i = 0; i < numberOfPasses; i++) {
                locations.add(new SimpleLocation(1, length, strand));
            }
            locations.add(new SimpleLocation(1, modEnd, strand));
            return new SimpleLocation(new SimplePoint(start),
                    new SimplePoint(end), strand, true, false, locations);
        }

        private static interface LocationPredicate {

            boolean accept(Location previous, Location current);
        }

        /**
         * Scans through a list of locations to find the Location with the
         * lowest start
         */
        public static Location getMin(List<Location> locations) {
            return scanLocations(locations, new LocationPredicate() {
                
                public boolean accept(Location previous, Location current) {
                    int res = current.getStart().compareTo(previous.getStart());
                    return res < 0;
                }
            });
        }

        /**
         * Scans through a list of locations to find the Location with the
         * highest end
         */
        public static Location getMax(List<Location> locations) {
            return scanLocations(locations, new LocationPredicate() {
                
                public boolean accept(Location previous, Location current) {
                    int res = current.getEnd().compareTo(previous.getEnd());
                    return res > 0;
                }
            });
        }

        /**
         * Used for scanning through a list of locations; assumes the
         * locations given will have at least one value otherwise
         * we will get a null pointer
         */
        private static Location scanLocations(List<Location> locations, LocationPredicate predicate) {
            Location location = null;
            for (Location l : locations) {
                if (location == null) {
                    location = l;
                } else {
                    if (predicate.accept(location, l)) {
                        location = l;
                    }
                }
            }
            return location;
        }

        /**
         * Takes a point on a circular location and moves it left until it falls
         * at the earliest possible point that represents the same base.
         *
         * @param index Index of the position to work with
         * @param seqLength Length of the Sequence
         * @return The shifted point
         */
        public static int modulateCircularIndex(int index, int seqLength) {
            // Dummy case
            if (seqLength == 0) {
                return index;
            }
            // Modulate
            while (index > seqLength) {
                index -= seqLength;
            }
            return index;
        }

        /**
         * Works in a similar way to modulateCircularLocation but returns
         * the number of complete passes over a Sequence length a circular
         * location makes i.e. if we have a sequence of length 10 and the
         * location 3..52 we make 4 complete passes through the genome to
         * go from position 3 to position 52.
         */
        public static int completeCircularPasses(int index, int seqLength) {
            int count = 0;
            while (index > seqLength) {
                count++;
                index -= seqLength;
            }
            return count - 1;
        }
    }
}

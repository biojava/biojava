/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava3.core.sequence.location;

import java.util.ArrayList;
import java.util.List;
import org.biojava3.core.sequence.Strand;
import org.biojava3.core.sequence.location.template.Location;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
 /**
     * Helper methods for use with the Location classes. Taking its
     * inspiration from the RichSequence.Tools class from the old BioJava
     */
    public class  LocationHelper {

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
                @Override
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
                @Override
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

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
package org.biojava3.core.sequence.location;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import org.biojava3.core.exceptions.ParserException;
import org.biojava3.core.sequence.AccessionID;
import org.biojava3.core.sequence.Strand;
import org.biojava3.core.sequence.location.template.Location;
import org.biojava3.core.sequence.location.template.Point;

/**
 * Helper methods for use with the Location classes. Taking its
 * inspiration from the RichSequence.Tools class from the old BioJava
 */
public class LocationHelper {

    /**
     * Used as a thin wrapper to the {@link #location(java.util.List, java.lang.String) }
     * method to bring the given location list together as a join (the default
     * type)
     */
    public static Location location(List<Location> subLocations) {
        return location(subLocations, "join");
    }

    /**
     * Builds a location from a List of locations; this can be circular or
     * linear joins. The code expects that these locations are in
     * a sensible format.
     *
     * @param subLocations The list of locations to use to build the location.
     * If given a list of size 1 we will return that location.
     * @param type The type of join for this location; defaults to join
     * @return
     */
    public static Location location(List<Location> subLocations, String type) {
        if (subLocations.size() == 1) {
            return subLocations.get(0);
        }

        boolean circular = detectCicular(subLocations);
        Strand strand = detectStrand(subLocations);
        Point start = detectStart(subLocations);
        Point end = detectEnd(subLocations, circular);
        Location l;
        if ("join".equals(type)) {
            l = new SimpleLocation(start, end, strand, circular, subLocations);
        }
        else if ("order".equals(type)) {
            l = new InsdcLocations.OrderLocation(start, end, strand, circular, subLocations);
        }
        else if ("one-of".equals(type)) {
            l = new InsdcLocations.OneOfLocation(subLocations);
        }
        else if ("group".equals(type)) {
            l = new InsdcLocations.GroupLocation(start, end, strand, circular, subLocations);
        }
        else if ("bond".equals(type)) {
            l = new InsdcLocations.BondLocation(subLocations);
        }
        else {
            throw new ParserException("Unknown join type " + type);
        }

        return l;
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
            }
            else {
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

    /**
     * Loops through the given list of locations and returns true if it looks
     * like they represent a circular location. Detection cannot happen if
     * we do not have consistent accessions
     */
    public static boolean detectCicular(List<Location> subLocations) {
        boolean isCircular = false;
        if(! consistentAccessions(subLocations))
            return isCircular;

        int lastMax = 0;
        for (Location sub : subLocations) {
            if (sub.getEnd().getPosition() > lastMax) {
                lastMax = sub.getEnd().getPosition();
            }
            else {
                isCircular = true;
                break;
            }
        }
        return isCircular;
    }

    /**
     * Scans a list of locations and returns true if all the given locations
     * are linked to the same sequence. A list of null accessioned locations
     * is the same as a list where the accession is the same
     *
     * @param subLocations The locations to scan
     * @return Returns a boolean indicating if this is consistently accessioned
     */
    public static boolean consistentAccessions(List<Location> subLocations) {
        Set<AccessionID> set = new HashSet<AccessionID>();
        for(Location sub: subLocations) {
            set.add(sub.getAccession());
        }
        return set.size() == 1;
    }

    /**
     * Loops through the given list of locations and returns the consensus
     * Strand class. If the class switches then we will return an undefined
     * strand
     */
    public static Strand detectStrand(List<Location> subLocations) {
        Strand strand = subLocations.get(0).getStrand();
        for (Location sub : subLocations) {
            if (strand != sub.getStrand()) {
                strand = Strand.UNDEFINED;
                break;
            }
        }
        return strand;
    }

    /**
     * Assumes that the first element is the start & clones it
     */
    public static Point detectStart(List<Location> subLocations) {
        return subLocations.get(0).getStart().clonePoint();
    }

    /**
     * This will attempt to find what the last point is and returns that
     * position. If the location is circular this will return the total length
     * of the location and does not mean the maximum point on the Sequence
     * we may find the locations on
     */
    public static Point detectEnd(List<Location> subLocations, boolean isCircular) {
        int end = 0;
        Point lastPoint = null;
        if(isCircular) {
            for (Location sub : subLocations) {
                lastPoint = sub.getEnd();
                end += lastPoint.getPosition();
            }
        }
        else {
            lastPoint = subLocations.get(subLocations.size()-1).getEnd();
            end = lastPoint.getPosition();
        }
        return new SimplePoint(end, lastPoint.isUnknown(), lastPoint.isUncertain());
    }
}

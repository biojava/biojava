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

import java.util.Arrays;
import java.util.List;
import org.biojava3.core.sequence.Strand;
import org.biojava3.core.sequence.location.template.AbstractLocation;
import org.biojava3.core.sequence.location.template.Location;
import org.biojava3.core.sequence.location.template.Point;

/**
 * A collection of locations which are used whenever we work with INSDC; some
 * of which could be deprecated (from INSDC's point of view) yet appear
 * in records.
 *
 * @author ayates
 */
public class InsdcLocations {

    /**
     * Used to represent bond locations equivalent to bond(7,8) or bond(7).
     * Bond locations are single position complex locations
     */
    public static class BondLocation extends AbstractLocation {

        public BondLocation(Location... subLocations) {
            this(Arrays.asList(subLocations));
        }

        public BondLocation(List<Location> subLocations) {
            super();
            Location min = Tools.getMin(subLocations);
            Location max = Tools.getMax(subLocations);
            setStart(min.getStart());
            setEnd(max.getEnd());
            setStrand(Strand.UNDEFINED);
            setSubLocations(subLocations);
            assertLocation();
        }

        @Override
        protected final void assertLocation() {
            for (Location l : this) {
                Point start = l.getStart();
                Point end = l.getEnd();
                if (!start.equals(end)) {
                    throw new IllegalStateException("The start "
                            + start + " is not equal to the end "
                            + end + ". bond locations must be a single "
                            + "compound long");
                }
            }
        }
    }

    /**
     * Used to describe a 5' to 3' ordering but no firm assurance it is correct
     */
    public static class OrderLocation extends SimpleLocation {

        public OrderLocation(Point start, Point end, Strand strand,
                boolean circular, Location... subLocations) {
            super(start, end, strand, circular, subLocations);
        }

        public OrderLocation(Point start, Point end, Strand strand,
                Location... subLocations) {
            this(start, end, strand, false, subLocations);
        }

        public OrderLocation(int start, int end, Strand strand,
                Location... subLocations) {
            this(new SimplePoint(start), new SimplePoint(end), strand, false, subLocations);
        }

        public OrderLocation(Point start, Point end, Strand strand,
                boolean circular, List<Location> subLocations) {
            super(start, end, strand, circular, subLocations);
        }

        public OrderLocation(Point start, Point end, Strand strand,
                List<Location> subLocations) {
            this(start, end, strand, false, subLocations);
        }

        public OrderLocation(int start, int end, Strand strand,
                List<Location> subLocations) {
            this(new SimplePoint(start), new SimplePoint(end), strand, false, subLocations);
        }
    }

    /**
     * Deprecated in INSDC yet still appears; equivalent to the order()
     * directive except no 5' to 3' ordering is defined. The constructor
     * reflects this relationship and only allows the creation of complex
     * locations
     */
    public static class GroupLocation extends SimpleLocation {

        public GroupLocation(Point start, Point end, Strand strand,
                boolean circular, Location... subLocations) {
            super(start, end, strand, circular, subLocations);
        }

        public GroupLocation(Point start, Point end, Strand strand,
                Location... subLocations) {
            this(start, end, strand, false, subLocations);
        }

        public GroupLocation(int start, int end, Strand strand,
                Location... subLocations) {
            this(new SimplePoint(start), new SimplePoint(end), strand, false, subLocations);
        }

        public GroupLocation(Point start, Point end, Strand strand,
                boolean circular, List<Location> subLocations) {
            super(start, end, strand, circular, subLocations);
        }

        public GroupLocation(Point start, Point end, Strand strand,
                List<Location> subLocations) {
            this(start, end, strand, false, subLocations);
        }

        public GroupLocation(int start, int end, Strand strand,
                List<Location> subLocations) {
            this(new SimplePoint(start), new SimplePoint(end), strand, false, subLocations);
        }
    }

    /**
     * Deprecated in INSDC; refers to a set of locations of which one
     * location could be valid e.g. one-of(location, location, location).
     * Originally used to describe split locations in alternative splicing
     * or variations. Now these are dealt with in their own feature fields.
     *
     * The default location is chosen to be the first however if you think
     * you need to work with this location you should use the sub-locations.
     */
    public static class OneOfLocation extends AbstractLocation {

        public OneOfLocation(Location... locations) {
            this(Arrays.asList(locations));
        }

        public OneOfLocation(List<Location> locations) {
            super();
            if (locations.isEmpty()) {
                throw new IllegalArgumentException("Need locations to build a OneOfLocation");
            }
            Location l = locations.get(0);
            setStart(l.getStart());
            setEnd(l.getEnd());
            setStrand(l.getStrand());
            setBetweenCompounds(l.isBetweenCompounds());
            setCircular(l.isCircular());
            setSubLocations(locations);
        }
    }
}

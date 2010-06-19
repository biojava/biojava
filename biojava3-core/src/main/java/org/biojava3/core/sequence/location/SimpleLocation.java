package org.biojava3.core.sequence.location;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.biojava3.core.sequence.Strand;
import org.biojava3.core.sequence.location.template.AbstractLocation;
import org.biojava3.core.sequence.location.template.Location;
import org.biojava3.core.sequence.location.template.Point;

/**
 * Very basic implementation of the Location interface which defines a series
 * of simple constructors.
 *
 * @author ayates
 */
public class SimpleLocation extends AbstractLocation {

    private static final List<Location> EMPTY_LOCS = Collections.emptyList();

    public SimpleLocation(int start, int end) {
        this(new SimplePoint(start), new SimplePoint(end));
    }

    public SimpleLocation(Point start, Point end) {
        this(start, end, Strand.POSITIVE);
    }

    public SimpleLocation(int start, int end, Strand strand) {
        this(new SimplePoint(start), new SimplePoint(end), strand);
    }

    public SimpleLocation(Point start, Point end, Strand strand) {
        super(start, end, strand, false, false, EMPTY_LOCS);
    }

    public SimpleLocation(int start, int end, Strand strand, Location... subLocations) {
        this(new SimplePoint(start), new SimplePoint(end), strand, subLocations);
    }

    public SimpleLocation(Point start, Point end, Strand strand, Location... subLocations) {
        super(start, end, strand, false, false, Arrays.asList(subLocations));
    }

    public SimpleLocation(Point start, Point end, Strand strand, boolean circular, Location... subLocations) {
        super(start, end, strand, circular, false, Arrays.asList(subLocations));
    }

    public SimpleLocation(Point start, Point end, Strand strand, boolean circular, List<Location> subLocations) {
        super(start, end, strand, circular, false, subLocations);
    }

    public SimpleLocation(Point start, Point end, Strand strand, boolean circular, boolean betweenBases, List<Location> subLocations) {
        super(start, end, strand, circular, betweenBases, subLocations);
    }
}

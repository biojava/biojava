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
package org.biojava.nbio.core.sequence.location;

import org.biojava.nbio.core.sequence.AccessionID;
import org.biojava.nbio.core.sequence.Strand;
import org.biojava.nbio.core.sequence.location.template.AbstractLocation;
import org.biojava.nbio.core.sequence.location.template.Location;
import org.biojava.nbio.core.sequence.location.template.Point;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Very basic implementation of the Location interface which defines a series
 * of simple constructors.
 *
 * @author ayates
 * @author Paolo Pavan
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

	public SimpleLocation(int start, int end, Strand strand, List<Location> subLocations) {
		this(new SimplePoint(start), new SimplePoint(end), strand, subLocations);
	}

	public SimpleLocation(Point start, Point end, Strand strand) {

		super(start, end, strand, false, false, new ArrayList<Location>());
	}

	public SimpleLocation(Point start, Point end, Strand strand, AccessionID accession) {
		super(start, end, strand, false, false, accession, EMPTY_LOCS);
	}

	public SimpleLocation(Point start, Point end, Strand strand, boolean betweenCompounds, AccessionID accession) {
		super(start, end, strand, false, betweenCompounds, accession, EMPTY_LOCS);
	}

	public SimpleLocation(Point start, Point end, Strand strand, boolean circular, boolean betweenBases) {
		super(start, end, strand, circular, betweenBases, EMPTY_LOCS);
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

	public SimpleLocation(Point start, Point end, Strand strand, List<Location> subLocations) {
		super(start, end, strand, false, false, subLocations);
	}

	public SimpleLocation(Point start, Point end, Strand strand, boolean circular, boolean betweenBases, List<Location> subLocations) {
		super(start, end, strand, circular, betweenBases, subLocations);
	}
}
